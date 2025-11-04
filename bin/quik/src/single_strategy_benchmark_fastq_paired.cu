//
// Created for Nextflow pipeline integration
// Single strategy barcode calling for paired FASTQ files with flexible parameters
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <chrono>
#include <cassert>
#include <map>

#include "read_sequence_file.h"
#include "sequence.h"
#include "barcode_assignment.h"
#include "constants.h"

// Include all possible strategies
#include "k_mer_filtered_calling_gpu_v1.cuh"
#include "k_mer_filtered_calling_gpu_v2.cuh"
#include "k_mer_filtered_calling_gpu_v3.cuh"
#include "k_mer_filtered_calling_gpu_v4.cuh"
#include "k_mer_filtered_calling_gpu_v5.cuh"
#include "k_mer_filtered_calling_host_v1.h"
#include "k_mer_filtered_calling_host_v2.h"
#include "two_step_k_mer_filtered_calling_gpu_v1.cuh"
#include "two_step_k_mer_filtered_calling_host_v1.h"

// Distance measure constants
#define LEVENSHTEIN_DISTANCE 0
#define SEQUENCE_LEVENSHTEIN_DISTANCE 1
#define PSEUDO_DISTANCE 2

// FASTQ record structure
struct fastq_record {
    std::string header;
    std::string sequence;
    std::string plus;
    std::string quality;
};

// Read a single FASTQ record
bool read_fastq_record(std::ifstream& file, fastq_record& record) {
    if (!std::getline(file, record.header)) return false;
    if (!std::getline(file, record.sequence)) return false;
    if (!std::getline(file, record.plus)) return false;
    if (!std::getline(file, record.quality)) return false;
    return true;
}

// Extract barcode from sequence with padding for short sequences
std::string extract_barcode(const std::string& seq_str, int barcode_start, int barcode_length) {
    std::string barcode_str;
    
    int seq_length = seq_str.length();
    
    if (seq_length >= barcode_start + barcode_length) {
        // Normal case: sequence is long enough
        barcode_str = seq_str.substr(barcode_start, barcode_length);
    } else if (seq_length > barcode_start) {
        // Sequence is too short: extract what we can and pad with 'A'
        int available = seq_length - barcode_start;
        barcode_str = seq_str.substr(barcode_start, available);
        barcode_str.append(barcode_length - available, 'A');
    } else {
        // Sequence is shorter than barcode_start: use all 'A's
        barcode_str = std::string(barcode_length, 'A');
    }
    
    return barcode_str;
}

// Pad or truncate sequence to specified length
std::string normalize_sequence(const std::string& seq_str, int target_length) {
    if (seq_str.length() >= target_length) {
        return seq_str.substr(0, target_length);
    } else {
        std::string padded = seq_str;
        padded.append(target_length - seq_str.length(), 'A');
        return padded;
    }
}

// Strategy execution function
barcode_assignment execute_strategy(
    const std::string& strategy,
    const std::vector<sequence>& barcodes,
    const std::vector<sequence>& reads,
    int distance_measure,
    int rejection_threshold
) {
    if (strategy == "4_mer_gpu_v1" || strategy == "4_7_mer_gpu_v1") {
        return k_mer_filtered_calling_gpu_v1(barcodes, reads, 4, distance_measure, rejection_threshold);
    }
    else if (strategy == "5_mer_gpu_v1" || strategy == "5_7_mer_gpu_v1") {
        return k_mer_filtered_calling_gpu_v1(barcodes, reads, 5, distance_measure, rejection_threshold);
    }
    else if (strategy == "6_mer_gpu_v1" || strategy == "6_7_mer_gpu_v1") {
        return k_mer_filtered_calling_gpu_v1(barcodes, reads, 6, distance_measure, rejection_threshold);
    }
    else if (strategy == "7_mer_gpu_v1" || strategy == "7_7_mer_gpu_v1") {
        return k_mer_filtered_calling_gpu_v1(barcodes, reads, 7, distance_measure, rejection_threshold);
    }
    else if (strategy == "4_mer_gpu_v2") {
        return k_mer_filtered_calling_gpu_v2(barcodes, reads, 4, distance_measure, rejection_threshold);
    }
    else if (strategy == "5_mer_gpu_v2") {
        return k_mer_filtered_calling_gpu_v2(barcodes, reads, 5, distance_measure, rejection_threshold);
    }
    else if (strategy == "4_mer_gpu_v3") {
        return k_mer_filtered_calling_gpu_v3(barcodes, reads, 4, distance_measure, rejection_threshold);
    }
    else if (strategy == "5_mer_gpu_v3") {
        return k_mer_filtered_calling_gpu_v3(barcodes, reads, 5, distance_measure, rejection_threshold);
    }
    else if (strategy == "4_mer_gpu_v4") {
        return k_mer_filtered_calling_gpu_v4(barcodes, reads, 4, distance_measure, rejection_threshold);
    }
    else if (strategy == "5_mer_gpu_v4") {
        return k_mer_filtered_calling_gpu_v4(barcodes, reads, 5, distance_measure, rejection_threshold);
    }
    else if (strategy == "6_mer_gpu_v4") {
        return k_mer_filtered_calling_gpu_v4(barcodes, reads, 6, distance_measure, rejection_threshold);
    }
    else if (strategy == "7_mer_gpu_v4") {
        return k_mer_filtered_calling_gpu_v4(barcodes, reads, 7, distance_measure, rejection_threshold);
    }
    else if (strategy == "4_mer_gpu_v5") {
        return k_mer_filtered_calling_gpu_v5(barcodes, reads, 4, distance_measure, rejection_threshold);
    }
    else if (strategy == "5_mer_gpu_v5") {
        return k_mer_filtered_calling_gpu_v5(barcodes, reads, 5, distance_measure, rejection_threshold);
    }
    else if (strategy == "two_step_gpu_v1") {
        // Two-step: use 5-mer for first pass, 7-mer for second pass
        return two_step_k_mer_filtered_calling_gpu_v1(barcodes, reads, 5, 7, distance_measure, rejection_threshold, rejection_threshold);
    }
    else if (strategy == "host_v1") {
        return k_mer_filtered_calling_host_v1(barcodes, reads, 5, distance_measure, rejection_threshold);
    }
    else if (strategy == "host_v2") {
        return k_mer_filtered_calling_host_v2(barcodes, reads, 5, distance_measure, rejection_threshold);
    }
    else if (strategy == "two_step_host_v1") {
        // Two-step: use 5-mer for first pass, 7-mer for second pass
        return two_step_k_mer_filtered_calling_host_v1(barcodes, reads, 5, 7, distance_measure, rejection_threshold, rejection_threshold);
    }
    else {
        std::cerr << "ERROR: Unknown strategy: " << strategy << std::endl;
        std::cerr << "Available strategies: 4_7_mer_gpu_v1, 5_7_mer_gpu_v1, 6_7_mer_gpu_v1, 7_7_mer_gpu_v1, "
                  << "4_mer_gpu_v2, 5_mer_gpu_v2, 4_mer_gpu_v3, 5_mer_gpu_v3, "
                  << "4_mer_gpu_v4, 5_mer_gpu_v4, 6_mer_gpu_v4, 7_mer_gpu_v4, "
                  << "4_mer_gpu_v5, 5_mer_gpu_v5, two_step_gpu_v1, host_v1, host_v2, two_step_host_v1" << std::endl;
        exit(1);
    }
}

int main(int argc, char** argv) {

    if (argc != 11) {
        std::cerr << "Usage: " << argv[0] << " <barcode_file> <r1_fastq> <r2_fastq> "
                  << "<barcode_start> <barcode_length> <strategy> <distance_measure> "
                  << "<rejection_threshold> <output_r1> <output_r2>" << std::endl;
        std::cerr << "\n  barcode_file: Text file with one barcode per line" << std::endl;
        std::cerr << "  r1_fastq: Forward reads FASTQ file" << std::endl;
        std::cerr << "  r2_fastq: Reverse reads FASTQ file" << std::endl;
        std::cerr << "  barcode_start: Starting position of barcode in R1 (0-indexed)" << std::endl;
        std::cerr << "  barcode_length: Length of barcode sequence" << std::endl;
        std::cerr << "  strategy: Calling strategy (e.g., 4_7_mer_gpu_v1)" << std::endl;
        std::cerr << "  distance_measure: LEVENSHTEIN or SEQUENCE_LEVENSHTEIN" << std::endl;
        std::cerr << "  rejection_threshold: Maximum distance for assignment" << std::endl;
        std::cerr << "  output_r1: Output file for filtered R1 reads" << std::endl;
        std::cerr << "  output_r2: Output file for filtered R2 reads" << std::endl;
        return 1;
    }

    // Parse command line arguments
    std::string barcode_file = argv[1];
    std::string r1_fastq = argv[2];
    std::string r2_fastq = argv[3];
    int barcode_start = std::stoi(argv[4]);
    int barcode_length = std::stoi(argv[5]);
    std::string strategy = argv[6];
    std::string distance_measure_str = argv[7];
    int rejection_threshold = std::stoi(argv[8]);
    std::string output_r1 = argv[9];
    std::string output_r2 = argv[10];

    // Parse distance measure
    int distance_measure = -1;
    if (distance_measure_str == "LEVENSHTEIN") {
        distance_measure = LEVENSHTEIN_DISTANCE;
    } else if (distance_measure_str == "SEQUENCE_LEVENSHTEIN") {
        distance_measure = SEQUENCE_LEVENSHTEIN_DISTANCE;
    } else if (distance_measure_str == "PSEUDO_DISTANCE") {
        distance_measure = PSEUDO_DISTANCE;
    } else {
        std::cerr << "ERROR: Unsupported distance measure: " << distance_measure_str << std::endl;
        std::cerr << "Supported measures: LEVENSHTEIN, SEQUENCE_LEVENSHTEIN, PSEUDO_DISTANCE" << std::endl;
        return 1;
    }

    // Read barcode file
    std::cout << "Reading barcodes from: " << barcode_file << std::endl;
    
    // Read barcodes as strings first
    std::ifstream barcode_stream(barcode_file);
    if (!barcode_stream.is_open()) {
        std::cerr << "ERROR: Cannot open barcode file: " << barcode_file << std::endl;
        return 1;
    }
    
    std::vector<std::string> barcode_strings;  // Normalized for matching
    std::vector<std::string> original_barcode_strings;  // Original for output headers
    std::string line;
    while (std::getline(barcode_stream, line)) {
        if (!line.empty()) {
            // Store original barcode (trimmed of whitespace only)
            std::string original = line;
            original.erase(original.find_last_not_of(" \n\r\t")+1);
            original_barcode_strings.push_back(original);
            
            // Normalize to barcode_length for GPU matching
            // barcode_length should match SEQUENCE_LENGTH for proper GPU processing
            std::string normalized = normalize_sequence(line, barcode_length);
            barcode_strings.push_back(normalized);
        }
    }
    barcode_stream.close();
    
    // Convert to sequence objects
    std::vector<sequence> barcodes;
    for (const auto& bc_str : barcode_strings) {
        barcodes.emplace_back(bc_str);
    }
    
    std::cout << "Loaded " << barcodes.size() << " barcodes" << std::endl;

    // Read FASTQ files and extract barcodes
    std::cout << "Reading FASTQ files..." << std::endl;
    std::ifstream r1_file(r1_fastq);
    std::ifstream r2_file(r2_fastq);
    
    if (!r1_file.is_open()) {
        std::cerr << "ERROR: Cannot open R1 file: " << r1_fastq << std::endl;
        return 1;
    }
    if (!r2_file.is_open()) {
        std::cerr << "ERROR: Cannot open R2 file: " << r2_fastq << std::endl;
        return 1;
    }

    // Storage for all records
    std::vector<fastq_record> r1_records;
    std::vector<fastq_record> r2_records;
    std::vector<std::string> extracted_barcode_strings;
    std::vector<sequence> barcode_sequences;

    // Read all records
    fastq_record r1_rec, r2_rec;
    int record_count = 0;
    int short_read_count = 0;
    
    while (read_fastq_record(r1_file, r1_rec) && read_fastq_record(r2_file, r2_rec)) {
        r1_records.push_back(r1_rec);
        r2_records.push_back(r2_rec);
        
        // Extract barcode from R1 with padding if needed
        std::string barcode_str = extract_barcode(r1_rec.sequence, barcode_start, barcode_length);
        extracted_barcode_strings.push_back(barcode_str);
        
        // Track short reads
        if (r1_rec.sequence.length() < barcode_start + barcode_length) {
            short_read_count++;
        }
        
        // Normalize barcode to barcode_length for GPU processing
        // barcode_length should match SEQUENCE_LENGTH for proper sequence object creation
        std::string normalized = normalize_sequence(barcode_str, barcode_length);
        barcode_sequences.emplace_back(normalized);
        
        record_count++;
        if (record_count % 100000 == 0) {
            std::cout << "Processed " << record_count << " read pairs..." << std::endl;
        }
    }

    r1_file.close();
    r2_file.close();

    std::cout << "Total read pairs processed: " << record_count << std::endl;
    std::cout << "Short reads padded: " << short_read_count << std::endl;

    // Run barcode calling
    std::cout << "\n=== Running Barcode Calling ===" << std::endl;
    std::cout << "Strategy: " << strategy << std::endl;
    std::cout << "Distance measure: " << distance_measure_str << std::endl;
    std::cout << "Rejection threshold: " << rejection_threshold << std::endl;
    std::cout << "Barcode start position: " << barcode_start << std::endl;
    std::cout << "Barcode length: " << barcode_length << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();
    auto assignments = execute_strategy(strategy, barcodes, barcode_sequences, distance_measure, rejection_threshold);
    auto end_time = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double> duration = end_time - start_time;
    std::cout << "Barcode calling completed in " << duration.count() << " seconds" << std::endl;
    std::cout << "Time per read: " << (duration.count() / record_count) * 1000.0 << " ms" << std::endl;

    // Count assignments and rejections
    std::map<unsigned, int> barcode_counts;
    int assigned_count = 0;
    int rejected_count = 0;

    for (unsigned i = 0; i < assignments.size(); i++) {
        if (assignments[i] == UINT_MAX) {
            rejected_count++;
        } else {
            assigned_count++;
            barcode_counts[assignments[i]]++;
        }
    }

    // Write filtered output files
    std::cout << "\n=== Writing Output Files ===" << std::endl;
    std::ofstream out_r1(output_r1);
    std::ofstream out_r2(output_r2);

    if (!out_r1.is_open() || !out_r2.is_open()) {
        std::cerr << "ERROR: Cannot open output files" << std::endl;
        return 1;
    }

    int written_count = 0;
    for (unsigned i = 0; i < assignments.size(); i++) {
        if (assignments[i] != UINT_MAX) {
            // Write assigned reads with barcode appended to header
            // Use the ORIGINAL barcode (36bp) not the normalized one (34bp)
            std::string assigned_barcode = original_barcode_strings[assignments[i]];
            
            // Append barcode to the original header using underscore separator
            // Format: @READID_BARCODE remaining_header
            std::string r1_header = r1_records[i].header;
            std::string r2_header = r2_records[i].header;
            
            // Find first space in header to insert barcode before it
            size_t space_pos = r1_header.find(' ');
            if (space_pos != std::string::npos) {
                // Insert _BARCODE before the space
                r1_header.insert(space_pos, "_" + assigned_barcode);
                r2_header.insert(space_pos, "_" + assigned_barcode);
            } else {
                // No space found, append at end
                r1_header += "_" + assigned_barcode;
                r2_header += "_" + assigned_barcode;
            }
            
            // Write with modified header
            out_r1 << r1_header << "\n"
                   << r1_records[i].sequence << "\n"
                   << r1_records[i].plus << "\n"
                   << r1_records[i].quality << "\n";
            
            out_r2 << r2_header << "\n"
                   << r2_records[i].sequence << "\n"
                   << r2_records[i].plus << "\n"
                   << r2_records[i].quality << "\n";
            
            written_count++;
        }
    }

    out_r1.close();
    out_r2.close();

    std::cout << "Written " << written_count << " filtered read pairs" << std::endl;

    // Print statistics
    std::cout << "\n=== Final Statistics ===" << std::endl;
    std::cout << "Total reads: " << record_count << std::endl;
    std::cout << "Assigned reads: " << assigned_count << " (" 
              << (100.0 * assigned_count / record_count) << "%)" << std::endl;
    std::cout << "Rejected reads: " << rejected_count << " (" 
              << (100.0 * rejected_count / record_count) << "%)" << std::endl;
    std::cout << "Short reads (padded): " << short_read_count << " (" 
              << (100.0 * short_read_count / record_count) << "%)" << std::endl;
    std::cout << "\nBarcode distribution (top 10):" << std::endl;
    
    // Sort barcodes by count
    std::vector<std::pair<unsigned, int>> sorted_barcodes(barcode_counts.begin(), barcode_counts.end());
    std::sort(sorted_barcodes.begin(), sorted_barcodes.end(), 
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    int displayed = 0;
    for (const auto& pair : sorted_barcodes) {
        std::cout << "  Barcode " << pair.first << ": " << pair.second << " reads" << std::endl;
        if (++displayed >= 10) break;
    }

    std::cout << "\n=== SUCCESS ===" << std::endl;

    return 0;
}
