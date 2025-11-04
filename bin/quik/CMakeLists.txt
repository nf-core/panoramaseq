cmake_minimum_required(VERSION 3.22)
project(barcode_calling LANGUAGES CXX CUDA)

set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CUDA_ARCHITECTURES 60;80)

# Allow parameters to be passed from command line or use defaults
set(SEQUENCE_LENGTH "36" CACHE STRING "Length of barcode sequences for GPU processing")
set(REJECTION_THRESHOLD "8" CACHE STRING "Default rejection threshold for barcode assignment")

# Pass these as compile-time definitions to all source files
add_compile_definitions(
    SEQUENCE_LENGTH=${SEQUENCE_LENGTH}
    REJECTION_THRESHOLD=${REJECTION_THRESHOLD}
)

message(STATUS "Compiling with SEQUENCE_LENGTH=${SEQUENCE_LENGTH}")
message(STATUS "Compiling with REJECTION_THRESHOLD=${REJECTION_THRESHOLD}")

find_package(OpenMP REQUIRED)

add_executable(simulate_errors
        src/simulate_errors.cu
)

add_executable(benchmark_barcode_calling
        src/benchmark_barcode_calling.cu
        src/two_step_k_mer_filtered_calling_gpu_v1.cuh
)

add_executable(benchmark_sl_distance
        src/benchmark_sl_distance.cu
)

add_executable(single_strategy_benchmark_fastq_paired
        src/single_strategy_benchmark_fastq_paired.cu
        src/two_step_k_mer_filtered_calling_gpu_v1.cuh
)

# Aktivieren von OpenMP für den Host-Compiler
if(OpenMP_CXX_FOUND)
    target_link_libraries(benchmark_barcode_calling PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(benchmark_sl_distance PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(single_strategy_benchmark_fastq_paired PUBLIC OpenMP::OpenMP_CXX)
    target_compile_options(benchmark_barcode_calling PUBLIC ${OpenMP_CXX_FLAGS})
    target_compile_options(benchmark_sl_distance PUBLIC ${OpenMP_CXX_FLAGS})
    target_compile_options(single_strategy_benchmark_fastq_paired PUBLIC ${OpenMP_CXX_FLAGS})
endif()

# Aktivieren von OpenMP für den CUDA-Compiler
set_target_properties(benchmark_barcode_calling PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
)

# Aktivieren von OpenMP für den CUDA-Compiler
set_target_properties(benchmark_sl_distance PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
)

# Aktivieren von OpenMP für den CUDA-Compiler
set_target_properties(single_strategy_benchmark_fastq_paired PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
)

target_compile_options(benchmark_barcode_calling
        PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${OpenMP_CXX_FLAGS}>
)

target_compile_options(benchmark_sl_distance
        PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${OpenMP_CXX_FLAGS}>
)

target_compile_options(single_strategy_benchmark_fastq_paired
        PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${OpenMP_CXX_FLAGS}>
)
