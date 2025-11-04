//
// Created by steffen on 17.07.24.
//

#ifndef INC_2OPT_CONSTANTS_H
#define INC_2OPT_CONSTANTS_H

#define MAX_INDEX_SIZE 100

// SEQUENCE_LENGTH and REJECTION_THRESHOLD can be set at compile time via CMake
// If not provided by CMake, use these defaults
#ifndef SEQUENCE_LENGTH
#define SEQUENCE_LENGTH 36
#endif

#ifndef REJECTION_THRESHOLD
#define REJECTION_THRESHOLD 8
#endif

#define PSEUDO_DISTANCE_WINDOW_SIZE 5

// options for the distance measure
#define LEVENSHTEIN_DISTANCE 0
#define SEQUENCE_LEVENSHTEIN_DISTANCE 1

#endif //INC_2OPT_CONSTANTS_H
