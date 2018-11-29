#include <stdio.h>
#include <stdlib.h>

/*
 * A simple utility which receives a graph (in g6 format) and outputs
 *  the g6 string of the complement graph
 * See the documentation for g6 format (from B. McKay and A. Piperino's nauty & traces)
 */

#define LINE_MAX_LEN 256
// I don't know how big a graph has to be before its g6 string exceeds 256 characters
//  but a quick estimate says that 20 will definitely be safe
#define MAX_VERTICES 20

// g6 format magic number (see documentation)
#define G6_START_CHAR 63

const char ONES[7] = {0,              // 000000
                      32,             // 100000
                      32+16,          // 110000
                      32+16+8,        // 111000
                      32+16+8+4,      // 111100
                      32+16+8+4+2,    // 111110
                      32+16+8+4+2+1}; // 111111


// Returns 1 if the character is a null terminator or a newline character, 0 otherwise
int isLineEnd(char c) {
    return (c == 0) || (c == '\n');
}

int main (void) {

    char line[LINE_MAX_LEN+1];
     // Make sure the string are null terminated
    line[0] = line[LINE_MAX_LEN] = 0;

    int numGraphs = 0;

    // Read the line and xor the bit appropriately
    while(fgets(line, LINE_MAX_LEN, stdin) != NULL) {
        if (isLineEnd(line[0])) {
            continue;
        }

        int n = line[0] - G6_START_CHAR;

        if (n > MAX_VERTICES) {
            fprintf(stderr, "ERROR: (g6compl) Found graph with %d vertices (greater than max of %d\n", n , MAX_VERTICES);
            continue;
        } else if (n < 0) {
            fprintf(stderr, "ERROR: (g6compl) Found graph with %d vertices\n", n);
            continue;
        }

        int bitsLeft = n * (n - 1) / 2;

        int linePos = 1;
        while (!isLineEnd(line[linePos])) {
            char c = line[linePos] - G6_START_CHAR;

            if (bitsLeft >= 6) {
                line[linePos] = (c ^ ONES[6]) + G6_START_CHAR;
                bitsLeft -= 6;
            } else if (bitsLeft > 0) {
                line[linePos] = (c ^ ONES[bitsLeft]) + G6_START_CHAR;
                bitsLeft = 0;
            } else {
                fprintf(stderr, "ERROR: (g6compl) Expected graph6 string to end, found %c\n", c);
            }

            linePos += 1;
        }

        if (bitsLeft > 0) {
            fprintf(stderr, "ERROR: (g6compl) reached EOF prematurely, expected %d more bits\n", bitsLeft);
        } else {
            fputs(line, stdout);
            numGraphs += 1;
        }
    }

    fprintf(stderr, ">%d graph complements generated\n", numGraphs);
}
