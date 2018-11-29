#include <stdio.h>
#include <stdlib.h>

/*
 * A simple utility which receives graphs (in g6 formatm one per line) and outputs the
 *  graphs (in g6 format) which are connected. (The graphs may have different numbers
 *  of vertices, but no more than MAX_VERTICES.)
 * See the documentation for g6 format (from B. McKay and A. Piperino's nauty & traces)
 */

// Uncomment to print the total number of graphs read and the number of connected graphs
//  to stderr once EOF is reached
// #define PRINT_STATS

#define LINE_MAX_LEN 256
// I don't know how big a graph has to be before its g6 string exceeds 256 characters
//  but a quick estimate says that 20 will definitely be safe
#define MAX_VERTICES 20

// g6 format magic number (see documentation)
#define G6_START_CHAR 63

struct _Graph {
    int nVerts;
    char *adjMat;
};
typedef struct _Graph Graph;


// Returns 1 if the character is a null terminator or a newline character, 0 otherwise
int isLineEnd(char c) {
    return (c == 0) || (c == '\n');
}


// Print contents of g (for debugging)
void printGraph (Graph *g) {
    printf("Graph on %d vertices:\n", g->nVerts);

    for (int i = 0; i < g->nVerts-1; i += 1) {
        for (int j = 1; j < g->nVerts; j += 1) {
            if (j <= i) {
                printf("  ");
            } else {
                printf("%d ", g->adjMat[(j * (j - 1)) / 2 + i]);
            }
        }
        printf("\n");
    }
}


// Returns 1 if the graph is connected, 0 otherwise
int isConnected (Graph *G) {
    int n = G->nVerts;
    int reachable[n];
    reachable[0] = 1;
    int pos = 0;
    for (int i = 1; i < n; i += 1) {
        reachable[i] = G->adjMat[pos];
        pos += i;
    }

    // Run Prim's Algorithm
    int reached_new_vertex;
    do {
        reached_new_vertex = 0;
        for (int v = 1; v < n; v += 1) {
            if (reachable[v] == 0) {
                int pos = (v * (v - 1)) / 2 + 1;
                int w = 1;
                // First we look down the v-th column until we hit the diagonal
                while(w != v) {
                    if (reachable[w] != 0 && G->adjMat[pos] != 0) {
                        reachable[v] = 1;
                        reached_new_vertex = 1;
                        break;
                    }
                    pos += 1;
                    w += 1;
                }
                if (reached_new_vertex != 0) {
                    break;
                }
                pos += v;
                v += 1;
                // Now we look along the w-th row
                //  (which is the rest of the column we were looking at before)
                while(v < n) {
                    if (reachable[w] != 0 && G->adjMat[pos] != 0) {
                        reachable[v] = 1;
                        reached_new_vertex = 1;
                        break;
                    }
                    pos += v;
                    v += 1;
                }
                if (reached_new_vertex != 0) {
                    break;
                }
            }
        }
    } while (reached_new_vertex != 0);


    // Are all vertices reachable
    for (int i = 0; i < n; i += 1) {
        if (reachable[i] == 0)
            return 0;
    }
    return 1;
}


int main (void) {
    char line[LINE_MAX_LEN+1];
    line[0] = line[LINE_MAX_LEN] = 0; // Make sure the string is null terminated
    Graph G;
    G.adjMat = NULL;
    int lastN = 0;
    #ifdef PRINT_STATS
        int totalGraphs = 0;
        int totalConnected = 0;
    #endif

    // Read the line and make a graph (carefully)
    while(fgets(line, LINE_MAX_LEN, stdin) != NULL) {
        if (isLineEnd(line[0])) {
            continue;
        }

        int n = line[0] - G6_START_CHAR;

        if (n > MAX_VERTICES) {
            fprintf(stderr, "ERROR: (g6connected) Found graph with %d vertices (greater than max of %d\n", n , MAX_VERTICES);
            continue;
        } else if (n < 0) {
            fprintf(stderr, "ERROR: (g6connected) Found graph with %d vertices\n", n);
            continue;
        }

        G.nVerts = n;
        const int adjMatSize = ((n * (n - 1)) / 2);

        if (n != lastN && G.adjMat != NULL) {
            free(G.adjMat);
            G.adjMat = NULL;
        }
        if (G.adjMat == NULL) {
            G.adjMat = (char*)malloc(sizeof(char) * adjMatSize);
            lastN = n;
        }

        int matPos = 0;
        int linePos = 1;
        while (!isLineEnd(line[linePos])) {
            if (matPos >= adjMatSize) {
                break;
            }
            char c = line[linePos++] - G6_START_CHAR;
            for (int i = 5; i >= 0 && matPos < adjMatSize; i -= 1) {
                G.adjMat[matPos++] = (c >> i) & 1;
            }
        }
        if (matPos >= adjMatSize) {
            if (!isLineEnd(line[linePos])) {
                fprintf(stderr, "ERROR: (g6connected) g6string \"%s\"is too long", line);
                free(G.adjMat);
                G.adjMat = NULL;
                continue;
            }
        } else if (matPos <= adjMatSize) {
            fprintf(stderr, "ERROR: (g6connected) g6string \"%s\" ended prematurely", line);
            free(G.adjMat);
            G.adjMat = NULL;
            continue;
        }


        // Now that we have a graph, let's check if it is connected.
        if (isConnected(&G) == 1) {
            fputs(line, stdout);
            #ifdef PRINT_STATS
                totalConnected += 1;
            #endif
        }
        #ifdef PRINT_STATS
            totalGraphs += 1;
        #endif
    }

    if (G.adjMat != NULL) {
        free(G.adjMat);
    }

    #ifdef PRINT_STATS
        fprintf(stderr, ">Found %d connected graphs out of %d.\n", totalConnected, totalGraphs);
    #endif
}
