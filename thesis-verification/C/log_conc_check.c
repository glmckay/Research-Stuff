#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include "gmp-6.1.0/gmp.h" // Change accordingly

/*
 * Verification used for section 6.1 (titled "Log-Concavity of $(P(G;k,\ell))_{\ell=0}^k$") of my master's thesis
 *
 * In the language of my master's thesis, the sequence of b's for the graph G and integer k is the sequence
 * (P(G;k,i))_{i=0}^k where P(G;x,y) is the two-variable chromatic polynomial of Dohmen, Poenitz, and Tittmann.
 * To compute these values, we run through all subsets of V(G) and record, for each i=0,...,k, the number
 * of stable sets of G cardinality i. We store these in the array 'results' and then compute the "b's" for
 * each k from 1 to MAX_COLOURS
 *
 * Requires the GNU Multiple Precision Arithmetic Library (GMP)
*/

#define N 9 // number of vertices
#define ADJ_MAT_SIZE ((N * (N - 1)) / 2)
#define RESULTS_ROWS ((N / 2) + 1)
#define RESULTS_COLS (N + 1)
#define RESULTS_SIZE (RESULTS_ROWS * RESULTS_COLS)
#define MAX_COLOURS (N + 2)

// Uncomment to write the number of stable sets of each size for each graph
// #define WRITE_RESULTS_TO_FILE

// g6 format magic number (see documentation)
#define G6_START_CHAR 63

// File permissions when we create a file (all can read, user can write)
#define FILE_PERMISSIONS S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH

// Min macro that avoids double evaluation and checks if types agree (thanks stackoverflow)
#define min(x, y) ({                \
    typeof(x) _min1 = (x);          \
    typeof(y) _min2 = (y);          \
    (void) (&_min1 == &_min2);      \
    _min1 < _min2 ? _min1 : _min2; })


// Our lists will be immutable and reference counted
// We represent an (ordered) set of integers by a list and an (ordered) partition by a list of lists

// List of integers
typedef struct ILst {
    int refcount;
    int first;
    struct ILst *rest;
} ILst;

// List of lists
typedef struct LLst {
    int refcount;
    ILst *first;
    struct LLst *rest;
} LLst;

// List of lists of lists
typedef struct LLLst {
    int refcount;
    LLst *first;
    struct LLLst *rest;
} LLLst;

/* Graph structure:
 * We store a graph as number of vertices and an adjacency matrix
 *  (for the purposes of this program this should always be N)
 * The adjacency matrix is stored as an array obtained by concatenating the columns of the upper triangular
 *  part of the adjacency matrix
 *  (i.e. for vertices v,w, adjMat[(u * (u - 1)) / 2 + v] is 1 if vw is an edge and 0 otherwise)
 */
typedef struct Graph {
    int nVerts;
    char *adjMat;
} Graph;



/*  (Muliple functions)
    "Constructors" for each lists
    Appends elem to the list
*/
ILst *consILst (int elem, ILst *lst) {
    ILst *newCell = (ILst*)malloc(sizeof(ILst));
    newCell->refcount = 0;
    newCell->first = elem;
    newCell->rest = lst;
    if (lst != NULL) {
        newCell->rest->refcount += 1;
    }
    return newCell;
}

LLst *consLLst (ILst *elem, LLst *lst) {
    LLst *newCell = (LLst*)malloc(sizeof(LLst));
    newCell->refcount = 0;
    newCell->first = elem;
    newCell->rest = lst;
    if (lst != NULL) {
        newCell->rest->refcount += 1;
    }
    return newCell;
}

LLLst *consLLLst (LLst *elem, LLLst *lst) {
    LLLst *newCell = (LLLst*)malloc(sizeof(LLLst));
    newCell->refcount = 0;
    newCell->first = elem;
    newCell->rest = lst;
    if (lst != NULL) {
        newCell->rest->refcount += 1;
    }
    return newCell;
}


/*  (Muliple functions)
    Print lists (for debugging)
*/
void printILst (ILst *lst) {
    printf("(");
    if (lst != NULL) {
        printf("%d", lst->first);
        lst = lst->rest;
        while (lst != NULL) {
            printf(" %d", lst->first);
            lst = lst->rest;
        }
    }
    printf(")");
}

void printLLst (LLst *lst) {
    printf("(");
    if (lst != NULL) {
        printILst(lst->first);
        lst = lst->rest;
        while (lst != NULL) {
            printf(" ");
            printILst(lst->first);
            lst = lst->rest;
        }
    }
    printf(")");
}

void printLLLst (LLLst *lst) {
    printf("(");
    if (lst != NULL) {
        printLLst(lst->first);
        lst = lst->rest;
        while (lst != NULL) {
            printf(" ");
            printLLst(lst->first);
            lst = lst->rest;
        }
    }
    printf(")");
}


/*  (Muliple functions)
    Appends elem to each element of the list
*/
LLst *consAllILst (int elem, LLst *lst) {
    LLst *pos = lst;
    while (pos != NULL) {
        pos->first = consILst(elem, pos->first);
        pos = pos->rest;
    }
    return lst;
}

LLLst *consAllLLst (ILst *elem, LLLst *lst) {
    LLLst *pos = lst;
    while (pos != NULL) {
        pos->first = consLLst(elem, pos->first);
        elem->refcount += 1;
        pos = pos->rest;
    }
    return lst;
}


/*  (Muliple functions)
    Functions to free lists
*/
void freeILst (ILst *lst) {
    while (lst != NULL) {
        lst->refcount -= 1;
        if (lst->refcount > 0) {
            break;
        }
        ILst *rest = lst->rest;
        free(lst);
        lst = rest;
    }
}

void freeLLst (LLst *lst) {
    while (lst != NULL) {
        lst->refcount -= 1;
        if (lst->refcount > 0) {
            break;
        }
        LLst *rest = lst->rest;
        freeILst(lst->first);
        free(lst);
        lst = rest;
    }
}

void freeLLLst (LLLst *lst) {
    while (lst != NULL) {
        lst->refcount -= 1;
        if (lst->refcount > 0) {
            break;
        }
        LLLst *rest = lst->rest;
        freeLLst(lst->first);
        free(lst);
        lst = rest;
    }
}


/*  (Muliple functions)
    Returns the length of the list
*/
int lengthLLst (LLst *lst) {
    int len = 0;
    while (lst != NULL) {
        len += 1;
        lst = lst->rest;
    }
    return len;
}


int lengthLLLst (LLLst *lst) {
    int len = 0;
    while (lst != NULL) {
        len += 1;
        lst = lst->rest;
    }
    return len;
}


// Returns the result of appending lst2 to the end of lst1 (both LLLsts)
LLLst *appendLLLst (LLLst *lst1, LLLst *lst2) {
    if (lst1 == NULL) {
        return lst2;
    } else if (lst2 == NULL) {
        return lst1;
    }

    LLLst *newLst = consLLLst(lst1->first, NULL);
    LLLst *pos = newLst;
    lst1 = lst1->rest;
    while (lst1 != NULL) {
        pos->rest = consLLLst(lst1->first, NULL);
        pos = pos->rest;
        lst1 = lst1->rest;
    }
    pos->rest = lst2;
    return newLst;
}


// Add elem tor prtn in all possible ways (i.e. add it to each part and as its own part)
LLLst *addToPrtn (int elem, LLst *prtn) {
    if (prtn == NULL) {
        return consLLLst(consLLst(consILst(elem, NULL), NULL), NULL);
    } else {
        return consLLLst(consLLst(consILst(elem, prtn->first), prtn->rest), consAllLLst(prtn->first, addToPrtn(elem, prtn->rest)));
    }
}


// Add elem to each partition in prtns (a list of partitions) in all possible ways
LLLst *addToAllPrtns (int elem, LLLst *prtns) {
    LLLst *result = NULL;
    while (prtns != NULL) {
        result = appendLLLst(addToPrtn(elem, prtns->first), result);
        prtns = prtns->rest;
    }
    return result;
}


// Produces a list of all set partitions of the set {0,...,n-1}
LLLst *partitionsList (int n) {
    if (n == 0) {
        return consLLLst(NULL, NULL);
    } else {
        return addToAllPrtns(n-1, partitionsList(n-1));
    }
}


// Returns 1 if the set is stable in the graph g, 0 otherwise
int isStable(ILst *set, Graph *g) {
    // For each vertex u in the set, check if u is adjacent to any remaining vertex in the set
    while(set != NULL) {
        int u = set->first;
        set = set->rest;
        int colStartIndx = (u * (u - 1)) / 2;
        for (ILst *pos = set; pos != NULL; pos = pos->rest) {
            if (g->adjMat[colStartIndx + pos->first]) {
                return 0;
            }
        }
    }
    return 1;
}


// Returns the number of sets in the given list which are stable in g
int numStableSets (LLst *sets, Graph *g) {
    int n = 0;
    while (sets != NULL) {
        n += isStable(sets->first, g);
        sets = sets->rest;
    }
    return n;
}


// Sets the i-th entry of 'results' to be the number of partitions (provided as prtns) of V(g)
// with exactly i stable parts
// Assumes 'results' is initialised to zero and has RESULTS_SIZE entries
void countPartitions (LLLst *prtns, Graph *g, int *results) {
    while (prtns != NULL) {
        int stbl = numStableSets(prtns->first, g);
        results[((N + 1) * (lengthLLst(prtns->first) - stbl)) + stbl] += 1;
        prtns = prtns->rest;
    }
}

// Read a graph (in g6 format) from the file descriptor 'in'
// Assumes the graph being read has no more than 62 vertices and that the adjacency matrix has ADJ_MAT_SIZE entries
// For more details see the documentation for g6 format (from B. McKay and A. Piperino's nauty & traces)
int readGraph (int in, Graph *g) {
    int c = 0;
    if (read(in, &c, sizeof(char)) == 0) {
        return 0; // Reached EOF
    }
    g->nVerts = c - G6_START_CHAR;

    int pos = 0;
    while (read(in, &c, sizeof(char)) != 0) {
        if (c == '\n') {
            break;
        }

        c = c - G6_START_CHAR;
        for (int i = 5; i >= 0 && pos < ADJ_MAT_SIZE; i -= 1) {
            g->adjMat[pos++] = (c >> i) & 1;
        }
    }
    return 1;
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


// Write the results to 'out'
// Assumes results has RESULTS_SIZE entries
void writeGraphResults (int out, int *results) {
    write(out, results, sizeof(int) * RESULTS_SIZE);
}


// Compute n! and return as z
void fact (mpz_t z, unsigned int n) {
    mpz_set_ui(z,1);
    while (n > 0) {
        mpz_mul_ui(z, z, n);
        n -= 1;
    }
}


// Compute the falling factorial (n)_k and return as z
void fallingFact (mpz_t z, unsigned int n, unsigned int k) {
    if (k > n) {
        mpz_set_ui(z, 0);
    } else {
        mpz_set_ui(z,1);
        for (unsigned int i = n; i > n - k; i -= 1) {
            mpz_mul_ui (z, z, i);
        }
    }
}


// For each i, returns the number of colour assignments with the first i colours stable
// as the i-th entry of b
// Assumes b has MAX_COLOURS entries
void computeBs (int k, int *results, mpz_t *b) {
    mpz_t z1, z2;
    mpz_init(z1);
    mpz_init(z2);
    for (unsigned int i = 0; i < MAX_COLOURS; i += 1) {
        mpz_set_ui(b[i], 0);
        if (i > k) {
            continue;
        }

        int resPos = 0;
        // Note: If the order of the loops is changed for some reason, don't keep the min
        //  in the exit condition for the ns for loop
        for (unsigned int ns = 0; ns < min((unsigned int) RESULTS_ROWS, k - i + 1); ns += 1) {
            for (unsigned int s = 0; s < RESULTS_COLS; s += 1) {
                if (ns + s <= k) {
                    fallingFact(z1, k, ns + s);
                    mpz_addmul_ui(b[i], z1, (unsigned int) results[resPos]);
                }
                resPos += 1;
            }
        }
    }
}



// Check if b's are log-concave
// In case they are not, print the line (g6 string) as well as the particular values of k and i
// Assumes b has MAX_COLOURS entries
void checkLogConc(int line, int k, mpz_t *b) {
    mpz_t bb, ac;
    mpz_init(bb);
    mpz_init(ac);
    for (int x = 0, y = 1, z = 2; z < MAX_COLOURS; x += 1, y += 1, z += 1) {
        mpz_mul(ac, b[x], b[z]);
        mpz_mul(bb, b[y], b[y]);
        if (mpz_cmp(bb,ac) < 0) {
            printf("line = %d -- k = %d -- i = %d\n", line, k, y);
        }
    }
}


// Print the b's (for debugging)
// Assumes b has MAX_COLOURS entries
void printBs(mpz_t *b) {
    printf("[");
    for (int i = 0; i < MAX_COLOURS; i += 1) {
        gmp_printf("%Zd", b[i]);
        if (i < MAX_COLOURS - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}


// Open the file of g6 strings
// Assumes graphs with n vertices are in the file 'graphs_data/connected/graphs_n.g6'
int openGraphDataFile(int n) {
    // length("graph_data/graphs_") = 18
    // length(".g6") = 3
    // We won't be worried about numbers of vertices larger than 99 so our filepath
    //  should be a string of length at most 24 = 18 + 2 + 3 + 1 (null terminator)
    char *filepath = (char*)malloc(sizeof(char) * 24);
    filepath[24-1] = filepath[24] = '\0';

    strcpy(filepath, "graph_data/connected/graphs_");
    if (0 <= n && n <= 99) {
        int pos = 18;
        if (n > 9) {
            filepath[pos++] = (n / 10) + '0';
        }
        filepath[pos++] = (n % 10) + '0';
        strcpy(&filepath[pos], ".g6");
    } else {
        printf("Error in openGraphDataFile: Value of 'n', %d, not in range 0..99.", n);
        free(filepath);
        return -1;
    }

    int fileDesc = open(filepath, O_RDONLY);
    if (fileDesc < 0) {
        printf("Failed to open intput file.\n");
        printf("%s\n", strerror(errno));
        printf("Filepath: %s", filepath);
    }
    free(filepath);
    return fileDesc;
}


int main () {
    // Initialise everything
    Graph g;
    g.adjMat = (char*)malloc(sizeof(char) * ADJ_MAT_SIZE);
    mpz_t *b = (mpz_t*)malloc(sizeof(mpz_t) * MAX_COLOURS);
    for (int i = 0; i < MAX_COLOURS; i += 1) {
        mpz_init(b[i]);
    }

    int in = openGraphDataFile(N);
    #ifdef WRITE_RESULTS_TO_FILE
        int out = open("graph_data/graphs_10_data_c.txt", O_WRONLY | O_CREAT, FILE_PERMISSIONS);
    #endif

    int *results = (int*)malloc(sizeof(int) * RESULTS_SIZE);

    if (in >= 0) {
        LLLst *prtns = partitionsList(N);
        int count = 0;
        int printCount = 0;
        while (readGraph(in, &g) != 0) {
            count += 1;

            for (int i = 0; i < RESULTS_SIZE; i += 1) {
                results[i] = 0;
            }

            countPartitions(prtns, &g, results);
            #ifdef WRITE_RESULTS_TO_FILE
                writeGraphResults(out, results)
            #endif
            for (int k = 0; k < MAX_COLOURS; k += 1) {
                computeBs(k, results, b);
                checkLogConc(count, k, b);
            }

            // Show progress in stdout
            printCount += 1;
            if (printCount >= 1000) {
                printf("At line: %d\n", count);
                fflush(stdout);
                printCount = 0;
            }
        }
        freeLLLst(prtns);
    }

    close(in);
    #ifdef WRITE_RESULTS_TO_FILE
        close(out);
    #endif
    free(b);
    free(results);
    free(g.adjMat);
}