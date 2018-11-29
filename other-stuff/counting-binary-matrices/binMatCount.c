#include "stdio.h"
#include "stdlib.h"
#include "stdbool.h"
#include "gmp-6.1.0/gmp.h" // Because numbers get big


/*
 * A friend told me he wanted to count mxn binary matrices with support s (each column
 *  has exactly s non-zero entries), so I wrote this.
 */


// Fancy min macro that avoid double evaluation and checks if types agree (credits to StackOverflow person)
#define min(x, y) ({                \
    typeof(x) _min1 = (x);          \
    typeof(y) _min2 = (y);          \
    (void) (&_min1 == &_min2);      \
    _min1 < _min2 ? _min1 : _min2; })


/* Matrix type
 *
 * Store columns as unsigned integers
 * Always order columns from smallest to largest integer (to avoid permutations of columns,
 *  we will count the number of possible permutations after)
 * use bitwise xor to subtract columns from one another
 * use bitwise and with an integer having only one bit to check for 1s
 */
typedef struct Matrix {
    int height;
    int width;
    int supp;
    unsigned int *colArr;
} Matrix;


// Notes:
//   All functions assume the matrices are not null (unless otherwise specified)
//   No matrix or column vector may have more than sizeof(unsigned int) rows


// Returns a vector of length m with s ones in the "bottom" (least significant coordinates)
// (our implementation doesn't need m, but it might catch bugs and makes it easier to
//  change the implementation later)
unsigned int leastVector (int m, int s) {
    if (s > m) {
        fprintf(stderr, "leastColumn: Warning! [s too large]\n");
        fprintf(stderr, "  s (value: %d) is greater than m (value: %d). Setting s = m.\n", s, m);
        s = m;
    }

    unsigned int v = 0;
    for (int i = 0; i < s; i += 1) {
        v |= (1 << i);
    }
    return v;
}


// Creates a m-by-n binary Matrix where the bottom s entries of each column are 1s and the rest are 0s
//  (effectively, the "smallest" matrix with all columns having support s)
Matrix *MatrixCreate (int m, int n, int s) {
    if (s > m) {
        fprintf(stderr, "CreateMatrix: Warning! [s too large]\n");
        fprintf(stderr, "  s (value: %d) is greater than m (value: %d). Setting s = m.\n", s, m);
        s = m;
    }

    Matrix *Mat = (Matrix*)malloc(sizeof(Matrix));
    Mat->height = m;
    Mat->width = n;
    Mat->supp = s;

    unsigned int *M = (unsigned int*)malloc(sizeof(unsigned int) * n);
    unsigned int v = leastVector(m, s);

    for (int c = 0; c < n; c += 1) {
        M[c] = v;
    }
    Mat->colArr = M;
    return Mat;
}

// Free's memory allocated for matrix
void MatrixFree (Matrix *Mat) {
    free(Mat->colArr);
    free(Mat);
}


// Prints an m-by-n binary matrix
// Each line has a leading space because it was simpler (and kinda looks nice?).
void MatrixPrint (Matrix *Mat) {
    unsigned int *M = Mat->colArr;

    for (unsigned int row = 1 << (Mat->height-1); row > 0; row >>= 1) {
        for (int col = 0; col < Mat->width; col += 1) {
            if (M[col] & row) {
                printf(" 1");
            } else {
                printf(" 0");
            }
        }
        printf("\n");
    }
}


// Copies the entries of the binary matrix M into the binary matrix N
// Assumes both matrices are m-by-n
void MatrixCopy (Matrix *Mat1, Matrix *Mat2) {
    if (Mat1->height != Mat2->height || Mat1->width != Mat2->width) {
        fprintf(stderr, "MatrixCopy: ERROR! [Order of arguments differ]");
        fprintf(stderr, "Arguments are %dx%d and %dx%d respectively.", Mat1->height, Mat1->width, Mat2->height, Mat2->width);
        fprintf(stderr, "Aborting MatrixCopy");
        return;
    }
    unsigned int *M = Mat1->colArr, *N = Mat2->colArr;

    for (int c = 0; c < Mat1->width; c += 1) {
        N[c] = M[c];
    }
}


// "Increment" the given vector v (which has m entries)
bool VectorIncrement (unsigned int *v, int m) {
    // Look for 01 among bits
    int numBits = 0;
    unsigned int i = 1; //   01
    unsigned int ii = 2; //  10 = ~01
    unsigned int iii = 3; // 11
    while (ii <= (1 << (m - 1))) {
        if (((*v ^ ii) & iii) == iii) {
            // Bump up bit i, reset lower bits
            *v |= ii;
            do {
                *v &= ~i;
                i >>= 1;
            } while (i > 0);
            *v |= leastVector(m,numBits);
            return true;
        } else if (*v & i) {
            numBits += 1;
        }
        i <<= 1;
        ii <<= 1;
        iii <<= 1;
    }
    return false;
}


// "Increment" given matrix
bool MatrixIncrement (Matrix *Mat) {
    unsigned int *M = Mat->colArr;

    // Find a column to increment (default to first column)
    int c = 0;
    for (c = Mat->width - 1; c > 0; c -= 1) {
        if (M[c] < M[c-1]) {
            break;
        }
    }

    // Attempt to increment column
    bool success = VectorIncrement(&M[c], Mat->height);
    if (success) {
        unsigned int v = leastVector(Mat->height, Mat->supp);
        for (c += 1; c < Mat->width; c += 1) {
            M[c] = v;
        }
    }
    return success;
}


// Computes rank of m-by-n binary matrix M
// Modifies Mat
int MatrixRank (Matrix *Mat) {
    unsigned int *M = Mat->colArr;
    unsigned int row = 1 << (Mat->height-1);
    int col = 0;

    // Row-reduce matrix
    while (row > 0 && col < Mat->width) {
        if (M[col] == 0) {
            col += 1;
        }else if ((M[col] & row) > 0) {
            for (int i = 0; i < Mat->width; i += 1) {
                if ((M[i] & row) > 0 && i != col) {
                    M[i] ^= M[col];
                }
            }
            row >>= 1;
            col += 1;
        } else {
            row >>= 1;
        }
    }

    int rank = 0;
    for (int c = 0; c < Mat->width; c += 1) {
        if (M[c] > 0) {
            rank += 1;
        }
    }
    return rank;
}


// Compute the number of matrices which can be obtained by permuting the columns of Mat
void MatrixConjClassSize (Matrix *Mat, mpz_t *result) {
    // colUnused indicates unique columns
    // colClass indicates the number of duplicate columns for each unique column
    int colClass[Mat->width];
    bool colUnused[Mat->width];
    for (int i = 0; i < Mat->width; i+= 1) {
        colClass[i] = 0;
        colUnused[i] = true;
    }

    unsigned int *M = Mat->colArr;
    for (int c = 0; c < Mat->width; c += 1) {
        if (colUnused[c]) {
            colClass[c] = 1;
            for (int i = c+1; i < Mat->width; i += 1) {
                if (colUnused[i] && M[c] == M[i]) {
                    colUnused[i] = false;
                    colClass[c] += 1;
                }
            }
        }
    }

    mpz_t denom;
    mpz_init(denom);
    mpz_set_ui(denom,1);
    for (int i = 0; i < Mat->width; i += 1) {
        if (colClass[i] > 1) {
            mpz_fac_ui(*result, colClass[i]); // using result for intermediate value
            mpz_mul(denom, denom, *result);
        }
    }
    mpz_fac_ui(*result, Mat->width);
    mpz_tdiv_q(*result, *result, denom);
    mpz_clear(denom);
}



int main (int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s\n", argv[0]);
        fprintf(stderr, "Count m-by-n binary matrices whose columns have support s\n");
        fprintf(stderr, "Expects three positive integers m, n, and s (no bigger than 30 just to be safe)\n.");
        return 1;
    }
    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    int s = atoi(argv[3]);
    printf("m = %d, n = %d, s = %d\n", m, n, s);

    if (s > m) {
        fprintf(stderr, "bmcount: Warning! [s too large]\n");
        fprintf(stderr, "  s (value: %d) is greater than m (value: %d). Setting s = m.\n", s, m);
        s = m;
    }

    int rs = min(m,n) + 1;
    mpz_t z; // Only used for intermediate computations
    mpz_t *ranks = (mpz_t*)malloc(sizeof(mpz_t) * rs);

    // Initialize
    mpz_init(z);
    for (int r = 0; r < rs; r += 1) {
        mpz_init(ranks[r]);
    }

    Matrix *M = MatrixCreate(m,n,s);
    Matrix *N = MatrixCreate(m,n,s);

    do {
        MatrixCopy(M, N);
        int r = MatrixRank(N);
        MatrixConjClassSize(M, &z);
        mpz_add(ranks[r], ranks[r], z);
    } while (MatrixIncrement(M));

    MatrixFree(M);
    MatrixFree(N);

    for (int r = 0; r < rs; r += 1) {
        gmp_printf("%d: %Zd\n", r, ranks[r]);
    }
}