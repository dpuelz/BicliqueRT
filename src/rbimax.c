/*
 * Bimax algorithm for finding binary inclusion-maximal biclusters
 *
 * Based on the algorithm by Prelic, Bleuler, Zimmermann, Wille, Buhlmann,
 * Gruissem, Hennig, Thiele, and Zitzler (2006).
 * "A systematic comparison and evaluation of biclustering methods for
 * gene expression data", Bioinformatics 22(9):1122-1129.
 *
 * This implementation is derived from the biclust R package (GPL-2 | GPL-3)
 * by Sebastian Kaiser, Rodrigo Santamaria, Roberto Theron, Ludwig Fabricius,
 * and Friedrich Leisch.
 *
 * License: GPL (>= 2)
 */

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <string.h>

/* --- Bimax algorithm core --- */

/*
 * The Bimax algorithm works on a binary matrix and finds inclusion-maximal
 * biclusters. The approach is a divide-and-conquer strategy:
 *
 * 1. Pick a column and divide rows into those with 1 in that column and those without
 * 2. Recursively find biclusters in the subproblems
 * 3. Merge overlapping results
 */

static void bimax_recurse(
    int *data,       /* binary matrix (nrow x ncol), stored column-major */
    int nrow,
    int ncol,
    int *row_set,    /* current row indices */
    int num_rows,    /* number of rows in row_set */
    int *col_set,    /* current column indices */
    int num_cols,    /* number of columns in col_set */
    int minr,        /* minimum rows in bicluster */
    int minc,        /* minimum cols in bicluster */
    int max_biclusters, /* max number of biclusters to find */
    int *result_rows,   /* output: rows per bicluster (nrow x max_biclusters) */
    int *result_cols,   /* output: cols per bicluster (ncol x max_biclusters) */
    int *num_found,     /* number of biclusters found so far */
    int total_nrow,     /* original nrow for matrix indexing */
    int total_ncol      /* original ncol for matrix indexing */
)
{
    if (*num_found >= max_biclusters) return;
    if (num_rows < minr || num_cols < minc) return;

    /* Check if current submatrix is all-ones */
    int all_ones = 1;
    for (int j = 0; j < num_cols && all_ones; j++) {
        for (int i = 0; i < num_rows && all_ones; i++) {
            if (data[row_set[i] + col_set[j] * total_nrow] == 0) {
                all_ones = 0;
            }
        }
    }

    if (all_ones) {
        /* Found a bicluster */
        if (num_rows >= minr && num_cols >= minc) {
            int idx = *num_found;
            for (int i = 0; i < num_rows; i++) {
                result_rows[row_set[i] + idx * total_nrow] = 1;
            }
            for (int j = 0; j < num_cols; j++) {
                result_cols[col_set[j] + idx * total_ncol] = 1;
            }
            (*num_found)++;
        }
        return;
    }

    /* Find the column with the most zeros (best split) */
    int best_col = -1;
    int best_zeros = -1;
    for (int j = 0; j < num_cols; j++) {
        int zeros = 0;
        for (int i = 0; i < num_rows; i++) {
            if (data[row_set[i] + col_set[j] * total_nrow] == 0) {
                zeros++;
            }
        }
        if (zeros > 0 && (best_col == -1 || zeros < best_zeros)) {
            best_col = j;
            best_zeros = zeros;
        }
    }

    if (best_col == -1) return; /* shouldn't happen since not all_ones */

    int split_col = col_set[best_col];

    /* Divide rows: those with 1 in split_col, and those with 0 */
    int *rows_one = (int *)R_alloc(num_rows, sizeof(int));
    int *rows_zero = (int *)R_alloc(num_rows, sizeof(int));
    int n_one = 0, n_zero = 0;

    for (int i = 0; i < num_rows; i++) {
        if (data[row_set[i] + split_col * total_nrow] == 1) {
            rows_one[n_one++] = row_set[i];
        } else {
            rows_zero[n_zero++] = row_set[i];
        }
    }

    /* Divide columns: those where all rows_one have 1, and the rest */
    int *cols_all = (int *)R_alloc(num_cols, sizeof(int));
    int *cols_some = (int *)R_alloc(num_cols, sizeof(int));
    int n_all = 0, n_some = 0;

    for (int j = 0; j < num_cols; j++) {
        int has_zero = 0;
        for (int i = 0; i < n_one; i++) {
            if (data[rows_one[i] + col_set[j] * total_nrow] == 0) {
                has_zero = 1;
                break;
            }
        }
        if (!has_zero) {
            cols_all[n_all++] = col_set[j];
        } else {
            cols_some[n_some++] = col_set[j];
        }
    }

    /* Recurse: rows_one with cols_all (the "all-ones" subproblem) */
    if (n_one >= minr && n_all >= minc) {
        bimax_recurse(data, nrow, ncol,
                      rows_one, n_one, cols_all, n_all,
                      minr, minc, max_biclusters,
                      result_rows, result_cols, num_found,
                      total_nrow, total_ncol);
    }

    /* Recurse: all rows with cols_some (columns not fully covered) */
    if (num_rows >= minr && n_some >= minc) {
        bimax_recurse(data, nrow, ncol,
                      row_set, num_rows, cols_some, n_some,
                      minr, minc, max_biclusters,
                      result_rows, result_cols, num_found,
                      total_nrow, total_ncol);
    }
}


/*
 * R-callable Bimax function
 *
 * Arguments:
 *   data_vec   - integer vector (binary matrix in column-major order)
 *   nrow_r     - number of rows
 *   ncol_r     - number of columns
 *   minr_r     - minimum number of rows in bicluster
 *   minc_r     - minimum number of columns in bicluster
 *   number_r   - maximum number of biclusters to find
 *   result_rows - output integer vector (nrow x number), 0/1 membership
 *   result_cols - output integer vector (ncol x number), 0/1 membership
 *   num_found  - output integer, number of biclusters found
 */
void bimax_c(int *data_vec, int *nrow_r, int *ncol_r,
             int *minr_r, int *minc_r, int *number_r,
             int *result_rows, int *result_cols, int *num_found)
{
    int nrow = *nrow_r;
    int ncol = *ncol_r;
    int minr = *minr_r;
    int minc = *minc_r;
    int max_biclusters = *number_r;

    *num_found = 0;

    /* Initialize result matrices to zero */
    memset(result_rows, 0, nrow * max_biclusters * sizeof(int));
    memset(result_cols, 0, ncol * max_biclusters * sizeof(int));

    /* Create initial row and column sets */
    int *row_set = (int *)R_alloc(nrow, sizeof(int));
    int *col_set = (int *)R_alloc(ncol, sizeof(int));

    for (int i = 0; i < nrow; i++) row_set[i] = i;
    for (int j = 0; j < ncol; j++) col_set[j] = j;

    /* Run the algorithm */
    bimax_recurse(data_vec, nrow, ncol,
                  row_set, nrow, col_set, ncol,
                  minr, minc, max_biclusters,
                  result_rows, result_cols, num_found,
                  nrow, ncol);
}
