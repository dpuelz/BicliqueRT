#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern void bimax_c(int *data_vec, int *nrow_r, int *ncol_r,
                    int *minr_r, int *minc_r, int *number_r,
                    int *result_rows, int *result_cols, int *num_found);

static const R_CMethodDef CEntries[] = {
    {"bimax_c", (DL_FUNC) &bimax_c, 9},
    {NULL, NULL, 0}
};

void R_init_BicliqueRT(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
