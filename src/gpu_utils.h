// =========================================================================== //
// Collection of utility functions for the GPU
// Includes linear solvers, card choosers and dot product
// =========================================================================== //

#ifndef _GPU_UTILS_H_
#define _GPU_UTILS_H_

void dnsspr_solve(float *L, float *b, int order);
void sparse_solve(float *valsL, int *rowPtrL, int *colIndL, float *b, int order, int nnz);
void dense_solve(float *L, float *b, int order);
void error_dot_prod(float *a, float *b, int n, float &x);

#endif
