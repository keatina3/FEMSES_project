// =========================================================================== //
// Collection of utility functions for the GPU.                                //
// Includes linear solvers, card choosers and dot product.                     //
// =========================================================================== //

#ifndef _GPU_UTILS_H_

#ifdef __cpluscplus

extern void dummy(float *dat, int n);

#endif

#define _GPU_UTILS_H_

__global__ void dummy_kernel(int n);

void dnsspr_solve(float *L, float *b, int order, cudaEvent_t start, cudaEvent_t end, float &tau);
void sparse_solve(float *valsL, int *rowPtrL, int *colIndL, float *b, int order, int nnz);
void dense_solve(float *L, float *b, int order);
void error_dot_prod(float *a, float *b, int n, float &x);
void array_max(double *a, int n, int &max);

#endif
