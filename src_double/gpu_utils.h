// =========================================================================== //
// Collection of utility functions for the GPU
// Includes linear solvers, card choosers and dot product
// =========================================================================== //

#ifndef _GPU_UTILS_H_

#ifdef __cpluscplus

extern void dummy(double *dat, int n);

#endif

#define _GPU_UTILS_H_

__global__ void dummy_kernel(int n);
__device__ double atomicAddDouble(double* address, double val);
__device__ double atomicExchDouble(double* address, double val);

void dnsspr_solve(double *L, double *b, int order, cudaEvent_t start, cudaEvent_t end, double &tau);
void sparse_solve(double *valsL, int *rowPtrL, int *colIndL, double *b, int order, int nnz);
void dense_solve(double *L, double *b, int order);
void error_dot_prod(double *a, double *b, int n, double &x);
void array_max(double *a, int n, int &max);

#endif
