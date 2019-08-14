// =========================================================================== //
// Functions to apply FEM Single Element Solition to solve a PDE on GPU
// Assembles element matrices and element vectors using standard approach
// Using a jacobi relaxation scheme, gets local solutions for each linear system
// Local solutions are then combines using a weighting to construct a global solution
// This is repeated until global solution converges
// Soluton then transferred back to host
// =========================================================================== //

#ifndef _FEMSES_H_

#ifdef __cplusplus

extern void gpu_femses(float *u, Mesh &M);

#endif

#define _FEMSES_H_

__device__ void calc_weights(float *we, float *cells, float *temp1, int idx, int idy);
__device__ void elems_glob_cpy(float *Le, float *be, float *temp1, int idx, int idy);
__device__ void elems_shared_cpy(float *Le, float *be, float *temp1, int idx, int idy);
__device__ void jacobi_iter(float *ue, float *temp1, int idx, int idy);

__global__ void assemble_elems_gpu(float *Le, float *be, float *we, float *vertices, int *cells, 
                        int *is_bound, float *bdry_vals);
__global__ void local_sols(float *Le, float *be, float *ue);
__global__ void glob_sols(float *Le, float *we, float *u, float *ue, int *cells);

#endif
