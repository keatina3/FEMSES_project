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

extern void gpu_femses(double *u, Mesh &M, Tau &t, int &count, int &reconfig);

#endif

#define _FEMSES_H_

__device__ void calc_weights(double *w, double *cells, double *temp1, int idx, int idy);
__device__ void elems_glob_cpy(double *Le, double *be, double *temp1, int idx, int idy);
__device__ void elems_shared_cpy(double *Le, double *be, double *temp1, int idx, int idy);
__device__ void jacobi_iter(double *ue, double *up_glob, int *cells, 
                                double *temp1, int idx, int idy);

__global__ void assemble_elems_gpu(double *Le, double *be, double *w, double *ue, double *vertices, 
                                int *cells, int *is_bound, double *bdry_vals, 
                                int order, int num_cells);
__global__ void local_sols(double *Le, double *be, double *ue, 
                                double *up_glob, int *cells, int num_cells);
__global__ void glob_sols(double *Le, double *w, double *u_glob, 
                                double *ue, int *cells, int num_cells);

#endif
