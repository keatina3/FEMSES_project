#ifndef _GPU_FEM_H_

#ifdef __cplusplus

extern void gpu_fem(Mesh &M);

#endif

#define _GPU_FEM_H_

__device__ float area(float *xi);
__device void elem_mat_gpu(float *vertices, float *cells, float *is_bound, float *bdry_vals, float *tmp1, int idx, int idy);
__global__ assemble_gpu(int num_cells, float *Le, float *be, float *vertices, float *cells, float *is_bound, float *bdry_vals);

#endif
