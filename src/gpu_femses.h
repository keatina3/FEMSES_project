#ifndef _FEMSES_H_
#define _FEMSES_H_

#ifdef __cplusplus

extern void gpu_femses(float *u, Mesh &M);

#endif

__device__ void calc_weights(float *we, float *cells, float *temp1, int idx, int idy);
__device__ void elems_glob_cpy(float *Le, float *be, float *temp1, int idx, int idy);

__global__ void assemble_elems_gpu(float *Le, float *be, float *we, float *vertices, int *cells, 
                        int *is_bound, float *bdry_vals);

#endif
