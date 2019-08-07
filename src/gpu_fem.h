#ifndef _GPU_FEM_H_

#ifdef __cplusplus

extern void gpu_fem(float *u, Mesh &M);

#endif

#define _GPU_FEM_H_

__device__ float area(float *xi);
__device__ void assemble_elem(float *vertices, int *cells, int *is_bound, float *bdry_vals, 
                        float *tmp1, int idx, int idy);
__device__ void assemble_mat(float *L, float *b, float *vertices, int *dof, float *temp1, 
                        int idx, int idy, int order);
__device__ void assemble_mat_csr(float *valsL, int *rowPtrL, int *colPtrL, float *b, 
                        float *vertices, int *dof, float *temp1, int idx, int idy, int order);


__global__ void assemble_gpu(float *L, float *b, float *vertices, int *cells, 
                        int *is_bound, float *bdry_vals, int order);
__global__ void assemble_gpu_csr(float *valsL, int *rowPtrL, int *colPtrL, float *b, 
                        float *vertices, int *cells, int *is_bound, float *bdry_vals, int order);


#endif
