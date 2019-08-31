// =========================================================================== //
// Functions to apply the FEM to solve a PDE on a GPU.                         //
// Uses decomposition of linear system and cuSparse/cuSolver libraries         //
// to solve the linear system.                                                 //
// Will solve dense/sparse systems. Will also assemble into dense and convert  //
// to CSR format and solve.                                                    //
// =========================================================================== //

#ifndef _GPU_FEM_H_

#ifdef __cplusplus

extern void gpu_fem(float *u, Mesh &M, Tau &t, int &reconfig);

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
                        int *is_bound, float *bdry_vals, int order, int num_cells, 
                        double *tau_d, int timing);
__global__ void assemble_gpu_csr(float *valsL, int *rowPtrL, int *colPtrL, float *b, 
                        float *vertices, int *cells, int *is_bound, float *bdry_vals, 
                        int order, int num_clles, double *tau_d, int timing);


#endif
