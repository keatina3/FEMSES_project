#include <cassert>
#include <iostream>
#include <vector>
#include <cuda_runtime.h>
#include <cusolverDn.h>     // need to change this for sparse //
#include "mesh.h"
#include "gpu_fem.h"

__device__ float area(float *xi){
    float tmp = 0.0;

    tmp += xi[0]*(xi[4]*xi[8] - xi[5]*xi[7]);
    tmp -= xi[1]*(xi[3]*xi[8] - xi[5]*xi[6]);
    tmp += xi[2]*(xi[3]*xi[7] - xi[4]*xi[6]);

    return 0.5*tmp;
}

__device__ void elem_mat_gpu(float *vertices, int *cells,  int *is_bound, float *bdry_vals, float *temp1, int idx, int idy){
    float *Le, *be, *xi, *consts;
    int v;
    float bound;

    Le = temp1;
    be = &temp1[9];
    xi = &temp1[12];
    consts = &temp1[21];
    
    // Potentially write function to return global vertex index //
    v = cells[(idx*3) + idy];
    
    xi[3*idy] = 1.0;
    xi[(3*idy) + 1] = vertices[2*v]; 
    xi[(3*idy) + 2] = vertices[(2*v)+1];

    if(idy==0)
        consts[9] = area(xi);
    __syncthreads();

    // alpha, beta, gamma //
    // consts[(3*idy)] = xi[(idy+1)%3][1] * xi[(i+2)%3][2] - xi[(i+2)%3][1] * xi[(i+1)%3][2];
    consts[(3*idy)+1] = xi[ 3*((idy+1)%3) +2] - xi[ 3*((idy+2)%3) + 2];
    consts[(3*idy)+2] = xi[ 3*((idy+2)%3) +1] - xi[ 3*((idy+1)%3) + 1];
    
    be[idy] = 0.0;      // or Int(fv) //
     
    for(int i=0; i<idy; i++){
        Le[(3*idy)+i] = 0.25*consts[9]*consts[9]*consts[9] * (consts[(3*idy)+1]*consts[(3*i)+1] 
                                    + consts[(3*idy)+2]*consts[(3*i)+2]);
        Le[(3*i)+idy] = Le[(3*idy)+i];
    }

    if(is_bound[v]){
        bound = bdry_vals[v];
        // change this appropriately for more DOF if necessary //
        for(int j=0; j<3; j++){
            if(idy != j){
                be[j] -= Le[(3*j) + idy]*bound;
                Le[(3*j) + idy] = 0.0;
                Le[(3*idy) + j] = 0.0;
            }
        }
        Le[(3*idy)+idy] = 1.0;
        be[idy] = bound;
    }                            
}

// CHANGE VERTICES ETC FROM FLOATS //
// WILL NEED TO CHANGE THIS FROM "CELLS" to "DOF" WHEN EXPANDING TO PN //
__device__ void assemble_mat(float *L, float *b, float *vertices, int *dof, float *temp1, int idx, int idy){
    float *Le, *be;
    int* dof_r;
 
    // size of these will change for non P1 //   
    Le = temp1;
    be = &temp1[9];
    dof_r = (int *)&temp1[12];

    // add in for loop for DOF // 
    if(idy==0){
        dof_r[0] = dof[(idx*3)];
        dof_r[1] = dof[(idx*3)+1];
        dof_r[2] = dof[(idx*3)+2];
    }
    __syncthreads();
    
    // ANY ALTERNATIVE TO ATOMICS ??? //
    // b[dof[idy]] += be[idy];
    atomicAdd(&b[dof_r[idy]], be[idy]);
     
    for(int i=0; i<3; i++){
        // L[dof_r[idy]][dof_r[(idy+i)%3]] += Le[idy][(idy+i)%3];
        // DEFINTELY CHECK THIS //
        atomicAdd(&L[((gridDim.x)*dof_r[idy]) + dof_r[(idy+i)%3]], 
                            Le[(3*idy) + (idy+i)%3]);
    }
}

__global__ void assemble_gpu(int num_cells, float *L, float *b, float *vertices, int *cells, int *is_bound, float *bdry_vals){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;
    extern __shared__ float temp1[];

    if(idx < num_cells && idy < 3){
        elem_mat_gpu(vertices, cells, is_bound, bdry_vals, temp1, idx, idy);
        assemble_mat(L, b, vertices, cells, temp1, idx, idy);
    }
} 

extern void gpu_fem(float *u, Mesh &M){
    int nr[2];
    int num_nodes, dim, order, num_cells;
    int block_size_X, block_size_Y;
    float *vertices_gpu, *vertices;
    int *cells_gpu, *cells;
    int *dof_gpu, *dof;
    int *is_bound_gpu, *is_bound;
    float *bdry_vals_gpu, *bdry_vals;
    float *L, *b;
    
    cusolverDnHandle_t handle = NULL;
    cudaStream_t stream = NULL;
    cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
    cudaError_t cudaStat1 = cudaSuccess;
    const cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
    const int nrhs = 1;
    int n, lda;
    float Lwork;
    int devInfo;
     
    M.get_recs(nr);

    num_nodes = (nr[0]+1)*(nr[1]+1);
    dim = 2+1;      // needs expansion here //
    order = num_nodes + 0;
    num_cells = 2*nr[0]*nr[1];
    n = lda = num_cells;
    // Sorting Mesh // 
    M.get_arrays(&vertices, &cells, &dof, &is_bound, &bdry_vals);

    std::cout << num_nodes << std::endl;
    
    std::cout << "test1\n";
    
    cudaMalloc( (void**)&vertices_gpu, 2*num_nodes*sizeof(float));
    cudaMalloc( (void**)&cells_gpu, num_cells*sizeof(int));
    cudaMalloc( (void**)&dof_gpu, num_cells*sizeof(int));
    cudaMalloc( (void**)&is_bound_gpu, num_nodes*sizeof(int));
    cudaMalloc( (void**)&bdry_vals_gpu, num_nodes*sizeof(float));
 
    cudaMemcpy(vertices_gpu, vertices, 2*num_nodes*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(cells_gpu, cells, num_cells*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dof_gpu, dof, num_cells*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(is_bound_gpu, is_bound, num_nodes*sizeof(bool), cudaMemcpyHostToDevice);
    cudaMemcpy(bdry_vals_gpu, bdry_vals, num_nodes*sizeof(float), cudaMemcpyHostToDevice);
    //////////////////

    std::cout << "test2\n";

    cudaMalloc( (void**)&L, order*order*sizeof(float));
    // cudaMalloc( (void**)&L0, order*order*sizeof(float));
    cudaMalloc( (void**)&b, order*sizeof(float));
    //cudaMalloc(b);
    //cudaMalloc(be);
    
    std::cout << "test3\n";

    status = cusolverDnCreate(&handle);
    assert(CUSOLVER_STATUS_SUCCESS == status);
    
    cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    assert(cudaSuccess == cudaStat1);
    
    status = cusolverDnSetStream(handle, stream);
    assert(CUSOLVER_STATUS_SUCCESS == status); 
    
    block_size_X = 1, block_size_Y = 3;
    // check these dimensions //
    // this will be scope for efficency experimentation //
    dim3 dimBlock(block_size_X, block_size_Y);
    dim3 dimGrid((order/dimBlock.x)+(!(order%dimBlock.x)?0:1),
                (1/dimBlock.y)+(!(1%dimBlock.y)?0:1));
    
    // shared memory will grow for more DOF //
    // this assumes 1 dof per triangle //
    assemble_gpu<<<dimGrid, dimBlock, 31*sizeof(float)>>>(order, L, b, vertices_gpu, cells_gpu, is_bound_gpu, bdry_vals_gpu);
    
    std::cout << "test4\n";

    status = cusolverDnSpotrf(handle, uplo, n, L, lda, &Lwork, 1, &devInfo);
    cudaStat1 = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == status);
    // assert(cudaSuccess == cudaStat1);
    
    status = cusolverDnSpotrs(handle, uplo, n, nrhs, L, lda, b, lda, &devInfo);
    cudaStat1 = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == status);
    assert(cudaSuccess == cudaStat1);

    std::cout << "test5\n";

    cudaMemcpy(u, b, order*sizeof(float), cudaMemcpyDeviceToHost);
    
    std::cout << "test6\n";
    
    cudaFree(vertices_gpu); cudaFree(cells_gpu); cudaFree(dof_gpu);
    cudaFree(is_bound_gpu); cudaFree(bdry_vals_gpu);
    cudaFree(L); cudaFree(b);
    
    std::cout << "test7\n";
    
}
