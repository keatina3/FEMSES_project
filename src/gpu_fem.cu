#include <cassert>
#include <cstdio>
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

    __syncthreads();
    
    if(idy==0)
        consts[9] = area(xi);

    // alpha, beta, gamma //
    // consts[(3*idy)] = xi[(idy+1)%3][1] * xi[(i+2)%3][2] - xi[(i+2)%3][1] * xi[(i+1)%3][2];
    consts[(3*idy)+1] = xi[ 3*((idy+1)%3) +2] - xi[ 3*((idy+2)%3) + 2];
    consts[(3*idy)+2] = xi[ 3*((idy+2)%3) +1] - xi[ 3*((idy+1)%3) + 1];
    
    be[idy] = 0.0;      // or Int(fv) //

    __syncthreads();
    /*
    if(idx==2 && idy==2){
        printf("b0 = %f, b1 = %f, b2 = %f\n", consts[1], consts[4], consts[7]);
        printf("g0 = %f, g1 = %f, g2 = %f\n", consts[2], consts[5], consts[8]);
        printf("del = %f\n",consts[9]);
    }    
    
    if(idx == 0 && idy == 2){
        for(int i=0;i<9;i++){
            printf("x = %f, y = %f, bdry = %d\n", vertices[2*i], vertices[(2*i)+1], is_bound[i]);
        }
        for(int i=0;i<8;i++){
            printf("cell = %d, cells[0] = %d , cells[1] = %d, cells[2] = %d\n", i, cells[i*3], cells[(i*3)+1], cells[(i*3)+2]);
        }
    }
    */
    for(int i=0; i<=idy; i++){
        Le[(3*idy)+i] = 0.25*consts[9]*consts[9]*consts[9] * (consts[(3*idy)+1]*consts[(3*i)+1] 
                                    + consts[(3*idy)+2]*consts[(3*i)+2]);
        Le[(3*i)+idy] = Le[(3*idy)+i];
    }
    __syncthreads();

    if(is_bound[v]){
        bound = bdry_vals[v];
        // change this appropriately for more DOF if necessary //
        for(int j=0; j<3; j++){
            if(idy != j){
                atomicAdd(&be[j], (-1)*Le[(3*j) + idy]*bound);
                Le[(3*j) + idy] = 0.0;
                Le[(3*idy) + j] = 0.0;
            }
        }
        __syncthreads();
        Le[(3*idy)+idy] = 1.0;
        atomicExch(&be[idy], bound);
    }                            
    /* 
    if(idx == 2 && idy == 2){
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                printf("%f ", Le[(3*i) + j]);
            }
            printf(" bi = %f\n",be[i]);
        }
    }
    */
}

__device__ void assemble_mat(float *L, float *b, float *vertices, int *dof, float *temp1, int idx, int idy, int order){
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
    
    /*
    if(idx == 0 && idy == 0){
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                printf("%f ", Le[(3*i) + j]);
            }
            printf(" bi = %f\n",be[i]);
        }
        printf("dof_0 = %d, dof_1 = %d, dof_2 = %d\n",dof_r[0], dof_r[1], dof_r[2]);
        printf("global b[2] = %f\n", b[2]);
        printf("L[0][0] = %f\n",L[0]);
    }
    */
    
    // ANY ALTERNATIVE TO ATOMICS ??? //
    // b[dof[idy]] += be[idy];
    // printf("dof_r = %d\n",dof_r[idy]);
    atomicAdd(&b[dof_r[idy]], be[idy]);
    
    for(int i=0; i<=idy; i++){
        // L[dof_r[idy]][dof_r[(idy+i)%3]] += Le[idy][(idy+i)%3];
        // DEFINITELY CHECK THIS //
        atomicAdd(&L[(order*dof_r[idy]) + dof_r[i]], 
                            Le[(3*idy) + i]);
        if(i==idy)
            continue;
        else {
            atomicAdd(&L[(order*dof_r[i]) + dof_r[idy]], 
                            Le[(3*idy) + i]);
        }
    }
}

__global__ void assemble_gpu(float *L, float *b, float *vertices, int *cells, int *is_bound, float *bdry_vals, int order){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;
    extern __shared__ float temp1[];

    if(idx < gridDim.x && idy < 3){
        elem_mat_gpu(vertices, cells, is_bound, bdry_vals, temp1, idx, idy);
        assemble_mat(L, b, vertices, cells, temp1, idx, idy, order);
        
        /* 
        if(idx == 2 && idy == 2){
            printf("\n\n");
            for(int i=0;i<order;i++){
                for(int j=0;j<order;j++){
                    printf("%f ", L[(i*order)+j]);
                }
                printf("bi = %f\n", b[i]);
            }
        }
        */
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
    float *Workspace;
    int Lwork, devInfo;
     
    M.get_recs(nr);

    num_nodes = (nr[0]+1)*(nr[1]+1);
    dim = 2+1;      // needs expansion here //
    order = num_nodes + 0;
    num_cells = 2*nr[0]*nr[1];
    n = lda = order;
    // Sorting Mesh // 
    M.get_arrays(&vertices, &cells, &dof, &is_bound, &bdry_vals);

    std::cout << num_nodes << std::endl;
    
    std::cout << "test1\n";
    
    cudaMalloc( (void**)&vertices_gpu, 2*num_nodes*sizeof(float));
    cudaMalloc( (void**)&cells_gpu, dim*num_cells*sizeof(int));
    cudaMalloc( (void**)&dof_gpu, dim*num_cells*sizeof(int));
    cudaMalloc( (void**)&is_bound_gpu, num_nodes*sizeof(int));
    cudaMalloc( (void**)&bdry_vals_gpu, num_nodes*sizeof(float));
 
    cudaMemcpy(vertices_gpu, vertices, 2*num_nodes*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(cells_gpu, cells, dim*num_cells*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dof_gpu, dof, dim*num_cells*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(is_bound_gpu, is_bound, num_nodes*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(bdry_vals_gpu, bdry_vals, num_nodes*sizeof(float), cudaMemcpyHostToDevice);
    //////////////////

    std::cout << "test2\n";

    cudaMalloc( (void**)&L, order*order*sizeof(float));
    cudaMalloc( (void**)&b, order*sizeof(float));
    
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
    dim3 dimGrid((num_cells/dimBlock.x)+(!(order%dimBlock.x)?0:1),
                (1/dimBlock.y)+(!(1%dimBlock.y)?0:1));
    
    // shared memory will grow for more DOF //
    // this assumes 1 dof per triangle //
    assemble_gpu<<<dimGrid, dimBlock, 31*sizeof(float)>>>(L, b, vertices_gpu, cells_gpu, is_bound_gpu, bdry_vals_gpu, order);
    
    std::cout << "test4\n";
    
    status = cusolverDnSpotrf_bufferSize(handle, uplo, order, L, order, &Lwork);
    assert(CUSOLVER_STATUS_SUCCESS == status);
    cudaMalloc( (void**)&Workspace, Lwork*sizeof(float));
    std::cout << "bufferSize " << Lwork << std::endl;

    status = cusolverDnSpotrf(handle, uplo, order, L, order, Workspace, Lwork, &devInfo);
    cudaStat1 = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == status);
    assert(cudaSuccess == cudaStat1);
    
    status = cusolverDnSpotrs(handle, uplo, order, nrhs, L, order, b, order, &devInfo);
    cudaStat1 = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == status);
    assert(cudaSuccess == cudaStat1);

    std::cout << "test5\n";

    cudaMemcpy(u, b, order*sizeof(float), cudaMemcpyDeviceToHost);
    
    std::cout << "test6\n";
    
    cusolverDnDestroy(handle);

    cudaFree(vertices_gpu); cudaFree(cells_gpu); cudaFree(dof_gpu);
    cudaFree(is_bound_gpu); cudaFree(bdry_vals_gpu);
    cudaFree(L); cudaFree(b);
    
    std::cout << "test7\n";
    
}
