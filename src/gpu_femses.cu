#include <cassert>
#include <cstdio>
#include <iostream>
#include <vector>
#include "mesh.h"
#include "gpu_utils.h"
#include "gpu_fem.h"
#include "gpu_femses.h"

__device__ void elems_glob_cpy(float *Le, float *be, float *temp1, int idx, int idy){
    float *Le_shrd, *be_shrd;

    Le_shrd = temp1;
    be_shrd = &temp1[9];

    be[(idx*blockDim.y) + idy] = be_shrd[idy];
    for(int i=0; i<blockDim.y; i++){
        Le[(idx*blockDim.y*blockDim.y) + (idy*blockDim.y) + i] = Le_shrd[(idy*blockDim.y) + i];
    }
}

// need to change to DOF for more than P1 //
__device__ void calc_weights(float *we, int*cells, float *temp1, int idx, int idy){
    float *Le;
    int v;

    Le = temp1;
    v = cells[(idx*blockIdx.y) + idy];

    atomicAdd(&we[v], Le[(idy*blockIdx.y) + idy]);
}

__global__ void assemble_elems_gpu(float *Le, float *be, float *we, float *vertices, int *cells, int *is_bound, float *bdry_vals){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;
    extern __shared__ float temp1[];

    // check blockDim.y //
    if(idx < gridDim.x && idy < blockDim.y){
        // __device__ fn taken from other header to avoid code-reuse //
        assemble_elem(vertices, cells, is_bound, bdry_vals, temp1, idx, idy);
        calc_weights(we, cells, temp1, idx, idy);
        elems_glob_cpy(Le, be, temp1, idx, idy);
    }
}

extern void gpu_femses(float *u, Mesh &M){
    int nr[2];
    int num_nodes, dim, order, num_cells;
    int block_size_X, block_size_Y, shared;
    float *vertices_gpu, *vertices;
    int *cells_gpu, *cells;
    int *dof_gpu, *dof;
    int *is_bound_gpu, *is_bound;
    float *bdry_vals_gpu, *bdry_vals;
    float *Le, *be, *u_gpu, *we;

    M.get_recs(nr);

    num_nodes = (nr[0]+1)*(nr[1]+1);
    dim = 2+1;      // needs expansion here //
    order = num_nodes + 0;
    num_cells = 2*nr[0]*nr[1];
    M.get_arrays(&vertices, &cells, &dof, &is_bound, &bdry_vals);

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

    cudaMalloc( (void**)&Le, num_cells*dim*dim*sizeof(float));
    cudaMalloc( (void**)&be, num_cells*dim*sizeof(float));
    cudaMalloc( (void**)&u_gpu, order*sizeof(float));
    cudaMalloc( (void**)&we, order*sizeof(float));

    block_size_X = 1, block_size_Y = 3;
    // check these dimensions //
    // this will be scope for efficency experimentation //
    dim3 dimBlock(block_size_X, block_size_Y);
    dim3 dimGrid((num_cells/dimBlock.x)+(!(order%dimBlock.x)?0:1),
            (1/dimBlock.y)+(!(1%dimBlock.y)?0:1));

    int dofs = 3;
    shared = (dofs*dofs) + dofs + 1;
    // MIGHT LEAVE SHARED AS 31 UNTIL LATER //
    if(true)
        shared += 18;

    assemble_elems_gpu<<<dimGrid, dimBlock, shared*sizeof(float)>>>(Le, be, we, vertices_gpu, cells_gpu, is_bound_gpu, bdry_vals_gpu);
     

    ////////// REST OF FEMSES GOES HERE //////////

    cudaFree(vertices_gpu); cudaFree(cells_gpu); cudaFree(dof_gpu);
    cudaFree(is_bound_gpu); cudaFree(bdry_vals_gpu);
    cudaFree(Le); cudaFree(be); cudaFree(u_gpu);

    std::cout << "test7\n";

}
