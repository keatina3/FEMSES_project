#include <cassert>
#include <cstdio>
#include <iostream>
#include <vector>
#include "mesh.h"
#include "gpu_utils.h"
#include "gpu_fem.h"
#include "gpu_femses.h"

// need to change to DOF for more than P1 //
__device__ void calc_weights(float *we, int*cells, float *temp1, int idx, int idy){
    float *Le;
    int v;

    Le = temp1;
    v = cells[(idx*3) + idy];
    
    atomicAdd(&we[v], Le[(idy*3) + idy]);
}

__device__ void elems_glob_cpy(float *Le, float *be, float *temp1, int idx, int idy){
    float *Le_shrd, *be_shrd;

    Le_shrd = temp1;
    be_shrd = &temp1[9];

    // CHANGE THIS FROM BLOCK DIM IN CASE OF UNEVEN BLOCK SIZES //
    be[(idx*blockDim.y) + idy] = be_shrd[idy];

    for(int i=0; i<3; i++){
        Le[(idx*9) + (idy*3) + i] = Le_shrd[(idy*3) + i];
    }
}

__device__ void elems_shared_cpy(float *Le, float *be, float *temp1, int idx, int idy){
    float *Le_shrd, *be_shrd;
    
    Le_shrd = temp1;
    be_shrd = &temp1[9];
    
    be_shrd[idy] = be[(idx*blockDim.y) + idy];
    for(int i=0; i<3; i++){
        Le_shrd[(idy*3) + i] = Le[(idx*9) + (idy*3) + i];
    } 
}

__device__ void jacobi_iter(float *ue, float *temp1, int idx, int idy){
    float *Le_shrd, *be_shrd, ue_loc;
    
    Le_shrd = temp1;
    be_shrd = &temp1[9];
    // for Other iteration types, make ue_loc shared //

    ue_loc = be_shrd[idy];
    for(int i=0; i<3; i++){
        if(i==idy) 
            continue;
        ue_loc -= Le_shrd[(idy*3) + i] * ue[(idx*3) + i];
    }
    ue_loc /= Le_shrd[(idy*3) + idy];

    ue[(idx*3) + idy] = ue_loc;
}

__global__ void assemble_elems_gpu(float *Le, float *be, float *we, float *vertices, int *cells, int *is_bound, float *bdry_vals){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;
    extern __shared__ float temp1[];

    // check blockDim.y //
    if(idx < gridDim.x && idy < 3){
        // __device__ fn taken from other header to avoid code-reuse //
        assemble_elem(vertices, cells, is_bound, bdry_vals, temp1, idx, idy);
        calc_weights(we, cells, temp1, idx, idy);
        elems_glob_cpy(Le, be, temp1, idx, idy);
    }
}

__global__ void local_sols(float *Le, float *be, float *ue){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;
    extern __shared__ float temp1[]; 

    // CAN I NOT DO ATMOMIC ADD OF ELEMENT SOLUTION VALES AT EACH ITERATION ?? //
    if(idx < gridDim.x && idy < 3){
        elems_shared_cpy(Le, be, temp1, idx, idy);
        jacobi_iter(ue, temp1, idx, idy);
    }
        
}

__global__ void glob_sols(float *Le, float *we, float *u, float *ue, int *cells){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;
    int v;
    float Lii, weight;
    
    if(idx < gridDim.x && idy < 3){
        v = cells[(idx*3) + idy];
        Lii = Le[(idx*9) + (idy*3) + idy];
        weight = we[(idx*3) + idy];
        weight /= Lii;

        atomicAdd(&u[v], weight * ue[(idx*3) + idy]);
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
    float *Le, *be, *ue, *u_gpu, *we;

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
    cudaMalloc( (void**)&ue, num_cells*dim*sizeof(float));
    cudaMalloc( (void**)&u_gpu, order*sizeof(float));
    cudaMalloc( (void**)&we, order*sizeof(float));

    block_size_X = 1, block_size_Y = 3;
    // check these dimensions //
    // this will be scope for efficency experimentation //
    dim3 dimBlock(block_size_X, block_size_Y);
    dim3 dimGrid((num_cells/dimBlock.x)+(!(num_cells%dimBlock.x)?0:1),
            (1/dimBlock.y)+(!(1%dimBlock.y)?0:1));

    int dofs = 3;
    shared = (dofs*dofs) + dofs + 1;
    // MIGHT LEAVE SHARED AS 31 UNTIL LATER //
    if(true)
        shared += 18;

    assemble_elems_gpu<<<dimGrid, dimBlock, shared*sizeof(float)>>>(Le, be, we, vertices_gpu, cells_gpu, is_bound_gpu, bdry_vals_gpu);

    // change shared memory accordingly //
    local_sols<<<dimGrid, dimBlock, shared*sizeof(float)>>>(Le, be, ue);
    
    glob_sols<<<dimGrid, dimBlock, shared*sizeof(float)>>>(Le, we, u, ue, cells_gpu); 

    ////////// REST OF FEMSES GOES HERE //////////

    cudaFree(vertices_gpu); cudaFree(cells_gpu); cudaFree(dof_gpu);
    cudaFree(is_bound_gpu); cudaFree(bdry_vals_gpu);
    cudaFree(Le); cudaFree(be); cudaFree(u_gpu);

    std::cout << "test7\n";

}
