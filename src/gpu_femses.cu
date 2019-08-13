#include <cmath>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <vector>
#include "mesh.h"
#include "utils.h"
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
    be[(idx*3) + idy] = be_shrd[idy];

    for(int i=0; i<3; i++){
        Le[(idx*9) + (idy*3) + i] = Le_shrd[(idy*3) + i];
    }
}

__device__ void elems_shared_cpy(float *Le, float *be, float *temp1, int idx, int idy){
    float *Le_shrd, *be_shrd;

    Le_shrd = temp1;
    be_shrd = &temp1[9];

    be_shrd[idy] = be[(idx*3) + idy];
    for(int i=0; i<3; i++){
        Le_shrd[(idy*3) + i] = Le[(idx*9) + (idy*3) + i];
    }
    __syncthreads();
    /* 
    if(idx == 1 && idy == 0){
        for(int i=0; i<3; i++){
            printf("femses: %f %f %f,  b[i] = %f\n", Le_shrd[i*3],Le_shrd[(3*i)+1],
                    Le_shrd[(3*i)+2], be_shrd[i]);
        }
    }
    */
}

__device__ void jacobi_iter(float *ue, float *temp1, int idx, int idy){
    float *Le_shrd, *be_shrd, ue_loc;

    Le_shrd = temp1;
    be_shrd = &temp1[9];
    // for Other iteration types, make ue_loc shared //
    
    /*
    if(idx == 0 && idy == 0){
        for(int i=0; i<3; i++){
            printf("femses: %f %f %f,  b[i] = %f\n", Le_shrd[i*3],Le_shrd[(3*i)+1],
                    Le_shrd[(3*i)+2], be_shrd[i]);
        }
    }
    */
    ue_loc = be_shrd[idy];
    for(int i=0; i<3; i++){
        if(i==idy) 
            continue;
        ue_loc -= Le_shrd[(idy*3) + i] * ue[(idx*3) + i];
    }
    ue_loc /= Le_shrd[(idy*3) + idy];

    __syncthreads();
    atomicExch(&ue[(idx*3) + idy], ue_loc);
    
    //ue[(idx*3) + idy] = ue_loc;
    //if(idx == 0 && idy == 0)
    //    printf("uloc = %f\n", ue[0]);
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
        
        /*
        if(idx == 0 && idy == 0){
            for(int i=0; i<3*gridDim.x; i++){
                printf("femses: %f %f %f,  b[i] = %f\n", Le[i*3],Le[(i*3)+1],Le[(i*3)+2],be[i]);
            }
        }
        */
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
        /*
        if(idx == 0 && idy == 1){
            for(int i=0; i<9; i++){
                printf("we[i] = %f\n", we[i]);
            }
            printf("Lii = %f\n", Lii);
        }
        */
        weight = Lii/we[v];

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
    float *Le, *be, *ue, *we;
    float *up_gpu, *un_gpu;
    float err = 1E16;

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
    cudaMalloc( (void**)&un_gpu, order*sizeof(float));
    cudaMalloc( (void**)&up_gpu, order*sizeof(float));
    cudaMalloc( (void**)&we, order*sizeof(float));
    
    cudaMemcpy(up_gpu, u, order*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(un_gpu, u, order*sizeof(float), cudaMemcpyHostToDevice);

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

    float *tmp;
    int count = 0;
    while(err > EPS && count < MAX_ITERS){
        // change shared memory accordingly //
        local_sols<<<dimGrid, dimBlock, shared*sizeof(float)>>>(Le, be, ue);
  
        cudaMemset(un_gpu, 0.0, order*sizeof(float)); 
        glob_sols<<<dimGrid, dimBlock, shared*sizeof(float)>>>(Le, we, un_gpu, ue, cells_gpu);

        error_dot_prod(un_gpu, up_gpu, order, err);
        std::cout << "error " << err <<std::endl;

        tmp = un_gpu;
        un_gpu = up_gpu;
        up_gpu = tmp;

        count++;
        std::cout << count << std::endl;
        if(count == MAX_ITERS){
            std::cout << "Maximum iterations reached!\n";
            exit(1);
        }
    }

    cudaMemcpy(u, up_gpu, order*sizeof(float), cudaMemcpyDeviceToHost);
    std::cout << "u[0] final = " << u[0] << std::endl;

    cudaFree(vertices_gpu); cudaFree(cells_gpu); cudaFree(dof_gpu);
    cudaFree(is_bound_gpu); cudaFree(bdry_vals_gpu);
    cudaFree(Le); cudaFree(be); cudaFree(un_gpu); cudaFree(up_gpu);
}
