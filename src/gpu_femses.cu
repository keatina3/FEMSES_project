#include <iostream>
#include <cmath>
#include <cassert>
#include <cstdio>
#include <vector>
#include "mesh.h"
#include "utils.h"
#include "gpu_utils.h"
#include "gpu_fem.h"
#include "gpu_femses.h"

//////////// Calculates weighting for assembling single element solution ///////////
// One weight is evaluated for each node
// Added back to global memory
__device__ void calc_weights(float *w, int *cells, float *temp1, int idx, int idy){
    float *Le;
    int v;

    Le = temp1;
    v = cells[(idx*3) + idy];

    atomicAdd(&w[v], Le[(idy*3) + idy]);
}
////////


/////////// Copies element matrices/element vector from global-shared memory //////////
__device__ void elems_glob_cpy(float *Le, float *be, float *temp1, int idx, int idy){
    float *Le_shrd, *be_shrd;

    Le_shrd = temp1;
    be_shrd = &temp1[9];

    be[(idx*3) + idy] = be_shrd[idy];

    for(int i=0; i<3; i++){
        Le[(idx*9) + (idy*3) + i] = Le_shrd[(idy*3) + i];
    }
}
////////


/////////// Copies element matrices/element vector from shared-global memory //////////
__device__ void elems_shared_cpy(float *Le, float *be, float *temp1, int idx, int idy){
    float *Le_shrd, *be_shrd;

    Le_shrd = temp1;
    be_shrd = &temp1[9];

    be_shrd[idy] = be[(idx*3) + idy];
    for(int i=0; i<3; i++){
        Le_shrd[(idy*3) + i] = Le[(idx*9) + (idy*3) + i];
    }
    __syncthreads();
}
////////


////////////// Performs Jacobi iteration to get updated approximation of u ////////////
__device__ void jacobi_iter(
                float *ue,
                float *up_glob,
                float *w,
                int *cells,
                float *temp1,
                int idx,
                int idy)
{
    float *Le_shrd, *be_shrd;
    float ue_new, *ue_old;
    int v;
    float weight;

    Le_shrd = temp1;
    be_shrd = &temp1[9];
    ue_old  = &temp1[12];
    // v = cells[(idx*3) + idy];
    // weight = Le_shrd[(idy*3) + idy] / w[v]; 

    ue_new = be_shrd[idy];
    ue_old[idy] = ue[(idx*3) + idy];
    // ue_old[idy] = weight * up_glob[v];
    // ue_old[idy] = weight * up_glob[v];

    __syncthreads();
    
    for(int i=0; i<3; i++){
        if(i==idy)  continue;

        ue_new -= Le_shrd[(idy*3) + i] * ue_old[i];
    }
    ue_new /= Le_shrd[(idy*3) + idy];

    /*
    if(idx==2 && idy ==0){
        printf("u_new = %f\n", ue_new);
    }
    */
    __syncthreads();
    atomicExch(&ue[(idx*3) + idy], ue_new); // transferring element solution of u to global mem
}
//////


/////////////////// Kernel to assemble element solutions ///////////////////////////
// Element solutions are calculated in shared memory
// Element solutions are then transferred to an array in global memory
__global__ void assemble_elems_gpu(
                float *Le, 
                float *be, 
                float *w,
                float *ue,
                float *vertices, 
                int *cells, 
                int *is_bound, 
                float *bdry_vals)
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;      // idx = cell number
    int idy = blockIdx.y*blockDim.y + threadIdx.y;      // idy = local node number
    extern __shared__ float temp1[];

    if(idx < gridDim.x && idy < blockDim.y){
        // __device__ fn taken from other header to avoid code-reuse //
        assemble_elem(vertices, cells, is_bound, bdry_vals, temp1, idx, idy);
        calc_weights(w, cells, temp1, idx, idy);
        elems_glob_cpy(Le, be, temp1, idx, idy);
        ue[(idx*3) + idy] = 1.0;
    }
}
//////


/////////////// Kernel to calculate local approximation of solution ue /////////////
// Each cell has its own local solution for its element matrix and element vector
// These are apprimated with a jacobi iteration
__global__ void local_sols(
                float *Le,
                float *be,
                float *ue,
                float *up_glob,
                float *w,
                int *cells)
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;
    extern __shared__ float temp1[]; 

    /*
    if(idx ==0 && idy == 0){
        for(int e=0; e<gridDim.x; e++){
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    printf("%f ", Le[(e*9) + (i*3) + j]);
                }
                printf(" b = %f,\n",be[(e*3) + i]);
            }
            printf("\n");
        }
    }
    */

    if(idx < gridDim.x && idy < blockDim.y){
        elems_shared_cpy(Le, be, temp1, idx, idy);
        jacobi_iter(ue, up_glob, w, cells, temp1, idx, idy);
    }
   
    /* 
    if(idx == 2 && idy == 0){
        for(int i=0; i<3; i++){
            for(int j=0;j<3;j++){
                printf("%f ", temp1[(i*3) + j]);
            }
            printf("\n");
        }
    }
    */

}
///////


///////////// Kernel to calculate global approximation of u //////////////////
// Calculated by combining all local solutions ue with a weighting
__global__ void glob_sols(float *Le, float *w, float *u, float *ue, int *cells){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;
    int v;
    float Lii, weight;

    if(idx < gridDim.x && idy < blockDim.y){
        v = cells[(idx*3) + idy];               // getting global vertex number
        Lii = Le[(idx*9) + (idy*3) + idy];      
        
        weight = Lii/w[v];
        
        atomicAdd(&u[v], weight * ue[(idx*3) + idy]);
        /*
        if(idy==0 && idx ==0){
            for(int i=0; i<16; i++)
                printf("%f\n", w[i]);
        }
        */
    }
}
///////


/////////////// C++ function invoked to apply FEM-SES to solve PDE /////////////////////
// Applies the novel FEM - Single Element Solution approach to solve PDE
// Calculates element matrices as standard in regular approach
// Gets local solution approximations to these using a jacobi iteration
// Combines these to a global solution using a weighing
// Repeats until convergence of global solution
extern void gpu_femses(float *u, Mesh &M, Tau &t){
    int nr[2];
    int order, num_cells;
    int block_size_X, block_size_Y, shared;
    float *vertices_gpu, *vertices;
    int *cells_gpu, *cells;
    int *dof_gpu, *dof;
    int *is_bound_gpu, *is_bound;
    float *bdry_vals_gpu, *bdry_vals;
    float *Le, *be, *ue, *w;
    float *up_gpu, *un_gpu;
    float err = 1E16;
    cudaEvent_t start, finish;
    cudaError_t stat = cudaSuccess;
    float tau = 0.0;

    std::cout << GREEN "\nFEMSES Solver...\n" RESET;
    
    cudaSetDevice(k);
    
    cudaEventCreate(&start);
    cudaEventCreate(&finish);
    
    //////////////////////////// Gathering info from mesh /////////////////////////////

    M.get_recs(nr);

    order = (nr[0]+1)*(nr[1]+1);
    num_cells = 2*nr[0]*nr[1];
    M.get_arrays(&vertices, &cells, &dof, &is_bound, &bdry_vals);

    ///////////////////////////////////////////////////////////////////////////////////
    
    
    ////////////// Allocating memory for mesh/stiffnesss matrix/stress vector//////////
    ///////////  /array of element matrics/array of stress vectors/weighting //////////

    cudaEventRecord(start,0);

    stat = cudaMalloc( (void**)&vertices_gpu, 2*order*sizeof(float));
    assert(stat == cudaSuccess);
    stat = cudaMalloc( (void**)&cells_gpu, 3*num_cells*sizeof(int));
    assert(stat == cudaSuccess);
    stat = cudaMalloc( (void**)&dof_gpu, 3*num_cells*sizeof(int));
    assert(stat == cudaSuccess);
    stat = cudaMalloc( (void**)&is_bound_gpu, order*sizeof(int));
    assert(stat == cudaSuccess);
    stat = cudaMalloc( (void**)&bdry_vals_gpu, order*sizeof(float));
    assert(stat == cudaSuccess);

    stat = cudaMalloc( (void**)&Le, num_cells*3*3*sizeof(float));
    assert(stat == cudaSuccess);
    stat = cudaMalloc( (void**)&be, num_cells*3*sizeof(float));
    assert(stat == cudaSuccess);
    stat = cudaMalloc( (void**)&ue, num_cells*3*sizeof(float));
    assert(stat == cudaSuccess);
    stat = cudaMalloc( (void**)&un_gpu, order*sizeof(float));
    assert(stat == cudaSuccess);
    stat = cudaMalloc( (void**)&up_gpu, order*sizeof(float));
    assert(stat == cudaSuccess);
    stat = cudaMalloc( (void**)&w, order*sizeof(float));
    assert(stat == cudaSuccess);
    
    cudaEventRecord(finish);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&t.alloc, start, finish);

    ///////////////////////////////////////////////////////////////////////////////////


    ///////////////// Copying data for Mesh from host to device ///////////////////////

    std::cout << "      Copying data from host...\n";
    
    cudaEventRecord(start,0);

    stat = cudaMemcpy(vertices_gpu, vertices, 2*order*sizeof(float), cudaMemcpyHostToDevice);
    assert(stat == cudaSuccess);
    stat = cudaMemcpy(cells_gpu, cells, 3*num_cells*sizeof(int), cudaMemcpyHostToDevice);
    assert(stat == cudaSuccess);
    stat = cudaMemcpy(dof_gpu, dof, 3*num_cells*sizeof(int), cudaMemcpyHostToDevice);
    assert(stat == cudaSuccess);
    stat = cudaMemcpy(is_bound_gpu, is_bound, order*sizeof(int), cudaMemcpyHostToDevice);
    assert(stat == cudaSuccess);
    stat = cudaMemcpy(bdry_vals_gpu, bdry_vals, order*sizeof(float), cudaMemcpyHostToDevice);
    assert(stat == cudaSuccess);
    
    cudaEventRecord(finish);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&t.transfer, start, finish);

    stat = cudaMemset(up_gpu, 0, order*sizeof(float));
    assert(stat == cudaSuccess);

    ///////////////////////////////////////////////////////////////////////////////////


    //////////// DIMENSIONS OF SYSTEM => block per cell, 1 thread per node ////////////
    
    block_size_X = 1, block_size_Y = 3;
    dim3 dimBlock(block_size_X, block_size_Y);
    dim3 dimGrid((num_cells/dimBlock.x)+(!(num_cells%dimBlock.x)?0:1),
            (1/dimBlock.y)+(!(1%dimBlock.y)?0:1));

    shared = 31;

    //////////////////////////////////////////////////////////////////////////////////


    ////// Kernel to assemble element matrices and store in an array on glob mem /////

    std::cout << "      Getting element matrices...\n";
    
    cudaEventRecord(start,0);
    
    assemble_elems_gpu<<<dimGrid, dimBlock, shared*sizeof(float)>>>
                    (Le, be, w, ue, vertices_gpu, cells_gpu, is_bound_gpu, bdry_vals_gpu);

    cudaEventRecord(finish);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&t.elem_mats, start, finish);

    //////////////////////////////////////////////////////////////////////////////////


    ///////////////// Iterates through kernels until convergence /////////////////////

    std::cout << "      Applying Jacobi relaxation scheme...\n";
    
    cudaEventRecord(start,0);

    float *tmp;
    int count = 0;
    while(err > EPS && count < MAX_ITERS){
        // getting local solutions ue and storing on global mem //
        local_sols<<<dimGrid, dimBlock, shared*sizeof(float)>>>(Le, be, ue, up_gpu, w, cells);
        
        // setting un_gpu to 0 //
        cudaMemset(un_gpu, 0.0, order*sizeof(float));
        
        // assembling global solution estimate from weightings //
        glob_sols<<<dimGrid, dimBlock, shared*sizeof(float)>>>(Le, w, un_gpu, ue, cells_gpu);

        // calculating error using 2-norm //
        error_dot_prod(un_gpu, up_gpu, order, err);
        // std::cout << err << std::endl;

        tmp = un_gpu;
        un_gpu = up_gpu;
        up_gpu = tmp;

        count++;
        if(count == MAX_ITERS){
            std::cerr << "FEMSES - maximum iterations reached.\n";
            exit(1);
        }
    }

    cudaEventRecord(finish);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&t.solve, start, finish);

    std::cout << "      Solved in " << count << " iterations...\n";
    
    //////////////////////////////////////////////////////////////////////////////////

    
    //////////////// Tranferring soln to host from device & tidy /////////////////////

    std::cout << "      Transferring result back to host...\n";
    
    cudaEventRecord(start,0);
    
    cudaMemcpy(u, up_gpu, order*sizeof(float), cudaMemcpyDeviceToHost);
    
    cudaEventRecord(finish);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&tau, start, finish);
    t.transfer += tau;

    cudaFree(vertices_gpu);     cudaFree(cells_gpu);    cudaFree(dof_gpu);
    cudaFree(is_bound_gpu);     cudaFree(bdry_vals_gpu);
    cudaFree(Le); cudaFree(be); cudaFree(un_gpu);       cudaFree(up_gpu);

    //////////////////////////////////////////////////////////////////////////////////
}
////////
