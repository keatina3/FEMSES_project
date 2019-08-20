#include <iostream>
#include <cassert>
#include <vector>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cusolverSp.h>
#include <cusparse.h>
#include "mesh.h"
#include "utils.h"
#include "gpu_utils.h"
#include "gpu_fem.h"

//////////// Calculates area of triangle, given coordinates ////////////////
__device__ float area(float *xi){
    float tmp = 0.0;

    tmp += xi[0]*(xi[4]*xi[8] - xi[5]*xi[7]);
    tmp -= xi[1]*(xi[3]*xi[8] - xi[5]*xi[6]);
    tmp += xi[2]*(xi[3]*xi[7] - xi[4]*xi[6]);

    return 0.5*tmp;
}
///////


///////////// Assembles element matrix and stress vector /////////////////////
// Calculated and stored on shared memory
__device__ void assemble_elem(
                float *vertices, 
                int *cells,  
                int *is_bound, 
                float *bdry_vals, 
                float *temp1, 
                int idx, 
                int idy)
{
    float *Le, *be, *xi, *consts;
    int v;
    float bound;

    Le = temp1;                 // element matrix
    be = &temp1[9];             // element stress vector
    xi = &temp1[12];            // matrix of global coordinates
    consts = &temp1[21];        // stores alpha, beta, gamma seen in serial version
    
    v = cells[(idx*3) + idy];   // global node number 
    
    //////////// Assigning global coordinates //////////////////
    
    xi[3*idy] = 1.0;
    xi[(3*idy) + 1] = vertices[2*v]; 
    xi[(3*idy) + 2] = vertices[(2*v)+1];

    __syncthreads();
    
    // using 1 thread to calculate area //
    if(idy==0)
        consts[9] = area(xi);
    
    ///////////////////////////////////////////////////////////


    //////////// calculating alpha, beta, gamma ///////////////
    
    // consts[(3*idy)] = xi[(idy+1)%3][1] * xi[(i+2)%3][2] - xi[(i+2)%3][1] * xi[(i+1)%3][2];
    consts[(3*idy)+1] = xi[ 3*((idy+1)%3) +2] - xi[ 3*((idy+2)%3) + 2];
    consts[(3*idy)+2] = xi[ 3*((idy+2)%3) +1] - xi[ 3*((idy+1)%3) + 1];
    
    __syncthreads();
    
    ////////////////////////////////////////////////////////////


    ///////////////////// Calculating LHS & RHS //////////////////

    be[idy] = 0.0;      // 0.0 intially until BCs enforced //
    for(int i=0; i<=idy; i++){
        Le[(3*idy)+i] = 0.25*consts[9]*consts[9]*consts[9] * (consts[(3*idy)+1]*consts[(3*i)+1] 
                                    + consts[(3*idy)+2]*consts[(3*i)+2]);
        Le[(3*i)+idy] = Le[(3*idy)+i];
    }
    __syncthreads();

    ///////////////////////////////////////////////////////////////


    ///////////////// Enforcing boundary conditions ///////////////
    
    if(is_bound[v]){
        bound = bdry_vals[v];
        // change this appropriately for more DOF if necessary //
        for(int j=0; j<3; j++){
            if(idy != j){
                atomicAdd(&be[j], (-1)*Le[(3*j) + idy]*bound);
                atomicExch(&Le[(3*j) + idy],  0.0);
                atomicExch(&Le[(3*idy) + j], 0.0);
            }
        }
        __syncthreads();
        Le[(3*idy)+idy] = 1.0;
        atomicExch(&be[idy], bound);
    }                            
    __syncthreads();
    
    /////////////////////////////////////////////////////////////////
}
////////


///////////////////// Assembles stiffness matrix and stress vector //////////
// Takes shared memory element matrices and vectors
// and maps them back to global memory as an assembled stiffness matrix
// Stored in dense format
__device__ void assemble_mat(
                float *L, 
                float *b, 
                float *vertices, 
                int *dof, 
                float *temp1, 
                int idx, 
                int idy, 
                int order)
{
    float *Le, *be;
    int* dof_r;
 
    Le = temp1;
    be = &temp1[9];
    dof_r = (int *)&temp1[12];      // stores in shared memory, global node numbers for 3 nodes

    ///////////////// Assigning global node numbers //////////////////

    if(idy==0){
        dof_r[0] = dof[(idx*3)];
        dof_r[1] = dof[(idx*3)+1];
        dof_r[2] = dof[(idx*3)+2];
    }
    __syncthreads();

    ////////////////////////////////////////////////////////////////
  

    ///////// Mapping values back to global stiffness matrix /////// 
    
    atomicAdd(&b[dof_r[idy]], be[idy]);
    
    for(int i=0; i<=idy; i++){
        atomicAdd(&L[(order*dof_r[idy]) + dof_r[i]], 
                            Le[(3*idy) + i]);
        if(i==idy)      // avoid double adding diagonal element
            continue;
        else {
            atomicAdd(&L[(order*dof_r[i]) + dof_r[idy]], 
                            Le[(3*idy) + i]);
        }
    }

    /////////////////////////////////////////////////////////////////
}
///////



///////////////////// Assembles stiffness matrix and stress vector //////////
// Takes shared memory element matrices and vectors
// and maps them back to global memory as an assembled stiffness matrix
// Stored in CSR format
__device__ void assemble_mat_csr(
                float *valsL, 
                int *rowPtrL, 
                int *colIndL, 
                float *b, 
                float *vertices, 
                int *dof, 
                float *temp1, 
                int idx, 
                int idy, 
                int order)
{
    float *Le, *be;
    int* dof_r;
    int row;
    int *tmp1, *tmp2;
    int off = 0;

    Le = temp1;
    be = &temp1[9];
    dof_r = (int *)&temp1[12];

    ///////////////// Assigning global node numbers //////////////////
    
    if(idy==0){
        dof_r[0] = dof[(idx*3)];
        dof_r[1] = dof[(idx*3)+1];
        dof_r[2] = dof[(idx*3)+2];
    }
    __syncthreads();
    
    //////////////////////////////////////////////////////////////////


    ///////// Mapping values back to global stiffness matrix ////////
    
    atomicAdd(&b[dof_r[idy]], be[idy]);
    
    row = rowPtrL[dof_r[idy]];
    tmp1 = &colIndL[row];
    for(int i=0; i<3; i++){
        tmp2 = tmp1;
        while(*tmp2 != dof_r[i]){
            off++;
            tmp2++;
        }
        atomicAdd(&valsL[row + off], Le[(3*idy) + i]);
        off = 0;
    }
    __syncthreads();

    ///////////////////////////////////////////////////////////////////
}
/////////


////////////// Kernel to calculate elements and assemble global stiffness matrix /////////////
// Assembled in dense format
__global__ void assemble_gpu(
                float *L, 
                float *b, 
                float *vertices, 
                int *cells, 
                int *is_bound, 
                float *bdry_vals, 
                int order
                float *tau_d)
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;      // idx = cell number
    int idy = blockIdx.y*blockDim.y + threadIdx.y;      // idy = local node number
    extern __shared__ float temp1[];            // shared mem to store elem mats & constants

    if(idx < gridDim.x && idy < blockDim.y){
        
        clock_t start = clock();
        assemble_elem(vertices, cells, is_bound, bdry_vals, temp1, idx, idy);
        clock_t end = clock();
        tau_d[(idx*blockDim.y) + idy] = (float)(stop - start);

        clock_t start = clock();
        assemble_mat(L, b, vertices, cells, temp1, idx, idy, order);
        clock_t end = clock();
        tau_d[(gridDim.x * blockDim.y) + (idx*blockDim.y) + idy] = (float)(stop - start);
    }
}
///////



////////////// Kernel to calculate elements and assemble global stiffness matrix /////////////
// Assembled in CSR format
// see comments above
__global__ void assemble_gpu_csr(
                float *valsL, 
                int *rowPtrL, 
                int *colIndL, 
                float *b, 
                float *vertices, 
                int *cells, 
                int *is_bound, 
                float *bdry_vals, 
                int order
                float *tau_d)
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;
    extern __shared__ float temp1[];

    if(idx < gridDim.x && idy < blockDim.y){
        assemble_elem(vertices, cells, is_bound, bdry_vals, temp1, idx, idy);
        assemble_mat_csr(valsL, rowPtrL, colIndL, b, vertices, cells, temp1, idx, idy, order);
    }
}    
////////


///////////////// C++ Function to be invoked from host to apply FEM to PDE ///////////////////
// Applies standard approach of assmebling stiffness matrix and
// decomposing the linear system on the GPU to solve
extern void gpu_fem(float *u, Mesh &M, Tau &t){
    int nr[2];
    int order, num_cells, nnz;
    int block_size_X, block_size_Y, shared;
    float *vertices_gpu, *vertices;
    int *cells_gpu, *cells;
    int *dof_gpu, *dof;
    int *is_bound_gpu, *is_bound;
    float *bdry_vals_gpu, *bdry_vals;
    float *L, *b, *valsL;
    int *rowPtrL, *colIndL;
    std::vector<float> valsLCPU;
    std::vector<int> rowPtrLCPU;
    std::vector<int> colIndLCPU;
    cudaError_t stat;
    cudaEvent_t start, finish;
    float tau = 0.0;
    float *tau_d;
    float tau_d_max[2];

    std::cout << GREEN "\nGPU Solver...\n" RESET;
    
    setCudaDevice(k);
    
    cudaEventCreate(&start);
    cudaEventCreate(&finish);

    //////////////////////// Gathering info from mesh ////////////////////////////////

    M.get_recs(nr);

    order = (nr[0]+1)*(nr[1]+1);
    num_cells = 2*nr[0]*nr[1];
    M.get_arrays(&vertices, &cells, &dof, &is_bound, &bdry_vals);
    
    cudaEventRecord(start, 0);
    if(!dense) {
        std::cout << "      Sparsity pass...\n";
        M.sparsity_pass(valsLCPU, rowPtrLCPU, colIndLCPU, nnz);
    }
    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&t.sparsity_scan, start, finish);
    
    //////////////////////////////////////////////////////////////////////////////////


    //////////// Allocating memory for mesh/stiffness matrix/stress vector //////////
    

    cudaEventRecord(start, 0);
    
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
    
    if(dense){
        stat = cudaMalloc( (void**)&L, order*order*sizeof(float));
        assert(stat == cudaSuccess);
    } else {
        stat = cudaMalloc( (void**)&valsL, nnz*sizeof(float)); 
        assert(stat == cudaSuccess);
        stat = cudaMalloc( (void**)&colIndL, nnz*sizeof(int)); 
        assert(stat == cudaSuccess);
        stat = cudaMalloc( (void**)&rowPtrL, (order+1)*sizeof(int)); 
        assert(stat == cudaSuccess);
    }
    stat = cudaMalloc( (void**)&b, order*sizeof(float));
    assert(stat == cudaSuccess);
    
    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&tau, start, finish);
    t.alloc += tau;

    if(timing){
        stat = cudaMalloc( (void**)&tau_d, num_cells*3*2);
        assert(stat == cudaSuccess);
    }

    ////////////////////////////////////////////////////////////////////////////////


    /////////////// Copying data of Mesh from Host to Device ///////////////////////
    
    std::cout << "      Copying data from host...\n";

    cudaEventRecord(start, 0);
    
    cudaMemcpy(vertices_gpu, vertices, 2*order*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(cells_gpu, cells, 3*num_cells*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dof_gpu, dof, 3*num_cells*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(is_bound_gpu, is_bound, order*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(bdry_vals_gpu, bdry_vals, order*sizeof(float), cudaMemcpyHostToDevice);
    
    ////////////////////////////////////////////////////////////////////////////////

   
    /////////////// Copying sparsity pattern if !dense ////////////////////////////
    
    if(!dense){
        cudaMemcpy(rowPtrL, &rowPtrLCPU[0], (order+1)*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(colIndL, &colIndLCPU[0], nnz*sizeof(int), cudaMemcpyHostToDevice);
    }

    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&t.transfer, start, finish);
    
    ///////////////////////////////////////////////////////////////////////////////


    /////////  DIMENSIONS OF SYSTEM => 1 block per cell, 1 thread per node ////////
    
    block_size_X = 1, block_size_Y = 3;
    dim3 dimBlock(block_size_X, block_size_Y);
    dim3 dimGrid((num_cells/dimBlock.x)+(!(num_cells%dimBlock.x)?0:1),
                (1/dimBlock.y)+(!(1%dimBlock.y)?0:1));
    
    shared = 31;
    
    ///////////////////////////////////////////////////////////////////////////////

    
    //////// Kernel to assemble Stiffness matrix and store in global mem //////////

    std::cout << "      Main stiffness matrix assembly kernel...\n";
    // cudaEventRecord(start, 0);
    
    if(dense) {
        assemble_gpu<<<dimGrid, dimBlock, shared*sizeof(float)>>>(L, b, vertices_gpu, 
                       cells_gpu, is_bound_gpu, bdry_vals_gpu, order, tau_d);
    } else {
        assemble_gpu_csr<<<dimGrid, dimBlock, shared*sizeof(float)>>>(valsL, rowPtrL, colIndL, 
                        b, vertices_gpu, cells_gpu, is_bound_gpu, bdry_vals_gpu, order, tau_d);
    }

    if(timing){
        array_max(tau_d, num_cells*3, tau.elem_mats);
        array_max(tau_d, num_cells*3, tau.assembly);
    }
    /*
    NEED TO PASS BACK 3 FLOATS FOR CONVERSION, ASSEMBLY & ELEM MATS
    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&, start, finish);
    */
    
    //////////////////////////////////////////////////////////////////////////////
     
    
    /////////// Solving linear system in dense, dense-sparse conversion, CSR format ////////

    std::cout << "      Solving linear system...\n";
    if(!dense) std::cout << "      (nnz = " << nnz << ")\n";
    
    cudaEventRecord(start, 0);
    
    if(dense){
        if(dnsspr)  dnsspr_solve(L, b, order);
        else        dense_solve(L, b, order);
    } else { 
        sparse_solve(valsL, rowPtrL, colIndL, b, order, nnz);
    }
    
    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&t.solve, start, finish);
    
    ///////////////////////////////////////////////////////////////////////////////////


    ////////////////// Transfer solution back to Host & tidy //////////////////////////

    std::cout << "      Transferring result back to host...\n";
    
    cudaEventRecord(start, 0);
    
    cudaMemcpy(u, b, order*sizeof(float), cudaMemcpyDeviceToHost);
    
    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&tau, start, finish);
    t.transfer += tau;
    
    cudaFree(vertices_gpu);      cudaFree(cells_gpu);   cudaFree(dof_gpu);
    cudaFree(is_bound_gpu);      cudaFree(bdry_vals_gpu);
    if(timing)  cudaFree(tau_d);
    if(dense)   cudaFree(L),     cudaFree(b);
    else        cudaFree(valsL), cudaFree(colIndL),     cudaFree(rowPtrL); 
    
    /////////////////////////////////////////////////////////////////////////////////
}
