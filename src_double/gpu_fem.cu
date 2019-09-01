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

#include <cstdio>

//////////// Calculates area of triangle, given coordinates ////////////////
__device__ double area(double *xi){
    double tmp = 0.0;

    tmp += xi[0]*(xi[4]*xi[8] - xi[5]*xi[7]);
    tmp -= xi[1]*(xi[3]*xi[8] - xi[5]*xi[6]);
    tmp += xi[2]*(xi[3]*xi[7] - xi[4]*xi[6]);

    return 0.5*tmp;
}
///////


///////////// Assembles element matrix and stress vector /////////////////////
// Calculated and stored on shared memory
__device__ void assemble_elem(
                double *vertices, 
                int *cells,  
                int *is_bound, 
                double *bdry_vals, 
                double *temp1, 
                int idx, 
                int idy)
{
    double *Le, *be, *xi, *consts;
    int v;
    double bound;
    int offset = 28*threadIdx.x;

    Le = &temp1[offset];                 // element matrix
    be = &temp1[offset + 9];             // element stress vector
    xi = &temp1[offset + 12];            // matrix of global coordinates
    consts = &temp1[offset + 21];        // stores beta, gamma and area seen in serial version
    
    v = cells[(idx*3) + idy];            // global node number 
    
    
    //////////// Assigning global coordinates //////////////////
    
    xi[3*idy] = 1.0;
    xi[(3*idy) + 1] = vertices[2*v]; 
    xi[(3*idy) + 2] = vertices[(2*v)+1];

    __syncthreads();
    
    // using 1 thread to calculate area //
    if(idy==0)
        consts[6] = area(xi);
    
    ///////////////////////////////////////////////////////////


    ///////////////// calculating beta, gamma //////////////////
    
    consts[idy] = xi[ 3*((idy+1)%3) + 2] - xi[ 3*((idy+2)%3) + 2];
    consts[idy + 3] = xi[ 3*((idy+2)%3) +1] - xi[ 3*((idy+1)%3) + 1];
    
    __syncthreads();
    
    ////////////////////////////////////////////////////////////


    ///////////////////// Calculating LHS & RHS //////////////////

    be[idy] = 0.0;      // 0.0 intially until BCs enforced //
    for(int i=0; i<=idy; i++){
        Le[(3*idy)+i] = 0.25*consts[6]*consts[6]*consts[6] * (consts[idy]*consts[i] 
                                    + consts[idy + 3]*consts[i+3]);
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
                atomicAddDouble(&be[j], (-1)*Le[(3*j) + idy]*bound);
                atomicExchDouble(&Le[(3*j) + idy], 0.0);
                atomicExchDouble(&Le[(3*idy) + j], 0.0);
            }
        }
        __syncthreads();
        Le[(3*idy)+idy] = 1.0;
        atomicExchDouble(&be[idy], bound);
    }                            
    
    /////////////////////////////////////////////////////////////////
}
////////


///////////////////// Assembles stiffness matrix and stress vector //////////
// Takes shared memory element matrices and vectors
// and maps them back to global memory as an assembled stiffness matrix
// Stored in dense format
__device__ void assemble_mat(
                double *L, 
                double *b, 
                double *vertices, 
                int *dof, 
                double *temp1, 
                int idx, 
                int idy, 
                int order)
{
    double *Le, *be;
    // int* dof_r;
    int offset = 28*threadIdx.x;
    int dof_r[3];

    Le = &temp1[offset];
    be = &temp1[offset + 9];
    // dof_r = (int*)&temp1[offset+12];  // stores in shared memory, global node numbers for 3 nodes

    if(idx==0 && idy == 0){
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                printf("%lf ", Le[(i*3) + j]);
            }
            printf("\n");
        }
    }

    ///////////////// Assigning global node numbers //////////////////
    
    dof_r[0] = dof[(idx*3)];
    dof_r[1] = dof[(idx*3)+1];
    dof_r[2] = dof[(idx*3)+2]; 
    
    // dof_r[idy] = dof[(idx*3) + idy];
    // __syncthreads();

    ////////////////////////////////////////////////////////////////
  

    ///////// Mapping values back to global stiffness matrix /////// 
    
    atomicAddDouble(&b[dof_r[idy]], be[idy]);
    
    for(int i=0; i<=idy; i++){
        atomicAddDouble(&L[(order*dof_r[idy]) + dof_r[i]], 
                            Le[(3*idy) + i]);
        if(i==idy)      // avoid double adding diagonal element
            continue;
        else {
            atomicAddDouble(&L[(order*dof_r[i]) + dof_r[idy]], 
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
                double *valsL, 
                int *rowPtrL, 
                int *colIndL, 
                double *b, 
                double *vertices, 
                int *dof, 
                double *temp1, 
                int idx, 
                int idy, 
                int order)
{
    double *Le, *be;
    // int* dof_r;
    int row;
    int *tmp1, *tmp2;
    int off = 0;
    int off_mem = 28*threadIdx.x;
    int dof_r[3];

    Le = &temp1[off_mem];
    be = &temp1[off_mem + 9];
    // dof_r = (int *)&temp1[off_mem + 12];

    if(idx==0 && idy == 0){
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                printf("%lf ", Le[(i*3) + j]);
            }
            printf("\n");
        }
    }

    ///////////////// Assigning global node numbers //////////////////
    
    
    dof_r[0] = dof[(idx*3)];
    dof_r[1] = dof[(idx*3)+1];
    dof_r[2] = dof[(idx*3)+2]; 
    
    // dof_r[idy] = dof[(idx*3) + idy];
    // __syncthreads();
    
    //////////////////////////////////////////////////////////////////


    ///////// Mapping values back to global stiffness matrix ////////
    
    atomicAddDouble(&b[dof_r[idy]], be[idy]);
    
    row = rowPtrL[dof_r[idy]];
    tmp1 = &colIndL[row];
    for(int i=0; i<3; i++){
        tmp2 = tmp1;
        while(*tmp2 != dof_r[i]){
            off++;
            tmp2++;
        }
        atomicAddDouble(&valsL[row + off], Le[(3*idy) + i]);
        off = 0;
    }

    ///////////////////////////////////////////////////////////////////
}
/////////


////////////// Kernel to calculate elements and assemble global stiffness matrix /////////////
// Assembled in dense format
__global__ void assemble_gpu(
                double *L, 
                double *b, 
                double *vertices, 
                int *cells, 
                int *is_bound, 
                double *bdry_vals, 
                int order,
                int num_cells,
                long long *tau_d,
                int timing)
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;      // idx = cell number
    int idy = blockIdx.y*blockDim.y + threadIdx.y;      // idy = local node number
    extern __shared__ double temp1[];            // shared mem to store elem mats & constants

    if(idx < num_cells && idy < blockDim.y){
        long long start = clock64();
        assemble_elem(vertices, cells, is_bound, bdry_vals, temp1, idx, idy);
        __syncthreads();
        long long end = clock64();
        if(timing)  tau_d[(idx*blockDim.y) + idy] = (end - start);

        start = clock64();
        assemble_mat(L, b, vertices, cells, temp1, idx, idy, order);
        end = clock64();
        if(timing)  tau_d[(num_cells*blockDim.y) + (idx*blockDim.y) + idy] = (end - start);
    }
}
///////


////////////// Kernel to calculate elements and assemble global stiffness matrix /////////////
// Assembled in CSR format
// see comments above
__global__ void assemble_gpu_csr(
                double *valsL, 
                int *rowPtrL, 
                int *colIndL, 
                double *b, 
                double *vertices, 
                int *cells, 
                int *is_bound, 
                double *bdry_vals, 
                int order,
                int num_cells,
                long long *tau_d,
                int timing)
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;
    extern __shared__ double temp1[];

    if(idx < num_cells && idy < blockDim.y){
        long long start = clock64();
        assemble_elem(vertices, cells, is_bound, bdry_vals, temp1, idx, idy);
        __syncthreads();
        long long end = clock64();
        if(timing)  tau_d[(idx*blockDim.y) + idy] = (end - start);

        start = clock64();
        assemble_mat_csr(valsL, rowPtrL, colIndL, b, vertices, cells, temp1, idx, idy, order);
        end = clock64();
        if(timing)  tau_d[(num_cells*blockDim.y) + (idx*blockDim.y) + idy] = (end - start);
    }
}    
////////


///////////////// C++ Function to be invoked from host to apply FEM to PDE ///////////////////
// Applies standard approach of assmebling stiffness matrix and
// decomposing the linear system on the GPU to solve
extern void gpu_fem(double *u, Mesh &M, Tau &t, int &reconfig){
    int nr[2];
    int order, num_cells, nnz;
    int block_size_Y, shared;
    double *vertices_gpu, *vertices;
    int *cells_gpu, *cells;
    int *dof_gpu, *dof;
    int *is_bound_gpu, *is_bound;
    double *bdry_vals_gpu, *bdry_vals;
    double *L, *b, *valsL;
    int *rowPtrL, *colIndL;
    std::vector<double> valsLCPU;
    std::vector<int> rowPtrLCPU;
    std::vector<int> colIndLCPU;
    cudaError_t stat;
    cudaEvent_t start, finish;
    double tau = 0.0, sprs_tau = 0.0, dummy=0.0;
    long long *tau_d;
    int threads, shrd_mem;

    std::cout << GREEN "\nGPU Solver...\n" RESET;
    
    stat = cudaSetDevice(k);
    assert(stat == cudaSuccess);
    
    cudaEventCreate(&start);
    cudaEventCreate(&finish);

    //////////////////////// Gathering info from mesh ////////////////////////////////

    M.get_recs(nr);

    order = (nr[0]+1)*(nr[1]+1);
    num_cells = 2*nr[0]*nr[1];
    M.get_arrays(&vertices, &cells, &dof, &is_bound, &bdry_vals);
    
    if(!dense) {
        std::cout << "      Sparsity pass...\n";
        M.sparsity_pass(valsLCPU, rowPtrLCPU, colIndLCPU, nnz, dummy, sprs_tau);
    }
    t.sparsity_scan = sprs_tau;

    //////////////////////////////////////////////////////////////////////////////////


    //////////// Allocating memory for mesh/stiffness matrix/stress vector ///////////
    

    cudaEventRecord(start, 0);
    
    stat = cudaMalloc( (void**)&vertices_gpu, 2*order*sizeof(double));
    assert(stat == cudaSuccess);
    stat = cudaMalloc( (void**)&cells_gpu, 3*num_cells*sizeof(int));
    assert(stat == cudaSuccess);
    stat = cudaMalloc( (void**)&dof_gpu, 3*num_cells*sizeof(int));
    assert(stat == cudaSuccess);
    stat = cudaMalloc( (void**)&is_bound_gpu, order*sizeof(int));
    assert(stat == cudaSuccess);
    stat = cudaMalloc( (void**)&bdry_vals_gpu, order*sizeof(double));
    assert(stat == cudaSuccess);
    
    if(dense){
        stat = cudaMalloc( (void**)&L, order*order*sizeof(double));
        assert(stat == cudaSuccess);
    } else {
        stat = cudaMalloc( (void**)&valsL, nnz*sizeof(double)); 
        assert(stat == cudaSuccess);
        stat = cudaMalloc( (void**)&colIndL, nnz*sizeof(int)); 
        assert(stat == cudaSuccess);
        stat = cudaMalloc( (void**)&rowPtrL, (order+1)*sizeof(int)); 
        assert(stat == cudaSuccess);
    }
    stat = cudaMalloc( (void**)&b, order*sizeof(double));
    assert(stat == cudaSuccess);
    
    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime((float*)&t.alloc, start, finish);

    if(timing){
        stat = cudaMalloc( (void**)&tau_d, num_cells*3*2*sizeof(long long));
        assert(stat == cudaSuccess);
    }

    ////////////////////////////////////////////////////////////////////////////////


    /////////////// Copying data of Mesh from Host to Device ///////////////////////
    
    std::cout << "      Copying data from host...\n";

    cudaEventRecord(start, 0);
    
    stat = cudaMemcpy(vertices_gpu, vertices, 2*order*sizeof(double), cudaMemcpyHostToDevice);
    assert(stat == cudaSuccess);
    stat = cudaMemcpy(cells_gpu, cells, 3*num_cells*sizeof(int), cudaMemcpyHostToDevice);
    assert(stat == cudaSuccess);
    stat = cudaMemcpy(dof_gpu, dof, 3*num_cells*sizeof(int), cudaMemcpyHostToDevice);
    assert(stat == cudaSuccess);
    stat = cudaMemcpy(is_bound_gpu, is_bound, order*sizeof(int), cudaMemcpyHostToDevice);
    assert(stat == cudaSuccess);
    stat = cudaMemcpy(bdry_vals_gpu, bdry_vals, order*sizeof(double), cudaMemcpyHostToDevice);
    assert(stat == cudaSuccess);
    
    ////////////////////////////////////////////////////////////////////////////////

   
    /////////////// Copying sparsity pattern if !dense ////////////////////////////
    
    if(!dense){
        stat =cudaMemcpy(rowPtrL, &rowPtrLCPU[0], (order+1)*sizeof(int), cudaMemcpyHostToDevice);
        assert(stat == cudaSuccess);
        stat = cudaMemcpy(colIndL, &colIndLCPU[0], nnz*sizeof(int), cudaMemcpyHostToDevice);
        assert(stat == cudaSuccess);
    }

    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime((float*)&t.transfer, start, finish);
    
    ///////////////////////////////////////////////////////////////////////////////


    /////////  DIMENSIONS OF SYSTEM => b cells per block, 1 thread per node ////////
    
    block_size_Y = 3;
    dim3 dimBlock(block_size_X, block_size_Y);
    dim3 dimGrid((num_cells/dimBlock.x)+(!(num_cells%dimBlock.x)?0:1),
                (1/dimBlock.y)+(!(1%dimBlock.y)?0:1));
    
    shared = 28*block_size_X;
   
    cudaDeviceGetAttribute(&shrd_mem, cudaDevAttrMaxSharedMemoryPerBlock, k);
    cudaDeviceGetAttribute(&threads, cudaDevAttrMaxThreadsPerBlock, k);
    
    // testing if shared memory is over the max amount on card //
    if(shared * sizeof(double) > shrd_mem){
        error_log();
        std::cerr << "      Not enough shared memory on device to continue..." << std::endl;
        std::cerr << "              Shared memory requested: " 
                                            << shared * sizeof(double) << std::endl;
        std::cerr << "              Shared memory available: " << shrd_mem << std::endl;
        std::cerr << "      Exiting." << std::endl;
        std::exit(1);
    }
    
    // testing if requested block size if over max amount allowed on card //
    if(block_size_X * block_size_Y > threads){
        error_log();
        std::cerr << "      Too many threads requested per block..." << std::endl;
        std::cerr << "              Threads requested: " 
                                            << block_size_X * block_size_Y << std::endl;
        std::cerr << "              Max threads available: " << threads << std::endl;
        std::cerr << "      Exiting." << std::endl;
        std::exit(1);
    }
    
    // reconfiguring memory if shared has spare, to allow more per thread registers //
    reconfig = 0; 
    if(mem_config){
        if(shared * sizeof(double) < shrd_mem / 3){
            cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
            std::cout << "      Changed cache to prefer L1...\n";
            reconfig = 1;
        } else if(shared * sizeof(double) < shrd_mem / 2){
            cudaDeviceSetCacheConfig(cudaFuncCachePreferEqual);
            std::cout << "      Set cache to equal shared memory...\n";
            reconfig = 2;
        } else {
            cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
            reconfig = 0;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////

    
    //////// Kernel to assemble Stiffness matrix and store in global mem //////////

    std::cout << "      Main stiffness matrix assembly kernel...\n";
    
    if(dense) {
        assemble_gpu<<<dimGrid, dimBlock, shared*sizeof(double)>>>(L, b, vertices_gpu, 
                       cells_gpu, is_bound_gpu, bdry_vals_gpu, order, num_cells, tau_d, timing);
    } else {
        assemble_gpu_csr<<<dimGrid, dimBlock, shared*sizeof(double)>>>(valsL, rowPtrL, colIndL, 
                        b, vertices_gpu, cells_gpu, is_bound_gpu, 
                        bdry_vals_gpu, order, num_cells, tau_d, timing);
    }

    //////// Getting timings from device functions ///////
    
    int ind = 0;
    long long tmp;
    int clocks_per_sec;
    if(timing){
        cudaDeviceGetAttribute(&clocks_per_sec, cudaDevAttrClockRate, k);
        
        array_max((double*)tau_d, num_cells*3, ind);
        stat = cudaMemcpy(&tmp, &tau_d[ind], sizeof(double), cudaMemcpyDeviceToHost);
        assert(stat == cudaSuccess);
        t.elem_mats = (double) tmp / (clocks_per_sec / 1000);

        array_max((double*)&tau_d[num_cells*3], num_cells*3, ind);
        stat = cudaMemcpy(&tmp, &tau_d[(num_cells*3)+ind],sizeof(double), cudaMemcpyDeviceToHost);
        assert(stat == cudaSuccess);
        t.assembly = (double) tmp / (clocks_per_sec / 1000);
    }
    
    //////////////////////////////////////////////////////////////////////////////
     
    
    /////////// Solving linear system in dense, dense-sparse conversion, CSR format ////////

    std::cout << "      Solving linear system...\n";
    if(!dense) std::cout << "      (nnz = " << nnz << ")\n";
    
    cudaEventRecord(start, 0); 
    if(dense){
        if(dnsspr)  dnsspr_solve(L, b, order, start, finish, t.convert);
        else        dense_solve(L, b, order);
    } else { 
        sparse_solve(valsL, rowPtrL, colIndL, b, order, nnz);
    }
    
    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime((float*)&t.solve, start, finish);
    t.solve -= t.convert;
    
    ///////////////////////////////////////////////////////////////////////////////////


    ////////////////// Transfer solution back to Host & tidy //////////////////////////

    std::cout << "      Transferring result back to host...\n";
    
    cudaEventRecord(start, 0);
    
    cudaMemcpy(u, b, order*sizeof(double), cudaMemcpyDeviceToHost);
    
    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime((float*)&tau, start, finish);
    t.transfer += tau;
    
    cudaFree(vertices_gpu);      cudaFree(cells_gpu);   cudaFree(dof_gpu);
    cudaFree(is_bound_gpu);      cudaFree(bdry_vals_gpu);
    if(timing)  cudaFree(tau_d);
    if(dense)   cudaFree(L);
    else        cudaFree(valsL), cudaFree(colIndL),     cudaFree(rowPtrL); 
    cudaFree(b);

    /////////////////////////////////////////////////////////////////////////////////
}