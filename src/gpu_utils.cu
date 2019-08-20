#include <cassert>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cusolverSp.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <vector>
#include "mesh.h"
#include "utils.h"
#include "gpu_utils.h"


////////////// Converts dense matrix to CSR and solves linear system ///////////////////
// function takes input of L, b and dim (order)
// Returns soln overwritten on b
// uses cuSolverSp and cuSparse
void dnsspr_solve(float *L, float *b, int order){
    cusparseHandle_t handle = NULL;
    cusolverSpHandle_t handleS = NULL;
    cudaStream_t stream = NULL;
    cusparseMatDescr_t desc = NULL;
    cusparseStatus_t status = CUSPARSE_STATUS_SUCCESS;
    cusolverStatus_t status2 = CUSOLVER_STATUS_SUCCESS;
    cudaError_t cudaStat1 = cudaSuccess;                        
    const cusparseDirection_t dir = CUSPARSE_DIRECTION_ROW;     // row major

    int* csrRowPtrL = NULL;
    int* csrColIndL = NULL;
    float* csrValL  = NULL;
    int* nnzLrow;               // number of non zeros per row //
    int nnzL;
    const float err = EPS;
    int reorder = 0;
    int singularity;

    // Setting up streams //
    cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    assert(cudaSuccess == cudaStat1);
    
    // creating cuSparse handle //
    status = cusparseCreate(&handle);
    assert(CUSPARSE_STATUS_SUCCESS == status);
    
    // creating cuSolver handle //
    status2 = cusolverSpCreate(&handleS);
    assert(CUSPARSE_STATUS_SUCCESS == status2);
    
    // setting stream to cuSparse //
    status = cusparseSetStream(handle, stream);
    assert(CUSPARSE_STATUS_SUCCESS == status);
    
    // creating matrixx description //
    status = cusparseCreateMatDescr(&desc);
    assert(CUSPARSE_STATUS_SUCCESS == status);
    
    /*  setting matrix description:
        0-base ordering, lower fill for Cholesky, general format */
    cusparseSetMatIndexBase(desc, CUSPARSE_INDEX_BASE_ZERO);
    cusparseSetMatType(desc, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatFillMode(desc, CUSPARSE_FILL_MODE_LOWER);
    
    // allocating memory for CSR arrays //
    cudaStat1 = cudaMalloc( (void**)&csrRowPtrL, sizeof(int)*(order+1));
    assert(cudaSuccess == cudaStat1);
    cudaStat1 = cudaMalloc( (void**)&nnzLrow, sizeof(int)*(order));
    assert(cudaSuccess == cudaStat1);
    
    // count number of no   n-zeros //
    status = cusparseSnnz(handle, dir, order, order, desc, L, order, nnzLrow, &nnzL);
    assert(CUSPARSE_STATUS_SUCCESS == status);
    
    // malloc remaining CSR array //
    cudaStat1 = cudaMalloc( (void**)&csrValL, nnzL*sizeof(float));
    assert(cudaSuccess == cudaStat1);
    cudaStat1 = cudaMalloc( (void**)&csrColIndL, nnzL*sizeof(float));
    assert(cudaSuccess == cudaStat1);

    // convert from dense to sparse //
    cusparseSdense2csr(handle,order,order,desc,L,order,nnzLrow,csrValL,csrRowPtrL,csrColIndL);
    assert(CUSPARSE_STATUS_SUCCESS == status);
    
    // set stream to cuSolver //
    status2 = cusolverSpSetStream(handleS, stream);
    assert(CUSOLVER_STATUS_SUCCESS == status2);

    // solver using Cholesky factorisation //
    status2 = cusolverSpScsrlsvchol(handleS, order, nnzL, desc, csrValL, csrRowPtrL,
                                            csrColIndL, b, err, reorder, b, &singularity);
    assert(CUSOLVER_STATUS_SUCCESS == status2);
    
    // destroy handles, desc & stream //
    cusparseDestroy(handle);
    cusolverSpDestroy(handleS);
    cudaStreamDestroy(stream);
    cusparseDestroyMatDescr(desc);
}
///////


///////////////////// Solves linear system in CSR format //////////////////////////
// NOTE: see comments from fn above //
void sparse_solve(float *valsL,int *rowPtrL, int *colPtrL, float *b, int order, int nnz){
    cusolverSpHandle_t handle = NULL;
    cudaStream_t stream = NULL;
    cusparseMatDescr_t desc = NULL;
    cusparseStatus_t status = CUSPARSE_STATUS_SUCCESS;
    cusolverStatus_t status2 = CUSOLVER_STATUS_SUCCESS;
    cudaError_t cudaStat1 = cudaSuccess;

    const float err = EPS;
    int reorder = 0;
    int singularity;

    cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    assert(cudaSuccess == cudaStat1);
    
    status2 = cusolverSpCreate(&handle);
    assert(CUSOLVER_STATUS_SUCCESS == status2);

    status = cusparseCreateMatDescr(&desc);
    assert(CUSPARSE_STATUS_SUCCESS == status);
    
    cusparseSetMatIndexBase(desc, CUSPARSE_INDEX_BASE_ZERO);
    cusparseSetMatType(desc, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatFillMode(desc, CUSPARSE_FILL_MODE_LOWER);
    
    status2 = cusolverSpSetStream(handle, stream);
    assert(CUSOLVER_STATUS_SUCCESS == status2);
    
    status2 = cusolverSpScsrlsvchol(handle, order, nnz, desc, valsL, rowPtrL,
                                            colPtrL, b, err, reorder, b, &singularity);
    //assert(CUSOLVER_STATUS_EXECUTION_FAILED  == status2);
    assert(CUSOLVER_STATUS_SUCCESS == status2);
    
    cusolverSpDestroy(handle);
    cudaStreamDestroy(stream);
    cusparseDestroyMatDescr(desc);
}
////////


//////////////////// Solves linear system stored in dense format ////////////////////
// function takes input of L, b and dim (order)
// Returns soln overwritten on b
// uses cuSolverDn
void dense_solve(float *L, float *b, int order){
    cusolverDnHandle_t handle = NULL;
    cudaStream_t stream = NULL;
    cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
    cudaError_t cudaStat1 = cudaSuccess; 
    const cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
    const int nrhs = 1;
    float *buffer = NULL;
    int bufferSize = 0; 
    int *info = NULL;
    int h_info = 0;
    
    // setting cuSolver (Dense) handle //
    status = cusolverDnCreate(&handle);
    assert(CUSOLVER_STATUS_SUCCESS == status);
    
    cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    assert(cudaSuccess == cudaStat1);
    
    status = cusolverDnSetStream(handle, stream);
    assert(CUSOLVER_STATUS_SUCCESS == status); 

    // calculating buffer size needed for factorisation fn //
    status = cusolverDnSpotrf_bufferSize(handle, uplo, order, L, order, &bufferSize);
    assert(CUSOLVER_STATUS_SUCCESS == status);
    
    // allocating space for buffer on GPU //
    cudaStat1 = cudaMalloc( (void**)&info, sizeof(int));
    assert(cudaSuccess == cudaStat1);
    cudaStat1 = cudaMalloc( (void**)&buffer, bufferSize*sizeof(float));
    assert(cudaSuccess == cudaStat1);
    cudaMemset(info, 0, sizeof(int));

    // applying Cholesky factorisation to matrix //
    status = cusolverDnSpotrf(handle, uplo, order, L, order, buffer, bufferSize, info);
    cudaStat1 = cudaDeviceSynchronize();        // sync needed since non-blocking streams
    assert(CUSOLVER_STATUS_SUCCESS == status);
    assert(cudaSuccess == cudaStat1);
    
    cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost);

    // solving linear system - overwrites existing b //
    status = cusolverDnSpotrs(handle, uplo, order, nrhs, L, order, b, order, info);
    cudaStat1 = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == status);
    assert(cudaSuccess == cudaStat1);
      
    cusolverDnDestroy(handle);
    cudaStreamDestroy(stream);
}
///////


/////////////// Function gets error using 2-norm /////////////////
// Calculated using cuBLAS 
void error_dot_prod(float *a, float *b, int n, float &x){
    cublasHandle_t handle;
    cublasStatus_t status = CUBLAS_STATUS_SUCCESS;
    const float alpha = -1.0;

    // creating cuBLAS handle //
    status = cublasCreate(&handle);
    assert(status == CUBLAS_STATUS_SUCCESS);
 
    // y = ax + y
    // function sets b = b - a
    status = cublasSaxpy(handle, n, &alpha, a, 1, b, 1); 
    assert(status == CUBLAS_STATUS_SUCCESS);

    // gets <b,b> //
    status = cublasSnrm2(handle, n, b, 1, &x);
    assert(status == CUBLAS_STATUS_SUCCESS);
    
    // destroys handle //
    status = cublasDestroy(handle);
    assert(status == CUBLAS_STATUS_SUCCESS);
}
///////


/////////////// Function gets max value of array /////////////////
// Calculated using cuBLAS
void array_max(float *a, int n, float &max){
    cublasHandle_t handle;
    cublasStatus_t status = CUBLAS_STATUS_SUCCESS;
    
    // creating cuBLAS handle //
    status = cublasCreate(&handle);
    assert(status == CUBLAS_STATUS_SUCCESS);

    // getting maximum value of array //
    status = cublasIsamax(handle, n, a, 0, &max);
    assert(status == CUBLAS_STATUS_SUCCESS);
    
    // destroys handle //
    status = cublasDestroy(handle);
    assert(status == CUBLAS_STATUS_SUCCESS);
}
////////
