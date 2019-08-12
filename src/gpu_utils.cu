#include <cassert>
#include <iostream>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cusolverSp.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include "gpu_utils.h"

void dnsspr_solve(float *L, float *b, int order){
    cusparseHandle_t handle = NULL;
    cusolverSpHandle_t handleS = NULL;
    cudaStream_t stream = NULL;
    cusparseMatDescr_t desc = NULL;
    cusparseStatus_t status = CUSPARSE_STATUS_SUCCESS;
    cusolverStatus_t status2 = CUSOLVER_STATUS_SUCCESS;
    cudaError_t cudaStat1 = cudaSuccess;
    const cusparseDirection_t dir = CUSPARSE_DIRECTION_ROW;

    int* csrRowPtrL = NULL;
    int* csrColIndL = NULL;
    float* csrValL  = NULL;
    int* nnzLrow;
    int nnzL;
    const float err = 1E-6;
    int reorder = 0;
    int singularity;

    cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    assert(cudaSuccess == cudaStat1);
    
    status = cusparseCreate(&handle);
    assert(CUSPARSE_STATUS_SUCCESS == status);
    std::cout << "lefence issue here\n";
    
    status2 = cusolverSpCreate(&handleS);
    assert(CUSPARSE_STATUS_SUCCESS == status2);
    
    status = cusparseSetStream(handle, stream);
    assert(CUSPARSE_STATUS_SUCCESS == status);
    
    status = cusparseCreateMatDescr(&desc);
    assert(CUSPARSE_STATUS_SUCCESS == status);
    
    cusparseSetMatIndexBase(desc, CUSPARSE_INDEX_BASE_ZERO);
    cusparseSetMatType(desc, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatFillMode(desc, CUSPARSE_FILL_MODE_LOWER);
    
    cudaStat1 = cudaMalloc( (void**)&csrRowPtrL, sizeof(int)*(order+1));
    assert(cudaSuccess == cudaStat1);
    cudaStat1 = cudaMalloc( (void**)&nnzLrow, sizeof(int)*(order));
    assert(cudaSuccess == cudaStat1);
    
    status = cusparseSnnz(handle, dir, order, order, desc, L, order, nnzLrow, &nnzL);
    assert(CUSPARSE_STATUS_SUCCESS == status);
    
    cudaStat1 = cudaMalloc( (void**)&csrValL, nnzL*sizeof(float));
    assert(cudaSuccess == cudaStat1);
    cudaStat1 = cudaMalloc( (void**)&csrColIndL, nnzL*sizeof(float));
    assert(cudaSuccess == cudaStat1);

    cusparseSdense2csr(handle,order,order,desc,L,order,nnzLrow,csrValL,csrRowPtrL,csrColIndL);

    status2 = cusolverSpSetStream(handleS, stream);
    assert(CUSOLVER_STATUS_SUCCESS == status);

    status2 = cusolverSpScsrlsvchol(handleS, order, nnzL, desc, csrValL, csrRowPtrL,
                                            csrColIndL, b, err, reorder, b, &singularity);
    
    cusparseDestroy(handle);
    cusolverSpDestroy(handleS);
    cudaStreamDestroy(stream);
    cusparseDestroyMatDescr(desc);
}

void sparse_solve(float *valsL,int *rowPtrL, int *colPtrL, float *b, int order, int nnz){
    cusolverSpHandle_t handle = NULL;
    cudaStream_t stream = NULL;
    cusparseMatDescr_t desc = NULL;
    cusparseStatus_t status = CUSPARSE_STATUS_SUCCESS;
    cusolverStatus_t status2 = CUSOLVER_STATUS_SUCCESS;
    cudaError_t cudaStat1 = cudaSuccess;

    const float err = 1E-6;
    int reorder = 0;
    int singularity;

    cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    assert(cudaSuccess == cudaStat1);
    
    status2 = cusolverSpCreate(&handle);
    std::cout << status2 << std::endl;
    assert(CUSPARSE_STATUS_SUCCESS == status2);

    status = cusparseCreateMatDescr(&desc);
    assert(CUSPARSE_STATUS_SUCCESS == status);
    
    cusparseSetMatIndexBase(desc, CUSPARSE_INDEX_BASE_ZERO);
    cusparseSetMatType(desc, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatFillMode(desc, CUSPARSE_FILL_MODE_LOWER);
    
    status2 = cusolverSpSetStream(handle, stream);
    assert(CUSOLVER_STATUS_SUCCESS == status);

    status2 = cusolverSpScsrlsvchol(handle, order, nnz, desc, valsL, rowPtrL,
                                            colPtrL, b, err, reorder, b, &singularity);
    
    cusolverSpDestroy(handle);
    cudaStreamDestroy(stream);
    cusparseDestroyMatDescr(desc);
}

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
    
    std::cout << "Dense testing not sparse\n";
    status = cusolverDnCreate(&handle);
    assert(CUSOLVER_STATUS_SUCCESS == status);
    std::cout << "lefence issue here\n";
    
    cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    //cudaStat1 = cudaStreamCreate(&stream);
    assert(cudaSuccess == cudaStat1);
    
    status = cusolverDnSetStream(handle, stream);
    assert(CUSOLVER_STATUS_SUCCESS == status); 

    status = cusolverDnSpotrf_bufferSize(handle, uplo, order, L, order, &bufferSize);
    assert(CUSOLVER_STATUS_SUCCESS == status);
    
    cudaMalloc( (void**)&info, sizeof(int));
    cudaMalloc( (void**)&buffer, bufferSize*sizeof(float));
    cudaMemset(info, 0, sizeof(int));

    status = cusolverDnSpotrf(handle, uplo, order, L, order, buffer, bufferSize, info);
    cudaStat1 = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == status);
    assert(cudaSuccess == cudaStat1);
    
    cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost);

    status = cusolverDnSpotrs(handle, uplo, order, nrhs, L, order, b, order, info);
    cudaStat1 = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == status);
    assert(cudaSuccess == cudaStat1);
      
    cusolverDnDestroy(handle);
    cudaStreamDestroy(stream);
}

void error_dot_prod(float *a, float *b, int n, float &x){
    cublasHandle_t handle;
    cublasStatus_t status = CUBLAS_STATUS_SUCCESS;
    const float alpha = -1.0;
    float *res;

    cudaMalloc( (void**)&res, n*sizeof(float));
    cudaMemcpy(res, b, n*sizeof(float), cudaMemcpyDeviceToDevice);

    status = cublasCreate(&handle);
    assert(status == CUBLAS_STATUS_SUCCESS);
 
    status = cublasSaxpy(handle, n, &alpha, a, 1, res, 1); 
    assert(status == CUBLAS_STATUS_SUCCESS);

    status = cublasSnrm2(handle, n, res, 1, &x);
    assert(status == CUBLAS_STATUS_SUCCESS);
    /*
    std::cout << "x " << x <<std::endl;
    status = cublasSdot(handle, n, a, 1, b, 1, &x); 
    assert(status == CUBLAS_STATUS_SUCCESS);
    std::cout << "x " << x <<std::endl;
    */
    status = cublasDestroy(handle);
    assert(status == CUBLAS_STATUS_SUCCESS); 
}
