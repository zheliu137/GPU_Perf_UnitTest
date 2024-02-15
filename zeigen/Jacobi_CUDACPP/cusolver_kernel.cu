#include <cstdio>
#include <cstdlib>
#include <vector>
#include <complex>
#include <algorithm>

#include <cuda_runtime.h>
#include <cusolverDn.h>
#include "cuda_settings.h"

//extern "C"
//{
//#define DEBUG
//#define SINGLERUN

int cusolver_c_batch(const int m, cuDoubleComplex *A_, const int nmat, const int batchSize, const double tol_) {

    const int lda = m;
    //const int nmat = 512;
    int nbatch = nmat/batchSize;

    cusolverDnHandle_t cusolverH;
    cudaStream_t stream;
    syevjInfo_t syevj_params;

    //std::vector<cuDoubleComplex> V(lda * m * nmat); // eigenvectors
    //std::vector<double> W(m*nmat);       // eigenvalues
    //printf("Allocate pinned memory.\n");
    printf("solving %d %dx%d matrices by Jacobi method, with tol = %g\n",nmat,m,m,tol_);
    cuDoubleComplex *A; // matrix stored in pinned memory
    CUDA_CHECK(cudaMallocHost((void **)&A,sizeof(cuDoubleComplex)*lda * m * batchSize));
    //cuDoubleComplex A[lda * m * nmat]; // matrix stored in pinned memory
    //cuDoubleComplex V[lda * m * nmat]; // eigenvectors
    //cuDoubleComplex AMV[m*nmat]; // A*V
    cuDoubleComplex *V; // eigenvectors
    cuDoubleComplex *AMV; // A*V
    V = (cuDoubleComplex *)malloc (lda * m * batchSize * sizeof (*V));
    AMV = (cuDoubleComplex *)malloc (m * batchSize * sizeof (*AMV));
    double W[m*batchSize];       // eigenvalues
    
    //printf("Copy matrix to pinned memory.\n");
    std::copy(A_,A_+lda * m * batchSize,A);
    cuDoubleComplex *d_A;
    double *d_W;
    int *devInfo;
    cuDoubleComplex *d_work;
    int lwork;
    int info_gpu[batchSize];

    /* configuration of syevj  */
    //const double tol = 1.e-12;
    double tol = tol_;
    const int max_sweeps = 40;
    const cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvectors.
    const cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;

    cudaEvent_t start, stop;
    float elapsed_time;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));

    /* numerical results of syevj  */
    //double residual = 0;
    //int executed_sweeps = 0;
    int nloop = 1;
    for (int i=0;i<nloop;i++){ 
    // step 0: allocate device memory
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_A), sizeof(cuDoubleComplex) * lda * m*batchSize));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_W), sizeof(double) * m * batchSize));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&devInfo), sizeof(int) * batchSize));

    /* step 1: create cusolver handle, bind a stream */
    CUSOLVER_CHECK(cusolverDnCreate(&cusolverH));

    CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    CUSOLVER_CHECK(cusolverDnSetStream(cusolverH, stream));

    /* step 2: configuration of syevj */
    CUSOLVER_CHECK(cusolverDnCreateSyevjInfo(&syevj_params));

    /* default value of tolerance is machine zero */
    CUSOLVER_CHECK(cusolverDnXsyevjSetTolerance(syevj_params, tol));

    /* default value of max. sweeps is 100 */
    CUSOLVER_CHECK(cusolverDnXsyevjSetMaxSweeps(syevj_params, max_sweeps));

    /* step 3: copy A to device */
    CUDA_CHECK(
        cudaMemcpyAsync(d_A, A, sizeof(cuDoubleComplex) * lda * m *batchSize, cudaMemcpyHostToDevice, stream));

    /* step 4: query working space of syevj */
    CUSOLVER_CHECK(
          cusolverDnZheevjBatched_bufferSize(cusolverH, jobz, uplo, m, 
          d_A, lda, d_W, &lwork, syevj_params,batchSize));

    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_work), sizeof(cuDoubleComplex) * lwork)
);
    /* step 5: compute eigen-pair   */

    CUDA_CHECK(cudaEventRecord(start,stream));


    for (int j=0; j<nbatch; j++){ 
    CUSOLVER_CHECK(cusolverDnZheevjBatched(cusolverH, jobz, uplo, m, 
                    d_A, lda, d_W, d_work, lwork, devInfo, syevj_params, batchSize));
    }

    CUDA_CHECK(cudaEventRecord(stop,stream));
    CUDA_CHECK(cudaEventSynchronize(stop));
    CUDA_CHECK(cudaEventElapsedTime(&elapsed_time, start, stop));

    printf("batchsize: %d, nbatch: %d, CUDA event time: %gs, avg time per diag : %g \n",batchSize, nbatch,elapsed_time/1000.0,elapsed_time/1000.0/double(batchSize*nbatch));

    // step 6: check status, show eigenvalues, and eigenvectors 
      
    CUDA_CHECK(cudaMemcpyAsync(V, d_A, sizeof(cuDoubleComplex) *lda*m *batchSize, 
                                                cudaMemcpyDeviceToHost, stream));
    
    CUDA_CHECK(cudaMemcpyAsync(W, d_W, sizeof(double)* m *batchSize, 
                                                cudaMemcpyDeviceToHost, stream));

    CUDA_CHECK(cudaMemcpyAsync(&info_gpu, devInfo, sizeof(int) *batchSize, 
                                                cudaMemcpyDeviceToHost, stream));
    for (int i = 0; i < batchSize; i++) {
        if (0 == info_gpu[i]) {
#ifdef SINGLERUN
            printf("matrix %d: syevj converges \n", i);
#endif
        } else if (0 > info_gpu[i]) {
            printf("Error: %d-th parameter is wrong \n", -info_gpu[i]);
            exit(1);
        } else { 
            printf("WARNING: matrix %d, info = %d : sygvj does not converge \n",
 i, info_gpu[i]);        
      }
#ifdef DEBUG
      printf("Eigenvalue =  (ascending order)\n");
      for (int j = i*lda; j < (i+1)*m; j++) {
        printf("W[%d] = %E\n", j - i*lda, W[j]);
      }

      printf("V = \n");
      print_matrix(m, m, &V[ lda * m * i ], lda);
      printf("=====\n");
#endif
}

    // step 7 free device memory

    CUDA_CHECK(cudaFree(d_A));
    CUDA_CHECK(cudaFree(d_W));
    CUDA_CHECK(cudaFree(devInfo));

    CUDA_CHECK(cudaStreamDestroy(stream));
    CUSOLVER_CHECK(cusolverDnDestroySyevjInfo(syevj_params));
    CUSOLVER_CHECK(cusolverDnDestroy(cusolverH));

    }
    CUDA_CHECK(cudaFreeHost(A));
    CUDA_CHECK(cudaDeviceReset());
    // step 8: check results
    double residual;
    for (int i=0; i < batchSize; i++ ) {
       residual = 0.0;
       for (int j=0;j < m; j++) { 
#ifdef DEBUG
           printf("A * V(%d), W(%d) * V (%d)\n",j,j,j);
#endif
           for (int k=0; k < m; k++) { 
               AMV[k] = {0.0,0.0};
               for (int l=0; l < m; l++) { 
                   AMV[k] = cuCadd(AMV[k],
                            cuCmul(A_[k+l*m], V[i*m*lda+l+j*m]));
               }
#ifdef DEBUG
               printf("%0.2f + %0.2fj ", AMV[k].x, AMV[k].y);
               printf("%0.2f + %0.2fj ", 
                     W[i*m+j]*V[i*m*lda+k+j*m].x, W[i*m+j]*V[i*m*lda+k+j*m].y);
               printf("\n");
#endif
               residual = residual + abs(AMV[k].x-W[i*m+j]*V[i*m*lda+k+j*m].x)+
                                     abs(AMV[k].y-W[i*m+j]*V[i*m*lda+k+j*m].y);
           }
       }
    }
#ifdef SINGLERUN
    printf("residual = %e \n", residual);
#endif
    
    return EXIT_SUCCESS;
}
