#include <math.h>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <complex>
#include <algorithm>
#include <sys/time.h>

#include <cuda_runtime.h>
#include <cusolverDn.h>
    
#ifdef DEBUG
#define CUSOLVER_CHECK(err) (HandlecusolverError(err, __FILE__, __LINE__))
#define CUDA_CHECK(err) (HandleError(err, __FILE__, __LINE__))
#define CUBLAS_CHECK(err) (HandleBlasError(err, __FILE__, __LINE__))
#else
#define CUSOLVER_CHECK(err) (err)
#define CUDA_CHECK(err) (err)
#define CUBLAS_CHECK(err) (err)
#endif

static void HandleBlasError(cublasStatus_t err, const char *file, int line)
{

    if (err != CUBLAS_STATUS_SUCCESS)
    {
        fprintf(stderr, "ERROR: %s in %s at line %d (error-code %d)\n",
                cublasGetStatusString(err), file, line, err);
        fflush(stdout);
        exit(-1);
    }
}



static void HandlecusolverError(cusolverStatus_t err, const char *file, int line )
{

    if (err != CUSOLVER_STATUS_SUCCESS)
    {
        fprintf(stderr, "ERROR: %d in %s at line %d, (error-code %d)\n",
                err, file, line, err);
        fflush(stdout);
        exit(-1);
    }
}

static void HandleError(cudaError_t err, const char *file, int line)
{

    if (err != cudaSuccess)
    {
        fprintf(stderr, "ERROR: %s in %s at line %d (error-code %d)\n",
                cudaGetErrorString(err), file, line, err);
        fflush(stdout);
        exit(-1);
    }
}

template <typename T> void print_matrix(const int &m, const int &n, const T *A, const int &lda);

template <> void print_matrix(const int &m, const int &n, const float *A, const int &lda) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            std::printf("%0.2f ", A[j * lda + i]);
        }
        std::printf("\n");
    }
}

template <> void print_matrix(const int &m, const int &n, const double *A, const int &lda) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            std::printf("%0.2f ", A[j * lda + i]);
        }
        std::printf("\n");
    }
}

template <> void print_matrix(const int &m, const int &n, const cuComplex *A, const int &lda) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            std::printf("%0.2f + %0.2fj ", A[j * lda + i].x, A[j * lda + i].y);
        }
        std::printf("\n");
    }
}

template <>
void print_matrix(const int &m, const int &n, const cuDoubleComplex *A, const int &lda) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            std::printf("%0.2f + %0.2fj ", A[j * lda + i].x, A[j * lda + i].y);
        }
        std::printf("\n");
    }
}

    cublasOperation_t char_to_cublas_trans(char trans)
    {
        cublasOperation_t cuTrans;
        switch (trans)
        {
        case 'n':
        case 'N':
            cuTrans = CUBLAS_OP_N;
            break;
        case 't':
        case 'T':
            cuTrans = CUBLAS_OP_T;
            break;
        case 'c':
        case 'C':
            cuTrans = CUBLAS_OP_C;
            break;
        default:
            exit(-1);
        }
        return cuTrans;
    }

    void matmul_strided_batched(int m, int n, int k, cuDoubleComplex *A, 
        cuDoubleComplex *B, cuDoubleComplex *C, int batch_count)
    {
        cuDoubleComplex alpha={1.0,0.0};
        cuDoubleComplex beta={0.0,0.0};

        long long int stridea=m*n;
        long long int strideb=n*k;
        long long int stridec=m*k;
        int lda=m;
        int ldb=n;
        int ldc=k;
        char transa = 'n', transb = 'n';
        cuDoubleComplex *dA = nullptr, *dB = nullptr, *dC = nullptr;
        cublasHandle_t blasHandle;
        cudaStream_t stream = NULL;

        CUDA_CHECK(cudaMallocAsync((void **)&dA, batch_count * m * k * sizeof(cuDoubleComplex), stream));
        CUDA_CHECK(cudaMallocAsync((void **)&dB, batch_count * k * n * sizeof(cuDoubleComplex), stream));
        CUDA_CHECK(cudaMallocAsync((void **)&dC, batch_count * m * n * sizeof(cuDoubleComplex), stream));

        CUBLAS_CHECK(cublasCreate(&blasHandle));
        CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
        CUBLAS_CHECK(cublasSetStream(blasHandle, stream));

        CUDA_CHECK(cudaMemcpyAsync(dA, A, batch_count * m * k * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice, stream));
        CUDA_CHECK(cudaMemcpyAsync(dB, B, batch_count * k * n * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice, stream));
        int batch_size=batch_count;
        int nstep = 1;
        int batch_step=batch_size/nstep;
        // CUDA timer
        cudaEvent_t start, stop;
        float elapsed_time;
        CUDA_CHECK(cudaEventCreate(&start));
        CUDA_CHECK(cudaEventCreate(&stop));

        for (int i=1;i<=nstep;i++){
            CUDA_CHECK(cudaEventRecord(start,stream));
            
            batch_count=i*batch_step;
            int nloop = 1000000;
            for (int j=1; j<=nloop;j++){
            CUBLAS_CHECK(cublasZgemmStridedBatched(blasHandle, char_to_cublas_trans(transa), char_to_cublas_trans(transb), m, n, k, &alpha, dA, lda, stridea, dB, ldb, strideb, &beta, dC, ldc, stridec, batch_count));
            }
            CUDA_CHECK(cudaEventRecord(stop,stream));
            CUDA_CHECK(cudaEventSynchronize(stop));
            CUDA_CHECK(cudaEventElapsedTime(&elapsed_time, start, stop));
            printf("batchsize: %d, CUDA event time: %gs, avg time per matmul : %g \n",batch_count*nloop,elapsed_time/1000.0,elapsed_time/1000.0/double(batch_count*nloop));
        }

        //CUDA_CHECK(cudaMemcpyAsync(C, dC, batch_count * m * n * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost, stream));

        CUDA_CHECK(cudaStreamSynchronize(stream));

        CUDA_CHECK(cudaFreeAsync(dA, stream));
        CUDA_CHECK(cudaFreeAsync(dB, stream));
        CUDA_CHECK(cudaFreeAsync(dC, stream));

        cublasDestroy(blasHandle);
        CUDA_CHECK(cudaStreamDestroy(stream));
    }
