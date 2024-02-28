//Example 1. Application Using C and cuBLAS: 1-based indexing
//-----------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <cuda_runtime.h>
#include "cublas_v2.h"
//#include "cuda_settings.h"

// void run_eig_wrapper_(const int N, cuDoubleComplex *x);
//void print_matrix(const int &m, const int &n, const cuDoubleComplex *A, const int &lda);
//void cusolver(int N); 
//void test(); 

int main (int argc, char* argv[]){
    cuDoubleComplex *A;
    //double *A;
    int nmat = 200000;
    int batchsize = 20000;
    if (argc > 1 ){
      N = strtol(argv[1],nullptr,0);
    }
    A = (cuDoubleComplex *)malloc(nmat*sizeof(cuDoubleComplex));


    // int device_count;

    // CUDA_CHECK(cudaGetDeviceCount(&device_count));

    // // CUDA_CHECK(cudaGetDeviceProperties());


    // current_device = 0;
    // while (current_device < device_count) {

    // }


    // findCudaDevice(argc, (const char **)argv)
    int device = 0;
    CUDA_CHECK(cudaSetDevice(device));
v
    cudaStream_t stream;

    // cuDoubleComplex *A;
    cuDoubleComplex *B, *C;
    cuDoubleComplex *d_A, *d_B, *d_C;
    int arraylen = 1000000;
    int nloop = 200;
    // int batchsize = 20000;
    if (argc > 1 ){
      arraylen = strtol(argv[1], nullptr, 0);
    }
    int N = arraylen;
    // double tol = exp10(-double(strtol(argv[1],nullptr,0)));
    // A = (cuDoubleComplex *)malloc(arraylen*sizeof(cuDoubleComplex));
    B = (cuDoubleComplex *)malloc(arraylen*sizeof(cuDoubleComplex));
    C = (cuDoubleComplex *)malloc(arraylen*sizeof(cuDoubleComplex));

    CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamDefault));

    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_A), N*sizeof(cuDoubleComplex)));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_B), N*sizeof(cuDoubleComplex)));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_C), N*sizeof(cuDoubleComplex)));


    cudaEvent_t start, stop;    
    float elapsed_time;
    double elapsed_time_sum;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));

    CUDA_CHECK(cudaMemset(d_A, 0, N*sizeof(cuDoubleComplex)));
    CUDA_CHECK(cudaMemcpy(d_B, B, N*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice ));
    CUDA_CHECK(cudaMemcpy(d_C, C, N*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice ));

    elapsed_time_sum=0.0;
    CUDA_CHECK(cudaEventRecord(start, stream));
    // for (int j=0; j<nloop; j++){ 
    dim3 block_dim = BLOCK_SIZE_1D;
    // dim3 block_dim = 32;
    dim3 grid_dim;
    grid_dim.x = (arraylen+block_dim.x-1)/block_dim.x;
    printf("%d %d\n", grid_dim.x,block_dim.x);
    ArrayAdd<<<grid_dim, block_dim>>>((ComplexD*)d_A, (ComplexD*)d_B, (ComplexD*)d_C, arraylen, nloop);
    // ArrayAdd<<<grid_dim, block_dim>>>(arraylen, nloop);
    CUDA_CHECK(cudaGetLastError());
    // }
    CUDA_CHECK(cudaEventRecord(stop, stream));
    CUDA_CHECK(cudaEventSynchronize(stop));
    CUDA_CHECK(cudaEventElapsedTime(&elapsed_time, start, stop));
    elapsed_time_sum+=elapsed_time;
    CUDA_CHECK(cudaDeviceSynchronize());

    double ms2s=0.001;
    double avgfac=1.0/double(arraylen)/double(nloop);
    // printf("Avgfac test : %g %g %g %g \n", double(long(arraylen)*long(nloop)), 1/double(arraylen*nloop), 1.0/double(arraylen*nloop), 1.0/double(arraylen)/double(nloop));
    printf("Array Length: %d, nloop: %d, CUDA event time: %gs, avg time for each sum : %gs \n", arraylen, nloop, elapsed_time_sum/1000.0, elapsed_time_sum*ms2s*avgfac);




    // printf("Success.\n");

    }
