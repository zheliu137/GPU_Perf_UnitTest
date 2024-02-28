//-----------------------------------------------------------
// CUDA CPP double complex Array Add performance test
//-----------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "cuda_settings.h"
#define IDX2F(i,j,ld) ((((j))*(ld))+((i)))
#define BLOCK_SIZE 32
#define BLOCK_SIZE_1D 1024

//void print_matrix(const int &m, const int &n, const cuDoubleComplex *A, const int &lda);
void createRandoms(int size, double *randomArray);

__global__ void ArrayAdd( ComplexD *A, ComplexD *B, ComplexD *C, int length, int nloop ){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if( idx<length ){
      // to be optimized
      for ( int iloop=0;iloop<nloop;iloop++){
        A[idx] += B[idx] + C[idx];
      }
    }
}

__global__ void ArrayAdd_v2( ComplexD *A, ComplexD *B, ComplexD *C, int length, int nloop ){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if( idx<length ){
      // to be optimized
      ComplexD AA=0.0;
      for ( int iloop=0;iloop<nloop;iloop++){
        AA += B[idx] + C[idx];
      }
      A[idx]=AA;
    }
}

int main (int argc, char* argv[]){

    // CUDA_CHECK(cudaGetDeviceProperties());

    int device = 0;
    CUDA_CHECK(cudaSetDevice(device));

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
    //A = (double *)malloc(pow(N,2)*sizeof(double));
    double *rand1;
    double *rand2;
    rand1 = (double *)malloc(arraylen*sizeof(double));
    rand2 = (double *)malloc(arraylen*sizeof(double));
    
    // printf("Generating %d by %d random matrix... \n",N,N);
    createRandoms(N, rand1);
    createRandoms(N, rand2);
    for (int i=0;i<N;i++){
      B[i] = {rand1[i],rand2[i]};
    } 
    createRandoms(N, rand1);
    createRandoms(N, rand2);
    for (int i=0;i<N;i++){
      C[i] = {rand1[i],rand2[i]};
      // C[i] = {1.5,1.5};
    } 

    CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamDefault));

    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_A), N*sizeof(cuDoubleComplex)));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_B), N*sizeof(cuDoubleComplex)));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_C), N*sizeof(cuDoubleComplex)));

    CUDA_CHECK(cudaMemset(d_A, 0, N*sizeof(cuDoubleComplex)));
    CUDA_CHECK(cudaMemcpy(d_B, B, N*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice ));
    CUDA_CHECK(cudaMemcpy(d_C, C, N*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice ));

    cudaEvent_t start, stop;    
    float elapsed_time;
    double elapsed_time_sum;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));

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
    // printf("Array Length: %d, nloop: %d, CUDA event time: %gs, avg time for each array 1 loop : %gs \n", arraylen, nloop, elapsed_time_sum/1000.0, elapsed_time_sum/1000.0/double(nloop));


    // version 2 
    elapsed_time_sum=0.0;
    CUDA_CHECK(cudaEventRecord(start, stream));
    // for (int j=0; j<nloop; j++){ 
    block_dim = BLOCK_SIZE_1D;
    // dim3 block_dim = 32;
    grid_dim = (arraylen+block_dim.x-1)/block_dim.x;
    printf("%d %d\n", grid_dim.x,block_dim.x);
    ArrayAdd_v2<<<grid_dim, block_dim>>>((ComplexD*)d_A, (ComplexD*)d_B, (ComplexD*)d_C, arraylen, nloop);
    // ArrayAdd<<<grid_dim, block_dim>>>(arraylen, nloop);
    CUDA_CHECK(cudaGetLastError());
    // }
    CUDA_CHECK(cudaEventRecord(stop, stream));
    CUDA_CHECK(cudaEventSynchronize(stop));
    CUDA_CHECK(cudaEventElapsedTime(&elapsed_time, start, stop));
    elapsed_time_sum+=elapsed_time;
    CUDA_CHECK(cudaDeviceSynchronize());

    printf("Ver2 Array Length: %d, nloop: %d, CUDA event time: %gs, avg time for each sum : %gs \n", arraylen, nloop, elapsed_time_sum/1000.0, elapsed_time_sum*ms2s*avgfac);

    // CUDA_CHECK(cudaMemcpy(A, d_A, N*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost ));

    // printf("The value of A[%d] is %30.15g + %30.15g i \n", 12000, A[11999].x, A[11999].y);
    /*
    */
}
