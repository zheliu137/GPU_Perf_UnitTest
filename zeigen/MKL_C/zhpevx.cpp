#include <complex>
#include <cstdlib>
#include <stdio.h>
#include "mkl.h"
#include <complex.h>
#include <ctime>

using namespace std;

//typedef struct{ double re; double im; } complex16;
// typedef complex<double> ComplexD;
typedef MKL_Complex16 ComplexD;
// typedef double complex ComplexD;
#define N 16
#define NLOOP 200000
int main()
{

int n, inca = N, incb = N, i, j;
int nmat=NLOOP;
// ComplexD a[N*(N+1)/2], b[N*(N+1)/2], c[N*N];
// ComplexD a[N*(N+1)/2];
ComplexD a[N*(N+1)/2];
const ComplexD cone = {1.0,0.00}, czero = {0.0,0.00};

n = N;

srand (static_cast <unsigned> (time(0)));


for( i = 0; i < n*(n+1)/2; i++ ){
  // for( j = 0; j < n; j++ ){
  double ar = (double) rand() / (double) RAND_MAX;
  double ai = (double) rand() / (double) RAND_MAX;
  a[i] = { ar, ai};
    // a[i] = { static_cast <double> (rand()) / static_cast <double> (RAND_MAX),
    //              static_cast <double> (rand()) / static_cast <double> (RAND_MAX)};
}

/*

lapack_int LAPACKE_zhpevx( int matrix_layout, char jobz, char range, char uplo,
lapack_int n, lapack_complex_double* ap, double vl, double vu, lapack_int il,
lapack_int iu, double abstol, lapack_int* m, double* w, lapack_complex_double* z,
lapack_int ldz, lapack_int* ifail );

*/

double clock_start = clock();
int neig;
double eig[n];
ComplexD eigvec[n*n];
int ifail = 0;
int info;
for (int i=0; i<nmat; i++){
  info = LAPACKE_zhpevx(LAPACK_COL_MAJOR, 'V', 'A', 'U', n, a, 0.0, 0.0, 0, 0, -1, &neig, eig, eigvec, n, &ifail );
  // printf("info : %d \n", info);
  // cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,  n, n, n, &cone, a, inca, b, incb, &czero, c, n );
}

double clock_stop = clock();

double time_elapsed = clock();

printf("lapacke_zhpevx elapsed time : %15.6fs CPU\t(%9d calls) %gs for each.\n",time_elapsed/CLOCKS_PER_SEC,nmat,time_elapsed/CLOCKS_PER_SEC/nmat);

// for( i = 0; i < n; i++ ){
//   for( j = 0; j < n; j++ ){
//     printf( "The complex matrix multiplication is: ( %6.2f, %6.2f)\n", c[j+i*n].real(), c[j+i*n].imag());
//   }
// }

// return 0;
}
