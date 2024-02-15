#include <cstdlib>
#include <stdio.h>
#include "mkl.h"
#include <complex.h>
#include <ctime>

using namespace std;

//typedef struct{ double re; double im; } complex16;
typedef complex<double> ComplexD;
//typedef double complex ComplexD;
#define N 8
#define NLOOP 100
int main()
{

int n, inca = N, incb = N, i, j;
int nmat=NLOOP;
ComplexD a[N*N], b[N*N], c[N*N];
const ComplexD cone={1.0,0.0}, czero={0.0,0.0};

n = N;

srand (static_cast <unsigned> (time(0)));
double r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

for( i = 0; i < n; i++ ){
  for( j = 0; j < n; j++ ){

    a[j+i*n] = { static_cast <double> (rand()) / static_cast <float> (RAND_MAX),
                 static_cast <double> (rand()) / static_cast <float> (RAND_MAX)};
    b[j+i*n] = { static_cast <double> (rand()) / static_cast <float> (RAND_MAX),
                 static_cast <double> (rand()) / static_cast <float> (RAND_MAX)};

    // printf( "The complex matrix multiplication is: ( %6.2f, %6.2f)\n", a[j+i*n].real(), a[j+i*n].imag());
    // printf( "The complex matrix multiplication is: ( %6.2f, %6.2f)\n", b[j+i*n].real(), b[j+i*n].imag());
  }
}

/*

lapack_int LAPACKE_zhpevx( int matrix_layout, char jobz, char range, char uplo,
lapack_int n, lapack_complex_double* ap, double vl, double vu, lapack_int il,
lapack_int iu, double abstol, lapack_int* m, double* w, lapack_complex_double* z,
lapack_int ldz, lapack_int* ifail );

*/

double clock_start = clock();
for (int i=0;i<nmat;i++){
  cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,  n, n, n, &cone, a, inca, b, incb, &czero, c, n );
}

double clock_stop = clock();

double time_elapsed = clock();

printf("cblas_zgemm elapsed time : %15.6fs CPU\t(%9d calls) %gs for each.\n",time_elapsed/CLOCKS_PER_SEC,nmat,time_elapsed/CLOCKS_PER_SEC/nmat);

// for( i = 0; i < n; i++ ){
//   for( j = 0; j < n; j++ ){
//     printf( "The complex matrix multiplication is: ( %6.2f, %6.2f)\n", c[j+i*n].real(), c[j+i*n].imag());
//   }
// }

return 0;
}
