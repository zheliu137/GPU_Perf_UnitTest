icc   -c cptimer.c
ifort -c zgemm.f90
ifort -mkl cptimer.o zgemm.o -o zgemm.x