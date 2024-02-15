icc   -c cptimer.c
ifort -c zhpevx.f90
ifort -qmkl cptimer.o zhpevx.o 
./a.out