icc  -traceback -c cptimer.c
ifort -traceback -c zhpevx.f90
ifort -mkl -traceback cptimer.o zhpevx.o -o zhpevx.x
