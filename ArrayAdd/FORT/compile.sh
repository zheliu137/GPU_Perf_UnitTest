icx  -traceback -c cptimer.c
ifx -traceback -c main.f90
ifx -qmkl -traceback cptimer.o main.o -o arrayadd.x
