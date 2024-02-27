LIBS+=-L/usr/local/cuda/lib64 -lcudart -lcurand -lcublas

# FLAGS=-g -O0 -gencode=arch=compute_70,code=sm_70
FLAGS=-g -O0 -gencode=arch=compute_86,code=sm_86

DFLAGS=-DDEBUG

#objects:=cusolver_kernel_wrapper.o cusolver_kernel.o main.o
objects_hk:=curand_kernel.o Hand_kernel.o

objects_cublas:=curand_kernel.o zgemv.o

all: Hand_kernel.x zgemv.x
Hand_kernel.x: $(objects_hk)
	nvcc -g -O0 -o $@ $^ $(LDFLAGS) $(DFLAGS) $(LIBS)
zgemv.x: $(objects_cublas)
	nvcc -g -O0 -o $@ $^ $(LDFLAGS) $(DFLAGS) $(LIBS)
curand_kernel.o: curand_kernel.cu
	nvcc $(FLAGS) $(DFLAGS) -c -o $@ $^
Hand_kernel.o: Hand_kernel.cu
	nvcc $(FLAGS) $(DFLAGS) -c -o $@ $^
zgemv.o: zgemv.cu
	nvcc $(FLAGS) $(DFLAGS) -c -o $@ $^

%.o : %.c
	icc -c -o $*.o $<
clean:
	rm -f ./cu*_kernel_f *.mod *.o
