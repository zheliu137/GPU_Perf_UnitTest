LIBS+=-L/usr/local/cuda/lib64 -lcudart -lcublas -lcurand

LIBS_CPU+=-mkl

objects:= cublas_kernel.o \
  		  curand_kernel.o \
		  main.o

objects_cpu=zgemm.o

all: main_gpu.x main_cpu.x
main_gpu.x: $(objects)
	nvcc -g -O3 -o $@ $^ $(LDFLAGS) $(LIBS)
main_cpu.x: $(objects_cpu)
	icc -g -O3 -o $@ $^ $(LDFLAGS) $(LIBS_CPU)

cublas_kernel.o: cublas_kernel.cu
	nvcc -g -O2 -gencode=arch=compute_70,code=sm_70 -c -o $@ $^
curand_kernel.o: curand_kernel.cu
	nvcc -g -O2 -gencode=arch=compute_70,code=sm_70 -c -o $@ $^
main.o: main.cu
	nvcc -g -O2 -gencode=arch=compute_70,code=sm_70 -c -o $@ $^
zgemm.o: zgemm.cpp
	icc -g -O2 -c -o $@ $^
%.o : %.c
	icc -c -o $*.o $<
clean:
	rm -f ./cu*_kernel_f *.mod *.o
