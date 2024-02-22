LIBS+=-L/usr/local/cuda/lib64 -lcudart -lcublas -lcurand

objects:=cublas_kernel.o curand_kernel.o cublas_zgemm.o

all: cublas_zgemm.x
cublas_zgemm.x: $(objects)
	nvcc -g -O3 -o $@ $^ $(LDFLAGS) $(LIBS)

cublas_kernel.o: cublas_kernel.cu
	nvcc -g -O2 -gencode=arch=compute_70,code=sm_70 -c -o $@ $^
curand_kernel.o: curand_kernel.cu
	nvcc -g -O2 -gencode=arch=compute_70,code=sm_70 -c -o $@ $^
cublas_zgemm.o: main.cu
	nvcc -g -O2 -gencode=arch=compute_70,code=sm_70 -c -o $@ $^

%.o : %.c
	icc -c -o $*.o $<
clean:
	rm -f *.mod *.o
