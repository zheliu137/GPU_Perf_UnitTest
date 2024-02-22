LIBS+=-L/usr/local/cuda/lib64 -lcudart -lcurand

FLAGS=-g -O2 -gencode=arch=compute_70,code=sm_70

#objects:=cusolver_kernel_wrapper.o cusolver_kernel.o main.o
objects:=curand_kernel.o main.o

all: main.x
main.x: $(objects)
	nvcc -g -O3 -o $@ $^ $(LDFLAGS) $(LIBS)

curand_kernel.o: curand_kernel.cu
	nvcc $(FLAGS) -c -o $@ $^
main.o: main.cu
	nvcc $(FLAGS) -c -o $@ $^

%.o : %.c
	icc -c -o $*.o $<
clean:
	rm -f ./cu*_kernel_f *.mod *.o
