LIBS+=-L/usr/local/cuda/lib64 -lcudart -lcublas -lcusolver -lcurand

#objects:=cusolver_kernel_wrapper.o cusolver_kernel.o main.o
objects:= main.o

FLAGS=-g -O0 -gencode=arch=compute_70,code=sm_70

DFLAGS=-DDEBUG

all: main.x
main.x: $(objects)
	nvcc -g -O3 -o $@ $^ $(LDFLAGS) $(DFLAGS) $(LIBS)
main.o: main.cu
	nvcc $(FLAGS) $(DFLAGS) -c -o $@ $^

%.o : %.c
	icc -c -o $*.o $<
clean:
	rm -f ./cu*_kernel_f *.mod *.o
