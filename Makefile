all: ./run
CUDADIR :=  /usr/local/cuda
LDFLAGS := -L$(CUDADIR)/lib

CPPFLAGS := -I$(CUDADIR)/include
CFLAGS := -O2 --std=c++11
NVCC := nvcc
CC := $(NVCC)

run.o: run.cu

%.o: %.cu
	$(NVCC) -g -G -cudart shared $(CPPFLAGS) $(CFLAGS) -c $< -o $@   

.PHONY: clean
clean:
	rm ./run run.o

