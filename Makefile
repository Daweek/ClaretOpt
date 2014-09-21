GCC = gcc
CFLAGS = -g -O2 -ffast-math -funroll-loops

CUDA 					= /usr/local/cuda-6.0
CUDA_SDK			= /usr/local/cuda-6.0/samples
CUDA_4.1_SDK 	= /usr/local/cuda-4.1/NVIDIA_GPU_Computing_SDK/C
DSCUDA		= /usr/local/DSCUDA/dscudapkg1.7.5.1
CUDA_LIB = -L$(CUDA)/lib -lcudart

CUDAINC   = -I. -I$(CUDA)/include -I$(CUDA_SDK)/common/inc -I$(DSCUDA)/include/common/inc -I/usr/local/include -I/opt/X11/include -I/opt/X11/lib
CUDALIB   = -lcudart
GLLIB  		= -L/usr/local/lib -framework OpenGL -lGLUT  
LIB 			= $(CUDA_LIB) $(CUDALIB) $(CUDAGLLIB) $(GLLIB) -lm

NVCC       = $(CUDA)/bin/nvcc 
NVCC_FLAGS = --compile -use_fast_math -gencode arch=compute_30,code=sm_30 -O -g

TARGET = claret

all: $(TARGET)

claret : cras36.c mr3.o
	$(GCC) $(CFLAGS) $(CUDAINC) $< -o $@ mr3.o  $(LIB) 

claret_org : cras36_org.c mr3_org.o  
	$(GCC) -DVTGRAPE $(CFLAGS) $(CUDAINC) $< -o $@ mr3_org.o  $(LIB) 

mr3.o : mr3.cu
	$(NVCC) $(NVCC_FLAGS) $(CUDAINC)  $< $(CUDAGLLIB) -o $@ 

mr3_org.o : mr3_org.cu
	$(NVCC) $(NVCC_FLAGS) $(CUDAINC)  $< $(CUDAGLLIB) -o $@ 

clean:
	rm -rf  *.o $(TARGET) claret.d* claret_org.d*
