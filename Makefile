INCLUDES = -lm -I$(CUDA_INCLUDE)
CFLAGS   = $(INCLUDES)
LDFLAGS  = -L$(CUDA_LIBRARY) -lcudart

MEX  = mex
CXX  = g++
NVCC = nvcc

#GPUFLAGS = -arch=sm_35 #Deneb
GPUFLAGS = -arch=sm_70 #v100

MEX2  = MyTemporal_MEX
CUDA1 = kernel
WRAP1 = wrapper

all: $(MEX2).mexa64 

$(MEX2).mexa64: $(WRAP1).o $(CUDA1).o
	$(MEX) -output $(MEX2) $(MEX2).cpp $^ $(LDFLAGS) $(CFLAGS)

%.o: %.cu
	$(NVCC) -Xptxas="-v" -Xptxas -dlcm=ca -Xcompiler -fPIC -c -lineinfo $(GPUFLAGS) $(INCLUDES) $<

clean:
	rm -rf *.o $(MEX2).mexa64
