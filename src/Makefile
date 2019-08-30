CC = g++
NVCC = nvcc

CFLAGS = -std=c++11 -O3 #-W -Wall -g
NVCCFLAGS = -O3 --use_fast_math --relocatable-device-code=true #-DNDEBUG #-lefence
LDFLAGS = -lm -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lcusolver -lcusparse -lcublas

OBJECTS = *.o
SOURCES = *.cc
CU_SOURCES = *.cu

TARGET = fem_solver

LIBS = -L/usr/local/cuda-10.1/lib64 -L/home/support/apps/intel/18.0.4/mkl/lib/intel64/
INCPATH = -I/usr/include/ -I/usr/local/cuda-10.1/include/ -I/home/support/apps/intel/18.0.4/mkl/include/

all: objs cu_objs
	mkdir -p results timings
	$(NVCC) $(OBJECTS) $(NVCCFLAGS) -o $(TARGET) $(LDFLAGS) $(INCPATH) $(LIBS)

cu_objs: $(CU_SOURCES)
	$(NVCC) $(NVCCFLAGS) -c $^

objs: $(SOURCES)
	$(CC) $(CFLAGS) -c $^ $(INCPATH)

.PHONY: test timings mem_tests clean clean_all

test: all
	./fem_solver

timings: all
	./test_scripts/timing_tests.sh
	cat error.log

mem_tests: all
	./test_scripts/mem_tests.sh
	cat mem_error.log

clean:
	$(RM) $(OBJECTS) $(TARGET)

clean_all:
	$(RM) $(OBJECTS) $(TARGET) results/* timings/* *.out *.log
