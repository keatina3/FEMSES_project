CC = g++
NVCC = nvcc

CFLAGS = -g -Wall -W #-O3
NVCCFLAGS = -g -G --use_fast_math #-lefence #-DNDEBUG
#LDFLAGS = -lgsl -lgslcblas -lcusolver
# CAN I PUT THESE IN AS INCLUDES AND NOT NEED TO LOAD MODULE ??
LDFLAGS = -lm -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lcusolver -lcusparse

OBJECTS = *.o
SOURCES = *.cc
CU_SOURCES = *.cu

TARGET = fem_solver

LIBS = -L/usr/local/cuda-10.1/lib64 -L/home/support/apps/intel/18.0.4/mkl/lib/intel64/
INCPATH = -I/usr/include/ -I/usr/local/cuda-10.1/include/ -I/home/support/apps/intel/18.0.4/mkl/include/

all: objs cu_objs
	$(NVCC) $(OBJECTS) $(NVCCFLAGS) -o $(TARGET) $(LDFLAGS) $(INCPATH) $(LIBS)

cu_objs: $(CU_SOURCES)
	$(NVCC) $(NVCCFLAGS) -c $^

objs: $(SOURCES)
	$(CC) $(CFLAGS) -c $^ $(INCPATH)

#.PHONY: test mem_test clean
#
#test: all
#	./$(TARGET) -m 95000 -t 350 -s 17 1> results.out 2> error.log
#
#mem_test: all
# ADD cuda-memcheck test here
#	valgrind ./$(TARGET) -m 900 -t 350 -s 17

clean:
	$(RM) $(OBJECTS) $(TARGET) *.csv results.out error.log
