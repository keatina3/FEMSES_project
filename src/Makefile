CC = g++
NVCC = nvcc

CFLAGS = -g -Wall -W #-O3
NVCCFLAGS = -g -G --use_fast_math
LDFLAGS = -lgsl -lgslcblas -lcusolver
# CAN I PUT THESE IN AS INCLUDES AND NOT NEED TO LOAD MODULE ??
#LDFLAGS = -lm -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

OBJECTS = *.o
SOURCES = *.cc
CU_SOURCES = *.cu

TARGET = fem_solver

LIBS = -L/usr/local/cuda-10.1/lib64/
INCPATH = -I/usr/include/ -I/usr/local/cuda-10.1/include/

all: objs cu_objs
	$(NVCC) $(OBJECTS) $(NVCCFLAGS) -o $(TARGET) $(LDFLAGS) $(LIBS) $(INCPATH)

cu_objs: $(CU_SOURCES)
	$(NVCC) $(NVCCFLAGS) -c $^

objs: $(SOURCES)
	$(CC) $(CFLAGS) -c $^

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
