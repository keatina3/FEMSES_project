#!/bin/bash

for i in 3 9 24 49; do
    for j in 3 9 24 49; do
        echo  "Testing NxM: $i x $j"
        for k in 1 8 16 27 32 164; do
            cuda-memcheck -n $i -m $j -M -d -b $k -c 1>> memcheck.out 2>>mem_error.log
            cuda-memcheck -n $i -m $j -d -b $k 1 -c >> memcheck.out 2>>mem_error.log
            cuda-memcheck -n $i -m $j -M -s -b $k -f -c 1>> memcheck.out 2>>mem_error.log
            cuda-memcheck -n $i -m $j -s -b $k 1 -f -c >> memcheck.out 2>>mem_error.log
            cuda-memcheck -n $i -m $j -M -b $k -f -c 1>> memcheck.out 2>>mem_error.log
            cuda-memcheck -n $i -m $j -b $k 1 -f -c >> memcheck.out 2>>mem_error.log
        done
    done
done
