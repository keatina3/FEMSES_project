#!/bin/bash

rm -f timings/* results/*

for i in 4 9 24 49 99 199 499 699; do
    echo "Testing NxM: $i x $i..."
    
    if [ $i -lt 500 ]; then
        echo "Testing dense solvers..."
        echo "      block_size_X = 1"
        ./fem_solver -n $i -m $i -d -v -t -f -b 1 1>> a.out 2>> error.log
        
        for i in {2..400}; do
            echo "      block_size_X = $j"
            ./fem_solver -n $i -m $i -d -v -t -c -f -b $j 1>>a.out 2>>error.log
        done

    fi
    
    echo "Testing sparse solvers..."
    echo "      block_size_X = 1"
    ./fem_solver -n $i -m $i -v -t -f -b 1 1>> a.out 2>> error.log
    for j in {2..400}; do
        echo "      block_size_X = $j"
        ./fem_solver -n $i -m $i -v -t -c -f -b $j 1>> a.out 2>> error.log
    done

    echo "Testing conversion solver..."
    echo "      block_size_X = 1"
    ./fem_solver -n $i -m $i -f -s -c -v -t -f -b 1 1>> a.out 2>> error.log
    for j in {2..400}; do
        echo "      block_size_X = $j"
        ./fem_solver -n $i -m $i -v -t -s -c -f -b $j 1>> a.out 2>> error.log
    done

done
