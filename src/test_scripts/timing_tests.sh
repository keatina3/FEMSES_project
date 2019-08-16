#!/bin/bash

rm -f timings/* results/*

for i in 5 10 25 50 100 200 500 750; do
    echo "Testing NxM: $i x $i..."
    if [ $i -lt 500 ]
    then
        echo "Testing dense solvers..."
        ./fem_solver -n $i -m $i -d -v -t 1>> a.out 2>> error.log
    fi
    echo "Testing sparse solvers..."

    ./fem_solver -n $i -m $i -v -t 1>> a.out 2>> error.log
    
    echo "Testing conversion solver..."
    ./fem_solver -n $i -m $i -f -s -c -v -t 1>> a.out 2>> error.log
done
