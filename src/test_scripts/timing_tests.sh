#!/bin/bash

#rm -f timings/* results/*

#for i in 4 9 24 49 99 199 499 699; do
for i in 24 49 99 199 499 699; do

    ############################################################################

    echo "Testing NxM: $i x $i..."
        
    if [ $i -lt 200 ]; then
        echo "Testing dense solvers..."
        echo "      block_size_X = 1"
        
        for x in {1..4}; do
            ./fem_solver -n $i -m $i -M -d -v -t -f -b 1 1>> a.out 2>> error.log
        done
        for x in {1..4} ; do
            ./fem_solver -n $i -m $i -d -v -t -c -f -b 1 1>> a.out 2>> error.log
        done

        for j in {2..300..6}; do
            echo "      block_size_X = $j"
            
            for x in {1..4}; do
            ./fem_solver -n $i -m $i -M -d -v -t -c -f -b $j 1>>a.out 2>>error.log
            done
            if [ $j -lt 220 ]; then
                #echo "  ...testing with Memory reconfiguration on"
                for x in {1..4}; do
                    ./fem_solver -n $i -m $i -d -v -t -c -f -b $j 1>>a.out 2>>error.log
                done
            fi
        done

    fi

    #############################################################################

    echo "Testing sparse solvers..."
    echo "      block_size_X = 1"
    
    for x in {1..4}; do
        ./fem_solver -n $i -m $i -M -v -t -f -b 1 1>> a.out 2>> error.log
    done
    for x in {1..4}; do
        ./fem_solver -n $i -m $i -v -t -f -c -b 1 1>> a.out 2>> error.log
    done

    for j in {2..300..6}; do
        echo "      block_size_X = $j"
        
        for x in {1..4}; do
            ./fem_solver -n $i -m $i -M -v -t -c -f -b $j 1>> a.out 2>> error.log
        done
        #echo "  ...testing with Memory reconfiguration on"
        if [ $j -lt 220 ]; then
            for x in {1..4}; do
                ./fem_solver -n $i -m $i -v -t -c -f -b $j 1>> a.out 2>> error.log
            done
        fi
    done
    
    ##############################################################################

    if [ $i -lt 200 ]; then 
        echo "Testing conversion solver..."
        echo "      block_size_X = 1"
        
        for x in {1..4}; do
            ./fem_solver -n $i -m $i -M -f -s -c -v -t -b 1 1>> a.out 2>> error.log
        done
        for x in {1..4}; do
            ./fem_solver -n $i -m $i -f -s -c -v -t -b 1 1>> a.out 2>> error.log
        done
        
        for j in {2..300..6}; do
            echo "      block_size_X = $j"
            
            for x in {1..4}; do
                ./fem_solver -n $i -m $i -M -v -t -s -c -f -b $j 1>> a.out 2>> error.log
            done
            
            if [ $j -lt 220 ]; then
                #echo "  ...testing with Memory reconfiguration on"
                for x in {1..4}; do
                    ./fem_solver -n $i -m $i -v -t -s -c -f -b $j 1>> a.out 2>> error.log
                done
            fi
        done
    fi

    ##############################################################################

    ' 
    echo "Testing FEMSES solver..."
    echo "      block_size_X = 1"
    if [ $i -lt 200 ]; then 
        ./fem_solver -n $i -m $i -M -g -c -v -t -b 1 1>> a.out 2>> error.log
        ./fem_solver -n $i -m $i -g -c -v -t -f -b 1 1>> a.out 2>> error.log
        for j in {2..300..6}; do
            echo "      block_size_X = $j"
            ./fem_solver -n $i -m $i -M -g -c -v -t -b $j 1>> a.out 2>> error.log
            if [ $j -lt 220 ]; then
                #echo "  ...testing with Memory reconfiguration on"
                ./fem_solver -n $i -m $i -g -c -v -t -b $j 1>> a.out 2>> error.log
            fi
        done
    fi
    '

done
