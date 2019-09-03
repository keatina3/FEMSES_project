#!/bin/bash

#rm -f timings/* results/*

for i in 4 9 24 49 99 199 499 699; do

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
    
     
    if [ $i -lt 200 ]; then 
        echo "Testing FEMSES solver..."
        echo "      block_size_X = 1"
        
        for x in {1..4}; do
            ./fem_solver -n $i -m $i -M -c -g -c -v -t -b 1 1>> a.out 2>> error.log
        done
        for x in {1..4}; do
            ./fem_solver -n $i -m $i -g -c -v -t -b 1 1>> a.out 2>> error.log
        done
        
        for j in {2..300..6}; do
            echo "      block_size_X = $j"
            for x in {1..4}; do
                ./fem_solver -n $i -m $i -M -g -c -v -t -b $j 1>> a.out 2>> error.log
            done
            
            if [ $j -lt 220 ]; then
                #echo "  ...testing with Memory reconfiguration on"
                for x in {1..4}; do
                    ./fem_solver -n $i -m $i -g -c -v -t -b $j 1>> a.out 2>> error.log
                done
            fi
        done
    fi
    
    
    if [ $i -lt 500 ]; then
        echo "Testing RTX2080..."
        echo "Testing sparse solver..."    
        for x in {1..4}; do
            ./fem_solver -n $i -m $i -c -v -t -k -f 1>> a.out 2>> error.log
        done
        
        if [ $i -lt 200 ]; then
            echo "Testing dense solver..."    
            for x in {1..4}; do
                ./fem_solver -n $i -m $i -c -v -t -k -d -f 1>> a.out 2>> error.log
            done
            echo "Testing FEMSES..."
            for x in {1..4}; do
                ./fem_solver -n $i -m $i -c -g -v -t -k 1>> a.out 2>> error.log
            done
        fi
    fi

done
