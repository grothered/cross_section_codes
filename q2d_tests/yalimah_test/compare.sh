#!/bin/bash
# Script to compare files here with files in reference_run 

for i in *;
    do
    x=$(diff $i reference_run/$i);
    lenx=$(echo ${#x});
        if [ $lenx -gt 0 ]; then
            echo $i ':' $x;
        else
            continue
            #echo $i ': fine'; 
        fi
    done
