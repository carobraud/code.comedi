#!/bin/bash

value_volt1=$1
value_volt2=$2

pi=$(echo "scale=10; 4*a(1)" | bc -l)

let value1=($value_volt1+10)*65535/20
let value2=($value_volt2+10)*65535/20

#((i=0 ; 10 - $i ; i++))

#for((i=0 ; 10 ; i++))

for ((;;))
do
./outp $value1 -c 0 -s 1
sleep 0.5
./outp $value2 -c 0 -s 1
sleep 0.5
done
