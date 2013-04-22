#!/bin/bash

#value_volt1=$1
#value_volt2=$2

pi=$(echo "scale=10; 4*a(1)" | bc -l)

f=10
resol=100
step=$(echo "scale=10; (1/$f)/$resol" | bc -l) 
#let value2=($value_volt2+10)*65535/20

#((i=0 ; 10 - $i ; i++))

#for((i=0 ; 10 ; i++))

input=0
for ((;;))
do
input=$(echo "scale=10; $input+$step" | bc -l) 

value_volt1=$(echo "scale=10; s(2*$pi*$f*$input)" | bc -l)
value1=$(echo "scale=10; ($value_volt1+10)*65535/20" | bc -l)

./outp $value1 -c 0 -s 1
#sleep 0.5
#./outp $value2 -c 0 -s 1
#sleep 0.5
done
