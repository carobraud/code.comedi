#!/bin/bash
# Date: 25/04/2013
# Author: Caroline Braud 
#
# Generate Triangular Signal for arduino AI (0-5volt)
# 

volt1=$1
volt2=$2
freq=$3
npoint=300

dt=$(echo "scale=14; 1/($npoint*$freq)" | bc -l)
dv=$(echo "scale=14; ($volt2-$volt1)/($npoint/2.)" | bc -l)


for ((;;))
do
for((i=0;i<npoint/2;i++))   
do
#echo "coucou" $i
volt1=$(echo "scale=14; $volt1+$dv" | bc -l)
int1=$(echo "scale=14;($volt1+10)*65535/20"| bc -l)
./outp $int1 -c 0 -s 1
done
volt1=$(echo "scale=14; $volt2" | bc -l)
for((i=0;i<npoint/2+1;i++))   
 do
volt1=$(echo "scale=14; $volt1-$dv" | bc -l)
int1=$(echo "scale=14;($volt1+10)*65535/20"| bc -l)
./outp $int1 -c 0 -s 1
done
done