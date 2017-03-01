#!/bin/bash
for i in *.fasta
do
NUM_A=`grep -A 1 ">$1" "$i" |grep -v ">" | grep -io "A"|grep -ic "A"`
NUM_T=`grep -A 1 ">$1" "$i" |grep -v ">" | grep -io "T"|grep -ic "T"`
NUM_G=`grep -A 1 ">$1" "$i" |grep -v ">" | grep -io "G"|grep -ic "G"`
NUM_C=`grep -A 1 ">$1" "$i" |grep -v ">" | grep -io "C"|grep -ic "C"`
LENGTH=`grep -A 1 ">$1" "$i" |grep -v ">" | awk '{print length}'`
GC_CONTENT=`echo "(($[NUM_C] +  $[NUM_G])/$LENGTH) * 100" |bc -l|sed "s|[0-9][0-9]$||g"`

echo "\nThis is the base composition of the $i gene for $1:"
echo "Total length: $LENGTH"
echo "Number of As: $NUM_A"
echo "Number of Ts: $NUM_T"
echo "Number of Gs: $NUM_G"
echo "Number of Cs: $NUM_C"
echo "Your GC% is: $GC_CONTENT"
done
