#!/bin/bash
pwdd=`pwd`
input="runlist.txt"
while IFS= read -r line
do
./new_Bishoy.sh $line 
echo $line
done < "$input"

