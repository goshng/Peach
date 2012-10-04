#!/bin/bash
rm -rf status; mkdir status
rm -rf output; mkdir output
rm -rf log; mkdir log
echo -n "How many computing nodes do you wish to use? (e.g., 3) "
read HOW_MANY_NODE
sed s/PBSARRAYSIZE/$HOW_MANY_NODE/g < batch.sh > tbatch.sh
nsub tbatch.sh 
rm tbatch.sh
