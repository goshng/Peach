PBS_ARRAYID=1
CORESPERNODE=8
PBS_ARRAYIDINDEX=$((PBS_ARRAYID - 1))
for (( i=1; i<=CORESPERNODE; i++))
do
  A=$((8 * $PBS_ARRAYIDINDEX + $i))   
  echo $A
done
