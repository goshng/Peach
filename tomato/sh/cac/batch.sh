#!/bin/bash
#PBS -l walltime=36:00:00,nodes=1
#PBS -A acs4_0001
#PBS -j oe
#PBS -N tomato
#PBS -q v4
#PBS -m e
# #PBS -M schoi@cornell.edu
#PBS -t 1-100

R=/home/fs01/sc2265/Downloads/r-devel/b/bin/Rscript

function copy-data {
  cd $TMPDIR
  # RData
  cp $PBS_O_WORKDIR/cor.RData .
  cp $PBS_O_WORKDIR/correlation-total.R .

  # Create directories
  mkdir output
}

function retrieve-data {
  cp output/* $PBS_O_WORKDIR/output
}

function process-data {
  cd $TMPDIR
  CORESPERNODE=8

  PBS_ARRAYIDINDEX=$((PBS_ARRAYID - 1))
  for (( i=1; i<=CORESPERNODE; i++))
  do
    A=$((8 * $PBS_ARRAYIDINDEX + $i))   
    $R correlation-total.R 800 $A &
  done
}

copy-data
process-data; wait
retrieve-data
cd $PBS_O_WORKDIR
rm -rf $TMPDIR
