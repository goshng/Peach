#!/bin/bash
#PBS -l walltime=236:00:00,nodes=1
#PBS -A acs4_0001
#PBS -j oe
#PBS -N phyml 
#PBS -q v4
#PBS -m e
# #PBS -M schoi@cornell.edu
#PBS -t 1-PBSARRAYSIZE

PHYML=$HOME/usr/bin/phyml

function copy-data {
  cd $TMPDIR
  # Copy excutables
  cp $PHYML .
  # Copy shell scripts
  cp $PBS_O_WORKDIR/batchjob.sh .
  # Create directories
  mkdir output
  mkdir log
  cp -r $PBS_O_WORKDIR/input .
  mkdir $PBS_O_WORKDIR/status/$PBS_ARRAYID
}

function retrieve-data {
  cp output/*_phyml_* $PBS_O_WORKDIR/output
  # Remove the status directory.
  rm -rf $PBS_O_WORKDIR/status/$PBS_ARRAYID
}

function process-data {
  cd $TMPDIR
  CORESPERNODE=8
  for (( i=1; i<=CORESPERNODE; i++))
  do
    bash batchjob.sh \
      $i \
      $PBS_O_WORKDIR/jobidfile \
      $PBS_O_WORKDIR/lockfile \
      $PBS_O_WORKDIR/status/$PBS_ARRAYID \
      PBSARRAYSIZE&
  done
}

copy-data
process-data; wait
retrieve-data
cd $PBS_O_WORKDIR
rm -rf $TMPDIR
