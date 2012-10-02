#!/bin/bash
#PBS -l walltime=136:00:00,nodes=1
#PBS -A acs4_0001
#PBS -j oe
#PBS -N tomato
#PBS -q v4
#PBS -m e
# #PBS -M schoi@cornell.edu

R=/home/fs01/sc2265/Downloads/r-devel/b/bin/Rscript
cd $PBS_O_WORKDIR
$R correlation-total-adjust.R
rm -rf $TMPDIR
