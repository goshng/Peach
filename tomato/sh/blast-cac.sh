#!/bin/bash
#PBS -l walltime=12:00:00,nodes=1
#PBS -A acs4_0001
#PBS -j oe
#PBS -N sm-BLAST
#PBS -q v4
#PBS -m e
#PBS -t 1-100

NUM=$(printf "%03d" $PBS_ARRAYID)
cd $TMPDIR
cp $PBS_O_WORKDIR/input/blastp .
mkdir -p output/tomato/1/data/t
mkdir -p output/tomato/1/run-analysis
cp $PBS_O_WORKDIR/input/output/tomato/1/data/uniref90.* output/tomato/1/data
cp $PBS_O_WORKDIR/input/output/tomato/1/data/t/frag$NUM output/tomato/1/data/t

./blastp -db output/tomato/1/data/uniref90 -query output/tomato/1/data/t/frag$NUM \
    -task blastp \
    -outfmt 6 \
    -num_threads 8 \
    -evalue 1e-3 \
    -out output/tomato/1/run-analysis/goseq$NUM.blast

cp output/tomato/1/run-analysis/goseq$NUM.blast $PBS_O_WORKDIR/output

cd
rm -rf $TMPDIR
