#!/bin/bash
#$ -cwd
# specify resources needed
#$ -l h_cpu=48:00:00
#$ -N fs_run
#$ -o fs.out
#$ -e fs.err
#$ -m a
#$ -q all.q

#Running fastsimcoal
./fastsimcoal -t 1PopDNArand.tpl -n 1 -e 1PopDNArand.est -E 1000 -q
