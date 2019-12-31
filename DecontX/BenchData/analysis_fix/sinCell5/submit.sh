#!/bin/bash -l

#$ -S /bin/bash

#$ -P camplab

#$ -cwd

#$ -j y

#$ -o  decon.log 

#$ -N  c5bench

#$ -m ea 

#$ -l h_rt=100:00:00

#$ -l mem_free=24g

#$ -pe omp 4

source ~/.bashrc
module load R/3.6.0
module load mpfr
module load gcc/5.1.0

Rscript 10x.R  & 
Rscript celseq2.R  & 
Rscript celseq2_p2.R &
Rscript celseq2_p3.R 

wait
