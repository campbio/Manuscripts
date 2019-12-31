#!/bin/bash -l

#$ -S /bin/bash

#$ -P camplab

#$ -cwd

#$ -j y

#$ -o  decon.log 

#$ -N  c3bench

#$ -m ea 

#$ -l h_rt=100:00:00

#$ -l mem_free=24g

#$ -pe omp 3

source ~/.bashrc
module load R/3.6.0
module load mpfr
module load gcc/5.1.0

Rscript 10x.R  & 
Rscript celseq2.R  & 
Rscript dropseq.R

wait
