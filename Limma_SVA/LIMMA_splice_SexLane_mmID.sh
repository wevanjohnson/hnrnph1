#!/bin/bash
#$ -S /bin/bash
#$ -N LIMMA_run
#$ -cwd
#$ -j Y
#$ -V
#$ -m be
#$ -M reeder@bu.edu
#$ -l h_rt=120:00:00

module load R/R-3.1.1
R CMD BATCH --vanilla LIMMA_splice_SexLane_mmID.R