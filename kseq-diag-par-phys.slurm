#!/bin/bash
#
#SBATCH --job-name=kseq-diag-par-phys
#SBATCH --output=diag-par-phys-res.txt
#
#SBATCH --partition=physical
#SBATCH --constraint=physg5
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=36G


module load gcc/10.1.0
make
./kseqalign < mseq-big13-example.dat
