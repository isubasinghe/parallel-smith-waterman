#!/bin/bash
#
#SBATCH --job-name=kseq-diag-par-phys
#SBATCH --output=diag-par-phys-res.txt
#
#SBATCH --partition=physical
#SBATCH --constraint=physg5
#SBATCH --cores-per-socket=18
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=18
#SBATCH --mem=36G


module load gcc/10.1.0
make
./kseqalign < mseq-big13-example.dat
