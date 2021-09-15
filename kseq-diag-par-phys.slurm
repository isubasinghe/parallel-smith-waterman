#!/bin/bash
#
#SBATCH --job-name=kseq-diag-par-phys
#SBATCH --output=diag-par-phys-res.txt
#
#SBATCH --partition=interactive
#SBATCH --constraint=physg5
#SBATCH --cpus-per-task=22
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem=32G


module load gcc/10.3.0
make
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=8
./kseqalign < mseq-big13-example.dat
export OMP_NUM_THREADS=16
./kseqalign < mseq-big13-example.dat
export OMP_NUM_THREADS=22
./kseqalign < mseq-big13-example.dat
export OMP_NUM_THREADS=8
./kseqalign < mseq-big13-example.dat
export OMP_NUM_THREADS=16
./kseqalign < mseq-big13-example.dat
export OMP_NUM_THREADS=22
./kseqalign < mseq-big13-example.dat
