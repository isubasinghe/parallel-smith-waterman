
#!/bin/bash
#
#SBATCH --job-name=kseq-sin-phys
#SBATCH --output=diag-sin-res.txt
#
#SBATCH --partition=physical
#SBATCH --constraint=physg5
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=32G


module load gcc/10.3.0
make
./kseqalign < mseq-big13-example.dat
