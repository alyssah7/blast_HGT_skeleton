#!/bin/bash
#SBATCH --job-name widevariant.main
#SBATCH -p sched_mit_tami
#SBATCH -n 1
#SBATCH --time=2-00:00:00
#SBATCH --mem=128GB
#SBATCH -o mainout.txt
#SBATCH -e mainerr.txt
#SBATCH --mail-user=alyssah@mit.edu
#SBATCH --mail-type=ALL

# Activate conda environment (may need to change name of env)
source activate snakemake

snakemake --profile .

# Print "Done!!!" at end of main log file
echo Done!!!
