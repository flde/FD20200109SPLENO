#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=32  # Number of CPU cores per task (match with --num_workers)
#SBATCH --job-name grn
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute grn.p.ipynb