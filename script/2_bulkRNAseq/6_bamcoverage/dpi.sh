#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --cpus-per-task 32
#SBATCH --mem=200G
#SBATCH --job-name dpi
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute dpi.p.ipynb
