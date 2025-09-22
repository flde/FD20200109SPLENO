#!/bin/bash

#SBATCH --partition=mediumq
#SBATCH --qos=mediumq
#SBATCH --mem=200G
#SBATCH --cpus-per-task 8
#SBATCH --job-name 1_pseudotime
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute 1_pseudotime.r.ipynb