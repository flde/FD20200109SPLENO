#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --mem=10G
#SBATCH --cpus-per-task 8
#SBATCH --job-name annotation
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute annotation.p.ipynb --output annotation.p