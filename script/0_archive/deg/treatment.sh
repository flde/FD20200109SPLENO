#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --mem=100G
#SBATCH --job-name treatment
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute treatment.r.ipynb