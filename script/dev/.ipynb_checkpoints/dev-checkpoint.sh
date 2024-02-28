#!/bin/bash

#SBATCH --partition=mediumq
#SBATCH --qos=mediumq
#SBATCH --mem=200G
#SBATCH --cpus-per-task 12
#SBATCH --job-name dev
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute dev.r.ipynb