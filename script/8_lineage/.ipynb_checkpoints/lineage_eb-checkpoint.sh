#!/bin/bash

#SBATCH --partition=mediumq
#SBATCH --qos=mediumq
#SBATCH --mem=50G
#SBATCH --cpus-per-task 1
#SBATCH --job-name lineage_eb
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute lineage_eb.r.ipynb