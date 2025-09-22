#!/bin/bash

#SBATCH --partition=mediumq
#SBATCH --qos=mediumq
#SBATCH --mem=200G
#SBATCH --cpus-per-task 8
#SBATCH --job-name ery_dpt
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute ery_dpt.p.ipynb