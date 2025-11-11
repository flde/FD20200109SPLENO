#!/bin/bash

#SBATCH --partition=mediumq
#SBATCH --qos=mediumq
#SBATCH --mem=200G
#SBATCH --cpus-per-task 8
#SBATCH --job-name fetch
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute fetch.r.ipynb
