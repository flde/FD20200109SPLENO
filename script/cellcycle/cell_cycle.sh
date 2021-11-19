#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --mem=100G
#SBATCH --job-name cell_cycle
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute cell_cycle.r.ipynb