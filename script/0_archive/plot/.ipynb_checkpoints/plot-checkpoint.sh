#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --mem=200G
#SBATCH --job-name plot
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute plot.r.ipynb