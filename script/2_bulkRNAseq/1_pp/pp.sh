#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --job-name pp
#SBATCH --mem=200G
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute pp.r.ipynb