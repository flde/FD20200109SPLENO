#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --mem=200G
#SBATCH --job-name fetch_1
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute fetch_1.r.ipynb