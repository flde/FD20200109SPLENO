#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --mem=200G
#SBATCH --job-name fetch_2
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute fetch_2.r.ipynb