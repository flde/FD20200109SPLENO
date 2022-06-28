#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --mem=200G
#SBATCH --job-name scanorama_original_API
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute scanorama_original_API.p.ipynb