#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --mem=200G
#SBATCH --job-name integration_r
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute integration_r.r.ipynb