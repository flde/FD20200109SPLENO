#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=50G
#SBATCH --job-name singler
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute singler.r.ipynb