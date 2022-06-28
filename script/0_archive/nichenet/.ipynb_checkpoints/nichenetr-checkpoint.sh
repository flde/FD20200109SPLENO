#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=100G
#SBATCH --job-name nichenetr
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute nichenetr.r.ipynb