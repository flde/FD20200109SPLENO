#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name cluster_clean
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute cluster_clean.r.ipynb