#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name haemosphere
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute haemosphere.r.ipynb