#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name sct_int
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute sct_int.r.ipynb
