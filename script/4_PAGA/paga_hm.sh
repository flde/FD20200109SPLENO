#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name paga_hm
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute paga_hm.r.ipynb