#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name velocyto_qc
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute velocyto_qc.r.ipynb

