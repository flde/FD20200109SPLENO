#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name raceid_myeloid
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute raceid_myeloid.r.ipynb