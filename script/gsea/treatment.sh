#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name treatment_mast_scvi
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute treatment.r.ipynb --output treatment_mast_scvi.r