#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name raceid_treatment_batch
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute raceid_treatment_batch.r.ipynb