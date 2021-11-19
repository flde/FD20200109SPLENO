#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name raceid_progenitor
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute raceid_progenitor.r.ipynb