#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name scanvi_spleno_prog_cpg
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute scanvi_spleno_prog_cpg.p.ipynb