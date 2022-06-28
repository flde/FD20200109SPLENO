#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --cpus-per-task 32
#SBATCH --job-name cellrank_cpg
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute cellrank.p.ipynb --output cellrank_cpg