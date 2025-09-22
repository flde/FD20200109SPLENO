#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --mem=200G
#SBATCH --job-name run_2
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute run_2.r.ipynb