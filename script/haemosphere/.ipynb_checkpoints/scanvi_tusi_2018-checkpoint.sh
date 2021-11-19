#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name scanvi_tusi_2018
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute scanvi_tusi_2018.p.ipynb