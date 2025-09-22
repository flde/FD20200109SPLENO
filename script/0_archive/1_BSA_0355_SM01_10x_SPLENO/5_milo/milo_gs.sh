#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name milo_gs
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute milo_gs.r.ipynb