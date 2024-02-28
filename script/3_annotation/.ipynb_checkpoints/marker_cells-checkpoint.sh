#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name marker_cells
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute marker_cells.r.ipynb