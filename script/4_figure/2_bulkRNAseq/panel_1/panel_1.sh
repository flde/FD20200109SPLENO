#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --job-name panel_1
#SBATCH --mem=200G
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute panel_1.r.ipynb