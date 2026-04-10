#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --mem=10G
#SBATCH --job-name panel_1
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute panel_1.r.ipynb
