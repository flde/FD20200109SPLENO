#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --mem=10G
#SBATCH --job-name panel_2
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute panel_2.r.ipynb
