#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --mem=10G
#SBATCH --job-name stat
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute stat.r.ipynb
