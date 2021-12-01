#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --mem=100G
#SBATCH --job-name cluster_deg
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute cluster_deg.r.ipynb