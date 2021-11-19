#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name dim_clust_hvg_conc
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute dim_clust.r.ipynb --output dim_clust_hvg_conc