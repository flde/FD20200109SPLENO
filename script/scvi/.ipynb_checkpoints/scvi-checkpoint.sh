#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name scvi_hvg_catc
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute scvi.p.ipynb --output scvi_hvg_catc