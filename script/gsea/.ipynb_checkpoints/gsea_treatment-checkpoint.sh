#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name gsea_treatment_wilcox
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute gsea_treatment.r.ipynb --output gsea_treatment_wilcox