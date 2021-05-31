#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name sct_int
#SBATCH -o %x.out
#SBATCH -e %x.err

module load Pandoc
R -e "rmarkdown::render('norm_clust.Rmd',output_file='norm_clust.html')"