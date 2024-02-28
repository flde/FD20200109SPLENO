#!/bin/bash

#SBATCH --partition=mediumq
#SBATCH --qos=mediumq
#SBATCH --mem=200G
#SBATCH --job-name dea_cluster
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute dea_cluster.r.ipynb