#!/bin/bash

#SBATCH --partition=mediumq
#SBATCH --qos=mediumq
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name ery_gam
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute ery_gam.r.ipynb