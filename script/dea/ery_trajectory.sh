#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name ery_trajectory
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute ery_trajectory.r.ipynb