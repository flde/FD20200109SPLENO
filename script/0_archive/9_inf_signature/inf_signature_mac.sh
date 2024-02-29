#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --cpus-per-task 1
#SBATCH --job-name inf_signature_mac
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute inf_signature_mac.r.ipynb