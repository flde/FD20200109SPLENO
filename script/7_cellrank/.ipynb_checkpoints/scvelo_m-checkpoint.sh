#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --cpus-per-task 48
#SBATCH --job-name scvelo_m
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute scvelo_m.p.ipynb