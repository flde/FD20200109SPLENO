#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --job-name solo
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task 38
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute solo.p.ipynb