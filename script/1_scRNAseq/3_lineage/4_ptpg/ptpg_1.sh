#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --cpus-per-task 8
#SBATCH --job-name ptpg_1
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute ptpg_1.r.ipynb