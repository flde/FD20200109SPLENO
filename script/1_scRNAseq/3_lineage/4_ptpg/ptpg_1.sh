#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name ptpg_1
#SBATCH -o %x.out
#SBATCH -e %x.err

export MC_CORES=1

jupyter nbconvert --to html --execute ptpg_1.r.ipynb