#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=8
#SBATCH --job-name rss
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute rss.r.ipynb