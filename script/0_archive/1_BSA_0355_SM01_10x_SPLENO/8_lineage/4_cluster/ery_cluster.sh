#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --cpus-per-task 8
#SBATCH --job-name ery_cluster
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute ery_cluster.r.ipynb