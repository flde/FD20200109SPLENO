#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --cpus-per-task 1
#SBATCH --job-name seurat_cluster_m
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute seurat_cluster_m.r.ipynb