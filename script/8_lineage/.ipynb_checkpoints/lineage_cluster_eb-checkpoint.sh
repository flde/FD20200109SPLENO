#!/bin/bash

#SBATCH --partition=mediumq
#SBATCH --qos=mediumq
#SBATCH --mem-per-cpu=MaxMemPerCPU
#SBATCH --cpus-per-task 8
#SBATCH --job-name lineage_cluster_eb
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute lineage_cluster_eb.r.ipynb