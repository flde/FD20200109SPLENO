#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name gam_ery
#SBATCH -o %x.out
#SBATCH -e %x.err

export MC_CORES=1

jupyter nbconvert --to html --execute gam_ery.r.ipynb