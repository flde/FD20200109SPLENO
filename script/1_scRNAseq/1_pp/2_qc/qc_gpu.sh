#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:h100pcie:1
#SBATCH --mem=200G
#SBATCH --cpus-per-task 8
#SBATCH --time=2:00:00
#SBATCH --job-name qc
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute qc.p.ipynb --output qc.p