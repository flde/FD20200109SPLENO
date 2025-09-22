#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:h100pcie:2
#SBATCH --mem=200G
#SBATCH --cpus-per-task 8
#SBATCH --time=0:30:00
#SBATCH --job-name fetch
#SBATCH -o %x.out
#SBATCH -e %x.err

jupyter nbconvert --to html --execute fetch.p.ipynb --output fetch.p