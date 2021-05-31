#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --job-name bsf_transfer
#SBATCH -o %x.out
#SBATCH -e %x.err

cp -r /nobackup/lab_bsf/projects/BSA_0355_SM01_10x_SPLENO /research/peer/fdeckert/FD20200109SPLENO/data/