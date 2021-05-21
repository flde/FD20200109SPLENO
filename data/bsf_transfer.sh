#!/bin/bash

#SBATCH --partition=longq
#SBATCH --qos=longq
#SBATCH --job-name bsf_transfer
#SBATCH -o %j.out
#SBATCH -e %j.err

cp /nobackup/lab_bsf/projects/BSA_0501_SK_RG_LEUKO /research/peer/fdeckert/FD20210520LEUKO/data/