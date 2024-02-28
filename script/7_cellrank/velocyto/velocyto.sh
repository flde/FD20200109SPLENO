#!/bin/bash

# Velocyto run for spliced/unspliced loom files from 10x data 
# Run conda activate p.3.6.13-FD20200109SPLENO to activate velocyto 
# Run as array job sbatch --array=0-7 velocyto.sh

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name velocyto
#SBATCH --cpus-per-task 48
#SBATCH -o %x_%a.out
#SBATCH -e %x_%a.err

# Load modules for velocyto 
module load SAMtools

# Root directory
cd /research/peer/fdeckert/FD20200109SPLENO/

# Reference files 
ref_gtf=$"data/annotation/mm10/cellranger/mm10-2020-A_build/gencode.vM23.primary_assembly.annotation.gtf.filtered"

# 10x sample folder array 
sample_dir=($(ls -d data/BSA_0355_SM01_10x_SPLENO/OUT/COUNT/*))

# Select 10x sample folder based on array job id
sample_dir=${sample_dir[$SLURM_ARRAY_TASK_ID]}

barcode_tsv=$''$sample_dir'/filtered_feature_bc_matrix/barcodes.tsv.gz'
barcode_bam=$''$sample_dir'/possorted_genome_bam.bam'
loom_file=$'velocyto'

velocyto run -l "Permissive10X" -@ 48 --samtools-memory 100 -b $barcode_tsv -o $sample_dir -e $loom_file $barcode_bam $ref_gtf