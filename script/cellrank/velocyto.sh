#!/bin/bash

# Velocyto run for spliced/unspliced loom files from 10x data 
# Run conda activate p.3.6.13-FD20200109SPLENO to activate velocyto 
# Run as array job sbatch --array=1-8 velocyto.sh

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --mem=200G
#SBATCH --job-name velocyto
#SBATCH -o %x_%a.out
#SBATCH -e %x_%a.err

# Load modules for velocyto 
module load SAMtools

# Root directory 
cd /research/peer/fdeckert/FD20200109SPLENO/

# Reference files 
mask_gtf=$'data/mm10/mm10_rmsk.gtf'
ref_gtf=$'data/mm10/cellranger/Mus_musculus.GRCm38.93.filtered.gtf'

# 10x sample folder array 
sample_dir=($(ls -d data/BSA_0355_SM01_10x_SPLENO/OUT/COUNT/*))

# Select 10x sample folder based on array job id
sample_dir=${sample_dir[$SLURM_ARRAY_TASK_ID]}

barcode_tsv=$''$sample_dir'/filtered_feature_bc_matrix/barcodes.qc.tsv'
barcode_bam=$''$sample_dir'/possorted_genome_bam.bam'
loom_file=$'velocyto'

velocyto run -l "ValidatedIntrons10X" -@ 10 --samtools-memory 100 -b $barcode_tsv -o $sample_dir -e $loom_file -m $mask_gtf $barcode_bam $ref_gtf