#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --job-name BSF_1609_H57CKDRX5_2_samples
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH -o %x.out
#SBATCH -e %x.err

module load picard

IN_DIR="/nobackup/lab_bsf/samples/BSF_1609_H57CKDRX5/BSF_1609_H57CKDRX5_2_samples"
OUT_DIR="/nobackup/peer/fdeckert/FD20200109SPLENO/geo/bulkRNAseq/"

for bam_file in "${IN_DIR}"/*FD_*.bam; do

    sample_name=$(basename "${bam_file}" .bam)
    sample_name=${sample_name#BSF_1609_H57CKDRX5_2#}
    sample_name=${sample_name%%_S*}

    java -jar $EBROOTPICARD/picard.jar SamToFastq --INPUT "${bam_file}" --FASTQ "${OUT_DIR}/${sample_name}_R2.fastq.gz"

done