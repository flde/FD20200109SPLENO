#!/bin/bash

#SBATCH --partition=tinyq
#SBATCH --qos=tinyq
#SBATCH --mem=200G
#SBATCH --cpus-per-task 8
#SBATCH --job-name transfer
#SBATCH -o %x.out
#SBATCH -e %x.err

cp -r /research/peer/fdeckert/FD20200109SPLENO/data/scRNAseq/cellranger/BSA_0355_SM01_10x_SPLENO/M6_control_Prog/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
cp -r /research/peer/fdeckert/FD20200109SPLENO/data/scRNAseq/cellranger/BSA_0355_SM01_10x_SPLENO/M8_control_Prog/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
cp -r /research/peer/fdeckert/FD20200109SPLENO/data/scRNAseq/cellranger/BSA_0355_SM01_10x_SPLENO/M1_CpG_Prog/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
cp -r /research/peer/fdeckert/FD20200109SPLENO/data/scRNAseq/cellranger/BSA_0355_SM01_10x_SPLENO/M2_CpG_Prog/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
cp -r /research/peer/fdeckert/FD20200109SPLENO/data/scRNAseq/cellranger/BSA_0355_SM01_10x_SPLENO/M6_control_Mac/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
cp -r /research/peer/fdeckert/FD20200109SPLENO/data/scRNAseq/cellranger/BSA_0355_SM01_10x_SPLENO/M8_control_Mac/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
cp -r /research/peer/fdeckert/FD20200109SPLENO/data/scRNAseq/cellranger/BSA_0355_SM01_10x_SPLENO/M1_CpG_Mac/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
cp -r /research/peer/fdeckert/FD20200109SPLENO/data/scRNAseq/cellranger/BSA_0355_SM01_10x_SPLENO/M2_CpG_Mac/outs/filtered_feature_bc_matrix/barcodes.tsv.gz