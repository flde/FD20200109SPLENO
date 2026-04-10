#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --job-name BSF_1450_HGVJNDRX3_GEO
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH -o %x.out
#SBATCH -e %x.err

set -euo pipefail

GEO_USER="geoftp"
GEO_PASS="inAlwokhodAbnib5"
GEO_HOST="ftp-private.ncbi.nlm.nih.gov"

WORKSPACE="uploads/theknapplab_8xSZQMQq"
LOCAL_DIR="/nobackup/peer/fdeckert/FD20230822MAVI/data/GEO/bulkRNAseq/"
REMOTE_DIR="bulkRNAseq"

lftp -u "${GEO_USER},${GEO_PASS}" ftp://${GEO_HOST} <<EOF
set ftp:ssl-allow yes
set ftp:passive-mode yes
set net:max-retries 20
set net:reconnect-interval-base 15
set net:reconnect-interval-max 60
set xfer:parallel 4

cd ${WORKSPACE}
mkdir -p ${REMOTE_DIR}
mirror -R --continue --verbose ${LOCAL_DIR} ${REMOTE_DIR}
bye
EOF