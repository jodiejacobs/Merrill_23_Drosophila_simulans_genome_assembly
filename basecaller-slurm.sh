#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --mail-user=jomojaco@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --exclude=phoenix-09
#SBATCH --mem=100gb
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=25
#SBATCH --output=logs/%x.%j.log
#SBATCH --time=06:00:00

set -o pipefail
set -e
set -u
set -o xtrace

POD5DIR=$1
OUTDIR=$2
REF=$3

# Ensure output directory exists
mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR/demux"

# Basecalling with alignment
/private/home/jomojaco/dorado-0.7.3-linux-x64/bin/dorado basecaller hac $POD5DIR \
    --device cuda:all \
    --reference $REF \
    --kit-name SQK-NBD114-24 \
    > $OUTDIR/aligned_basecalled.bam

# Demux the specific BAM file
/private/home/jomojaco/dorado-0.7.3-linux-x64/bin/dorado demux \
    "$OUTDIR/aligned_basecalled.bam" \
    --output-dir "$OUTDIR/demux" \
    --kit-name SQK-NBD114-24