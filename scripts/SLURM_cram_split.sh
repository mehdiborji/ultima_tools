#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=2G
#SBATCH -t 0:30:00
#SBATCH -p short
#SBATCH --account=chen_fec176
#SBATCH -o ultima_cram_split_fq_%A.out

echo 'indir =' $1
echo 'sample =' $2

python ~/ultima_demux/cram_pipeline.py -c 20 -i $1 -s $2
#--limit
