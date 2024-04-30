#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=4G
#SBATCH -t 0:80:00
#SBATCH -p priority
#SBATCH -o ultima_cram_split_fq_%A.out

samtools index -@20 $1.cram
samtools view -@20 -c $1.cram
samtools flagstat -@20 $1.cram
samtools fastq -@20 $1.cram > $1.fastq
pigz $1.fastq
samtools view -@20 -c $1.fastq.gz