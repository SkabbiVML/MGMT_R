#!/bin/bash

mkdir -p /home/vml/MGMT/Analysis/Results/GuppY2Modbed/$1/$2

guppy_basecaller -i /encrypteddata/MGMT/Sequences/$1/$2 \
-s /home/vml/MGMT/Analysis/Results/GuppY2Modbed/$1/$2 \
  -x cuda:0 \
  --config dna_r9.4.1_450bps_modbases_5mc_cg_hac.cfg \
  --disable_qscore_filtering \
  --bam_out \
  --index \
  --align_ref /home/vml/MGMT/Analysis/Chrom10_reference/Chr10_reference.fasta \
  --chunks_per_runner 512


samtools merge - /home/vml/MGMT/Analysis/Results/GuppY2Modbed/$1/$2/*.bam | \
samtools sort -@ 16 -o /home/vml/MGMT/Analysis/Results/GuppY2Modbed/$1/$2/$2.sorted.bam

samtools index /home/vml/MGMT/Analysis/Results/GuppY2Modbed/$1/$2/$2.sorted.bam

modbam2bed --aggregate \
  -e \
  -m 5mC \
  --cpg \
  -t 12 \
  --prefix=/home/vml/MGMT/Analysis/Results/GuppY2Modbed/$1/$2/$2 \
  /home/vml/MGMT/Analysis/Chrom10_reference/Chr10_reference.fasta \
  /home/vml/MGMT/Analysis/Results/GuppY2Modbed/$1/$2/$2.sorted.bam > /home/vml/MGMT/Analysis/Results/GuppY2Modbed/$1/$2/$2.FULL.bed


