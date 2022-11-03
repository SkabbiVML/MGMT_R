#!/bin/bash

mkdir -p /home/vml/MGMT/Analysis/Results/$1/$2

megalodon /home/vml/MGMT/Analysis/Sequences/$1/$2 \
--output-directory /home/vml/MGMT/Analysis/Results/$1/$2  \
--guppy-server-path /home/vml/ont-guppy/ont-guppy/bin/guppy_basecall_server \
--ref-mods-all-motifs m 5mC CG 0 \
--remora-modified-bases dna_r10.4_e8.1 hac 0.0.0 5mc CG 0 \
--reference /home/vml/MGMT/Analysis/Chrom10_reference/Chr10_reference.fasta \
--outputs basecalls mappings mod_mappings mods \
--sort-mappings \
--mod-map-emulate-bisulfite \
--mod-map-base-conv C T \
--mod-map-base-conv m C \
--guppy-config dna_r10.4.1_e8.2_400bps_hac.cfg \
--devices 0 \
--processes 16 \
--overwrite  

mv /home/vml/MGMT/Analysis/Results/$1/$2/modified_bases.5mC.bed /home/vml/MGMT/Analysis/Results/$1/$2/"$2".5mC.bed 
mv /home/vml/MGMT/Analysis/Results/$1/$2/mod_mappings.5mC.sorted.bam /home/vml/MGMT/Analysis/Results/$1/$2/"$2".5mC.sorted.bam
mv /home/vml/MGMT/Analysis/Results/$1/$2/mod_mappings.5mC.sorted.bam.bai /home/vml/MGMT/Analysis/Results/$1/$2/"$2".5mC.sorted.bam.bai
mv /home/vml/MGMT/Analysis/Results/$1/$2/basecalls.fastq /home/vml/MGMT/Analysis/Results/$1/$2/"$2".basecalls.fastq