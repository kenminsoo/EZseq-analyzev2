#!/bin/bash

#have this file loaded in by a config option.

#usage: export var=value3
#1 = true
#0 = false

#we will build new genome for each experiment in a specific place. check if exist
export unalign_dir=/Users/kenminsoo/Desktop/raw_data/*.fastq.gz
export annotation_dir=/home/kenmn/bioinfo_tools/smrnaseq_annotation/
export alignment_tool='bowtie'
export genome_path=/home/kenmn/bioinfo_tools/hg38/hg38.fasta
export multilane='no'
export threads=12
export user_adapter1="AACTGTAGGCACCATCAAT"
export user_adapter2=''
export minlen=20
export reg=1
export sm_other=1
export lnc=1
export interest_all=1
export mirna=1
export pirna=1
