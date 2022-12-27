#!/bin/bash

#SBATCH --job-name=
#SBATCH --partition=shared
#SBATCH --time=2-00:00:00

#SBATCH --cpus-per-task=12
#SBATCH --mem=128G

#SBATCH --error=sm_process%A.err
#SBATCH --output=sm_process%A.out

#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --mail-user=knakats@emory.edu

module load lang/Anaconda3/2022.05
module load bio/SAMtools/1.9-intel-2018.5.274
eval "$(conda shell.bash hook)"
conda activate smRNAseq

for f in *.bam; do
file="$f"
sampname="${file%%.sa*}"
umi_tools dedup -I $f --output-stats="deduplicated"_$sampname -S "dedup"_"$f"
done
