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

#this script will do:
# 1) run a quality check
# 2) remove adapters based upon a given adapter seq
# 3) run a post-trim quality check
# 4) run alignment with bowtie 2 against hg38
# 5) run annotation with four different files
# 5.1) all, all_smrna, mirna, pirna
#start the script in a dir with two files, annotation and unaligned files with ref. genome
#ensure that you have conda installed with all packages

export unalign_dir=/Users/kenminsoo/Desktop/raw_data/*.fastq.gz
export annotation_dir=/home/kenmn/bioinfo_tools/smrnaseq_annotation/
export alignment_tool='bowtie'
export genome_path=/home/kenmn/bioinfo_tools/hg38/hg38.fasta
export genome_name=hg38
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

module load lang/Anaconda3/2022.05

eval "$(conda shell.bash hook)"
#conda activate ezanalyze
conda activate smRNAseq

bash /home/kenmn/scripts/sm_rnaseq/raw_process/config.sh
mkdir unaligned

mv $unalign_dir unaligned

mkdir output
#all output files will be stored here
mkdir output/fastqc_analysis
mkdir output/fastqc_analysis/unprocessed
#fastq
mkdir output/fastqc_analysis/docs
#zips and logs
mkdir output/fastqc_analysis/results
#html files

#run fastqc with 8 threads. move to respective directory
fastqc unaligned/*.fastq.gz -t $threads -o output/fastqc_analysis
mv output/fastqc_analysis/*.html output/fastqc_analysis/results
mv output/fastqc_analysis/*.zip output/fastqc_analysis/docs

#user can input a new adapter or add another one with --adapter2 '$user_adapter2'
#this will remove adapters and move it to its respective dir

#add UMIs
mkdir UMI_fastq
cd unaligned
for f in *.fastq.gz ; do
	file="$f"
	sampname="${file%%.*}"
	umi_tools extract --extract-method=regex \
                  --bc-pattern=".+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12}).+" \
                   -I $file \
                   -S "processed".${sampname}".fastq.gz"
done

cd ..
mv unaligned processed* UMI_fastq


cd UMI_fastq

#var%%.x will remove everything after . i.e. y.x -> y

for f in *.fastq.gz ; do
    file="$f"
    sampname="${file%%.*}"
    AdapterRemoval --file1 "${file}" --adapter1 $user_adapter1 --basename "trimmed"_"${sampname}" --gzip --threads $threads --trimqualities --minlength $minlen
done

cd ..

mkdir output/adaptertrim
mkdir output/adaptertrim/discarded_ar
mkdir output/adaptertrim/settings

mv unaligned/*.discarded.gz output/adaptertrim/discarded_ar
mv unaligned/*.settings output/adaptertrim/settings
mv unaligned/*.truncated.gz output/adaptertrim

#change into fastq file
#note: ${var//x/y} will take the variable, remove x from it, and replace it with y
for f in output/adaptertrim/*.gz ; do
    mv $f ${f//truncated/fastq};
done
#passes ~800,000-900,000 reads per second

#run post trimming QC
mkdir output/fastqc_post_trim
fastqc output/adaptertrim/*.fastq.gz -t $threads -o output/fastqc_post_trim

#now lets align and move those aligned files to a new folder

mkdir output/aligned

bowtie-build --threads $threads $genome_dir "unaligned/"$genome_name

for f in output/adaptertrim/*.fastq.gz; do
    file="$f"
    sampname="${file%%.*}"
    bowtie -x 'unaligned/'$genome_name $f -v 1 -S "${sampname}".sam -p $threads --reorder
done

for f in output/adaptertrim/*.sam; do
    mv $f output/aligned
done

#run featureCounts for each one of the annotation files

if (($reg==1)) #genes
then
    featureCounts -T $threads -a $annotation_dir"hsa_all.gff3" -F 'GTF' -g 'gene_name' -o "hsa_all.tsv" output/aligned/*.sam -O -M --fraction
fi

if (($interest_all==1)) #wide
then
    featureCounts -T $threads -a $annotation_dir"hsa_alt.gff" -F 'GTF' -t 'mir-3p,mir-5p,miRNAprimary,mRNA,mRNA_exon,mRNA_intron,piRNA,rRNA,scRNA,snoRNA,snRNA,tRF3,tRF5,tRNA,UTR3_exon,UTR3_intron,UTR5_exon,UTR5_intron' -g 'ID'  -o "hsa_most.tsv" output/aligned/*.sam -O -M --fraction
fi

if (($lnc==1)) #lncrna
then
    featureCounts -T $threads -a $annotation_dir"hsa_alt.gff" -F 'GTF' -t 'lnc_RNA' -g 'ID'  -o "hsa_lnc.tsv" output/aligned/*.sam -O -M --fraction
fi

if (($mirna==1)) #mirna
then
    featureCounts -T $threads -a $annotation_dir"hsa_mirna.gff3" -F 'GTF' -t 'miRNA' -g 'Name' -o "hsa_mirna.tsv" output/aligned/*.sam -s 0 -O -M --fraction
fi

if (($pirna==1)) #pirna
then
    featureCounts -T $threads -a $annotation_dir"hsa_pirna.gtf" -F 'GTF' -t 'piRNA' -g 'piRNA_code' -o "hsa_pirna.tsv" output/aligned/*.sam -O -M --fraction
fi

if (($sm_other==1)) #other ncRNA
then
    featureCounts -T $threads -a $annotation_dir"hsa_smrna.gff" -F 'GTF' -t 'snRNA' -g 'ID' -o 'hsa_snrna.tsv' output/aligned/*.sam -O -M --fraction
fi

if (($sm_other==1)) #other ncRNA
then
    featureCounts -T $threads -a $annotation_dir"hsa_smrna.gff" -F 'GTF' -t 'tRF3,tRF5,tRNA' -g 'ID' -o 'hsa_trna.tsv' output/aligned/*.sam -O -M --fraction
fi

if (($sm_other==1)) #other ncRNA
then
    featureCounts -T $threads -a $annotation_dir"hsa_smrna.gff" -F 'GTF' -t 'snoRNA' -g 'ID' -o 'hsa_snorna.tsv' output/aligned/*.sam -O -M --fraction
fi
    
mkdir rawcount-analysis
    
mkdir output/counts
mv *.tsv output/counts
for f in *.summary ; do
    mv $f ${f//.tsv.summary/.summary.tsv};
done
mv *.tsv output/counts
mkdir output/counts/summary
mv output/counts/*.summary.tsv output/counts/summary

mkdir rawcount-analysis/counts
mv output/counts/*.tsv rawcount-analysis/counts

touch output/counts/note.txt
echo 'raw counts are now in the rawcount-analysis/counts folder' > output/counts/note.txt


