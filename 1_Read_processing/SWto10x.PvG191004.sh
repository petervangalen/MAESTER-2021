#!/bin/bash

# Peter van Galen, 191004
# Move cell barcode (CB) and unique molecular identifier (UMI) from read identifier to sam tags.
# Resulting bam file will have tags for cell barcode (CB) and UMI (UB) as per 10X convention https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam)
# These tags are required to run MAEGATK software

# Example execution:
#cd /broad/vangalenlab/vangalen/SW_pipeline/4_star
#qsub -N OCIAML3-10X -l h_rt=4:00:00,h_vmem=16g -j y -b y -cwd /broad/vangalenlab/vangalen/packages.sh /broad/vangalenlab/vangalen/SW_pipeline/scripts/SWto10X.PvG190524.sh /broad/vangalenlab/vangalen/SW_pipeline/4_star/190606.190418.OCIAML3.star/190606.190418.OCIAML3_Aligned.toTranscriptome.out.bam

source /broad/software/scripts/useuse
reuse -q Samtools

# First variable is bam file to convert
INPUT=$1
# Second variable is bam file to write
OUTPUT="$(echo "${INPUT/bam/10x.bam}")"

echo "Converting $INPUT into $OUTPUT..."

samtools view -h $INPUT | awk 'BEGIN{FS="\t"; OFS="\t"} {
	if (substr($1,1,1) == "@") {
		print $0
	} else {
		split($1, a, "_")
		$1=""
		print a[1]"_"a[2]$0"\tCB:Z:"a[3]"-1\tUB:Z:"a[4]
	} }' | samtools view -bh > $OUTPUT

echo "Done!"

date

exit 0