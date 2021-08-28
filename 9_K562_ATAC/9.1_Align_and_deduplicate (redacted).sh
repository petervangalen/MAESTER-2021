# Peter van Galen, 210827
# Align K562 ATAC fastq for MAESTER revision
# This was run on the Broad cluster

# Fastqs are available on GEO. Alignment was performed with STAR version 2.7.2a as follows:

STAR \
--readFilesIn K562_bulk_ATAC_R1.fastq.gz K562_bulk_ATAC_R2.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix K562_ATAC_ \
--genomeDir ~/GRCh38 \
--runThreadN 4 \
--genomeLoad NoSharedMemory \
--outSAMtype BAM SortedByCoordinate

# Subset for chrM (script is available upon request, it basically applies "samtools view $INPUT $REF -b > $OUTPUT")
./SubsetChr.PvG191004.sh K562_ATAC_Aligned.sortedByCoord.out.bam MT

# Remove duplicates
use Picard-Tools
java -Xmx2g -jar /seq/software/picard/1.802/bin/MarkDuplicates.jar \
REMOVE_DUPLICATES=true \
INPUT=K562_ATAC_Aligned.sortedByCoord.out.MT.bam \
OUTPUT=K562_ATAC_Aligned.sortedByCoord.out.MT.dedup.bam \
METRICS_FILE=picard_dedup_metrics.txt

# Save MT reference for bam-readcount
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa MT > MT.fa
