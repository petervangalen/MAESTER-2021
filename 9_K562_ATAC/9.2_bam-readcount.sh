# Peter van Galen, 210704
# To call variant allele frequences from K562 ATAC-seq, Leif performed standard bulk ATAC-seq
# After demultiplexing, I used STAR for alignment and PicardTools for deduplication, as you can see in 9.1_Align_and_deduplicate_edit.sh
# I copied the resulting bam file for local processing using bam-readcount and R

#### INSTALL bam-readcount ###
cd ~
conda create --name py2 python=2.7
conda activate py2
conda install -c bioconda bam-readcount
# That worked.


#### RUN BAM-READCOUNT ####
cd ~/DropboxMGB/Projects/Maester/AnalysisPeter/K562_ATAC

# Run bam-readcount with standard options. I added -w to prevent a warning for every read. bam-readcount does not take into account paired reads; reads from the same fragment will be counted independently.
bam-readcount -w 5 -f MT.fa K562_ATAC_Aligned.sortedByCoord.out.MT.dedup.bam > bam-readcount_results.txt
