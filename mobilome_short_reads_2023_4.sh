# The aim of this script is to run eccFinder on mobilome-seq short read samples and 

# actiave eccFinder enviroment
conda activate ecc_finder (version 2/11/2021)

# mapping of reads
minimap2 -x sr -d reference.sr.idx reference.fa -I 4G

# index file
samtools idxstats ecc.sr.sorted.bam

# run eccfinder
python ecc_finder.py map-sr xxx.idx xxxR1.fastq.gz xxxR2.fastq.gz -r xxx.fasta  --aligner minimap2

# generate bedfile from csv output where we have position of each read
awk '{ print $1,$2,$3,$4}' *.csv > ecc_sites_185M.bed

# concatenate discordant and split read files, sum number of splir reads and discordant reads per position and select regions with atleast 5 reads coverage
cat *num.bed | awk -v OFS="\t" '{a[$1"\t"$2"\t"$3]+=$4}END{for(i in a) print i,a[i]}' | sort | uniq | awk -v OFS="\t" '($4>=5){print $1,$2,$3}' | bedtools sort | uniq > ecc.sr.disc.split.num.bed

# select regions with atleast 5 discordant/split reads and count total number of read mapping, and select locus with atleast 20 reads coverage and of atleast 100 bp in size
bedtools intersect -a ecc.sr.disc.split.num.bed -b align*/ecc.sr.sorted_185M.bed  -wa -wb | bedtools groupby -g 1,2,3 -c 7 -o count_distinct |  awk -v OFS="\t" '($4>=20 && ($3-$2)>=100){print $0}' > ecc.sr.disc.split.num.rd.bed

# intersection with TE annotation
bedtools intersect -a ecc.sr.disc.split.num.rd.bed -b '/home/lgdp/Desktop/FsylCur4_refTEs_wclassif_wreliable.frags.no_microsat.bed' -wa -wb -f 0.2 -r | awk -v OFS="\t" ' ($5!=".") {print $0}' | sort |  uniq  | grep -v "Unclassified"   | awk -v OFS="\t" '{print $1,$2,$3,$4,$8}' | bedtools sort > ecc.sr.sorted_185M-TEs.bed
