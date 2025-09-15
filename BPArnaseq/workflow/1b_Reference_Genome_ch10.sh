# NOTE: this code builds the reference genome of chromosome 10 in data/mm10_UCSC/STAR_index_chr10 using data/mm10_UCSC/genome_chr10.fa and data/mm10_UCSC/genes_chr10.gtf

# Set up the directory where the genome_chr10.fa and genes_chr10.gtf files are located:
cd data/mm10_UCSC

mkdir -p STAR_index_chr10

STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir STAR_index_chr10 \
     --genomeFastaFiles genome_chr10.fa \
     --sjdbGTFfile genes_chr10.gtf \
     --sjdbOverhang 99 \
     --genomeSAindexNbases 5 

