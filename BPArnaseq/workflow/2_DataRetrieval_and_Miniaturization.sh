#!/bin/bash

# NOTE: this code generates the miniaturized Fastq files, that will be placed in ./BPArnaseq/data/PreProcessing_Miniaturized/Fastq

mkdir -p ./data/PreProcessing/{Raw,Fastq,Aligned,Counts, Trimmed}
mkdir -p ./data/PreProcessing_Miniaturized/{Fastq,Aligned,Counts, QC}
mkdir -p ./data/mm10_UCSC

SRR_FILES=(SRR28395030 SRR28395031 SRR28395040 SRR28395041)
GROUP=(PlaceboE13.5 PlaceboE13.5  BPAE13.5 BPAE13.5)
ID=(PlaceboE13.5.1 PlaceboE13.5.2 BPAE13.5.1 BPAE13.5.2)
THREADS=32

cd ./data/PreProcessing/Raw

prefetch ${SRR_FILES[@]}

for SRR in "${SRR_FILES[@]}"; do
    fasterq-dump $SRR \
        --split-files \
        --threads $THREADS \
        -O ../Fastq

    pigz -p $THREADS ../Fastq/${SRR}_*.fastq
done


cd ./data/PreProcessing/


for SRR in "${SRR_FILES[@]}"; do
    echo "Subsampling 1 Fastq files for ${SRR}..."
    seqtk sample -s100 Fastq/${SRR}_1.fastq.gz 0.005 > Fastq/${SRR}_1sub.fastq
    echo "Subsampling 2 Fastq files for ${SRR}..."
    seqtk sample -s100 Fastq/${SRR}_2.fastq.gz 0.005 > Fastq/${SRR}_2sub.fastq
done


for SRR in "${SRR_FILES[@]}"; do
    fastp \
        --in1 Fastq/${SRR}_1.fastq.gz \
        --in2 Fastq/${SRR}_2.fastq.gz \
        --out1 Trimmed/${SRR}_1.trimmed.fastq.gz \
        --out2 Trimmed/${SRR}_2.trimmed.fastq.gz \
        --trim_tail1 5 \
        --trim_tail2 5 \
        --detect_adapter_for_pe \
        --length_required 50 \
        --thread $THREADS \
        --html Trimmed/${SRR}_fastp.html \
        --json Trimmed/${SRR}_fastp.json
done

STAR --genomeDir ../mm10_UCSC/STAR_index \
     --genomeLoad LoadAndExit

for SRR in "${SRR_FILES[@]}"; do
    STAR \
        --runThreadN $THREADS \
        --genomeDir ../mm10_UCSC/STAR_index \
        --genomeLoad LoadAndKeep \
        --readFilesIn Trimmed/${SRR}_1.trimmed.fastq.gz Trimmed/${SRR}_2.trimmed.fastq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix Aligned/${SRR}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterMismatchNmax 3 \
        --outFilterMultimapNmax 1 \
        --alignSJoverhangMin 8
done

STAR --genomeDir ../mm10_UCSC/STAR_index \
     --genomeLoad Remove

for SRR in "${SRR_FILES[@]}"; do
    samtools view -@ $THREADS -b Aligned/${SRR}_Aligned.sortedByCoord.out.bam 10 \
        > Miniaturized/${SRR}_chr10.bam
    samtools index Miniaturized/${SRR}_chr10.bam
    samtools fastq Miniaturized/${SRR}_chr10.bam \
        -1 ../PreProcessing_Miniaturized/Fastq/${SRR}_chr10_1.fastq.gz \
        -2 ../PreProcessing_Miniaturized/Fastq/${SRR}_chr10_2.fastq.gz \
        -@ $THREADS 
done