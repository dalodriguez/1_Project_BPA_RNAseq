#!/bin/bash

# NOTE: this code replicates the paper's pipeline using miniaturized (chromosome 10) Fastq files

# SET UP DIRECTORIES PROPERLY
cd "data"
GenRefDir="mm10_UCSC"
FastqDir="PreProcessing_Miniaturized/Fastq"
AlignmentDir="PreProcessing_Miniaturized/Aligned"
CountsDir="PreProcessing_Miniaturized/Counts"

SRR_FILES=(SRR28395030 SRR28395031 SRR28395040 SRR28395041)
GROUP=(PlaceboE13.5 PlaceboE13.5  BPAE13.5 BPAE13.5)
ID=(PlaceboE13.5.1 PlaceboE13.5.2 BPAE13.5.1 BPAE13.5.2)

STAR --genomeDir $GenRefDir/STAR_index_chr10 \
     --genomeLoad LoadAndExit

for SRR in "${SRR_FILES[@]}"; do
    STAR \
        --runThreadN 16 \
        --genomeDir $GenRefDir/STAR_index_chr10 \
        --readFilesIn $FastqDir/${SRR}_chr10_1.fastq $FastqDir/${SRR}_chr10_2.fastq \
        --readFilesCommand zcat \
        --outFileNamePrefix $AlignmentDir/${SRR}_ \
        --outSAMtype BAM SortedByCoordinate

    samtools index "$AlignmentDir/${SRR}_Aligned.sortedByCoord.out.bam"
done

STAR --genomeDir $GenRefDir/STAR_index_chr10 \
     --genomeLoad Remove

featureCounts \
    -T 32 \
    -p -B -C \
    -a $GenRefDir/genes_chr10.gtf \
    -o $CountsDir/gene_counts_chr10.txt \
    --extraAttributes gene_name \
    $AlignmentDir/*.bam

cut -f 1,8- $CountsDir/gene_counts_chr10.txt | tail -n +3 \
| awk '{ sum=0; for(i=2;i<=NF;i++) sum+=$i; if(sum>0) print $0 }' > $CountsDir/count_matrix.txt

echo -e "SampleID\tCondition\tSampleName" > $CountsDir/metadata.txt
for i in "${!SRR_FILES[@]}"; do
    echo -e "${SRR_FILES[$i]}_Aligned.sortedByCoord.out.bam\t${GROUP[$i]}\t${ID[$i]}" >> $CountsDir/metadata.txt
done

