#  README

## 1. Data sources

This task is based on publicly available sequencing data from a study of mouse female germline exposed to the environmental chemical Bisphenol A (BPA) or placebo. The dataset includes multiple samples under different conditions (e.g., BPA exposed vs. control) and was originally sequenced using the Illumina NovaSeq 6000.

Study: Effect of Brief Maternal Exposure to Bisphenol A on the Fetal Female Germline in a Mouse Model.
DOI: 10.1289/EHP15046
SRA repo: GSE261980

Four samples of the original study were subsampled (miniaturized) to contain only reads from chromosome 10 in a FASTQ format, stored in  `data/PreProcessing_Miniaturized/Fastq` and are used as the inputs for the workflow. 

## 2. How to download

The following samples were retrieved SRA and converted to fastq:  

| SRA Run      | Group        |
|--------------|--------------|
| SRR28395030 | PlaceboE13.5 |
| SRR28395031 | PlaceboE13.5 |
| SRR28395040 | BPAE13.5     |
| SRR28395041 | BPAE13.5     |

using 

```bash
prefetch ${SRR_FILES[@]}

for SRR in "${SRR_FILES[@]}"; do
    fasterq-dump $SRR \
        --split-files \
        --threads $THREADS \
        -O ../Fastq

    pigz -p $THREADS ../Fastq/${SRR}_*.fastq
done
```

Initial processing includes trimming, alignment and miniaturization, and count matrix generation. Tools required to implement the pipeline can be found in the metadata.yaml file. UCSC mm10 was used as reference genome, as in the original publication. 

Data analysis can be started directly using the miniaturized FASTQ files in  `data/PreProcessing_Miniaturized/Fastq`.


## 3. Pre-processing / subsampling

The four original selected SRA files...

| SRA Run      | Group        |
|--------------|--------------|
| SRR28395030 | PlaceboE13.5 |
| SRR28395031 | PlaceboE13.5 |
| SRR28395040 | BPAE13.5     |
| SRR28395041 | BPAE13.5     |

were miniaturized to contain only reads from chromosome 10 (ch10) using the aligned .bam files. The resulting subsample .bam files were transformed back to FASTQ files to replicate the original workflow. 

```bash
SRR_FILES=(SRR28395030 SRR28395031 SRR28395040 SRR28395041)

for SRR in "${SRR_FILES[@]}"; do
    samtools view -@ $THREADS -b Aligned/${SRR}_Aligned.sortedByCoord.out.bam 10 \
        > Miniaturized/${SRR}_chr10.bam
    samtools index Miniaturized/${SRR}_chr10.bam
    samtools fastq Miniaturized/${SRR}_chr10.bam \
        -1 ../PreProcessing_Miniaturized/Fastq/${SRR}_chr10_1.fastq.gz \
        -2 ../PreProcessing_Miniaturized/Fastq/${SRR}_chr10_2.fastq.gz \
        -@ $THREADS 
done
```

## 4. How the workflow works

The workflow files is stored in `workflow/` and it is divided into different steps:

### Step 1 – mm10 Reference Setup

**Purpose:** Generate the reference genome files 
**Tools:**   `STAR`
**Inputs:**  None
**Outputs:** reference genome, stored in `data/mm10`

The genome index was built using STAR using genome FASTA and GFT annotations. Genome index is stored in `data/mm10_UCSC` and workflow in `workflow/1_Reference_Genome.sh`

```bash
STAR --runThreadN ${THREADS} \
     --runMode genomeGenerate \
     --genomeDir "${STAR_INDEX_DIR}" \
     --genomeFastaFiles "${GENOME_FASTA}" \
     --sjdbGTFfile "${GTF_FILE}" \
     --sjdbOverhang $((READ_LENGTH - 1)) \
     --genomeSAindexNbases 14 \
```

The `workflow/1b_Reference_Genome_ch10.sh` allows to generate a reference genome only for chromosome 10 for practical purposes. 

### Step 2 – Original data QC, alignment and miniaturization.

**Purpose:** Generate miniaturized FASTQ files
**Tools:**   `prefetch` `fasterq-dump` `pigz` `fastp` `fastqc` `STAR` `samtools` 
**Inputs:**  SRA list files and mm10 reference genome from `data/mm10`
**Outputs:** subsampled FASTQ, stored in `data/PreProcessing_Miniaturized`

Following the download of SRA files, they were FASTQ converted, trimmed, QC-checked, aligned and miniaturized. The resulting putput files are stored in `data/PreProcessing` and the corresponding workflow is in `workflow/2_DataRetrieval_and_Miniaturization.sh`. 

Output miniaturized FASTQ files used for downstream processing are in `data/PreProcessing_Miniaturized`

### Step 3 - Miniaturized data workflow

**Purpose:** Subsample workflow for QC, alignment and feature counts
**Tools:**  `fastqc` `STAR` `samtools` `featureCounts`
**Inputs:**  subsampled FASTQ files from `data/PreProcessing_Miniaturized/Fastq` and mm10 reference genome from `data/mm10`
**Outputs:** aligned .bam files, .txt count matrix and metadata 

FASTQ miniaturized files were QC-checked and aligned. Counts were obtained using the following code. Given the miniaturized nature of the files, only features having at least 1 count in any sample were retained in the count_matrix.

```bash
featureCounts \
    -T 32 \
    -p -B -C \
    -a $GenRefDir/genes_chr10.gtf \
    -o $CountsDir/gene_counts_chr10.txt \
    --extraAttributes gene_name \
    $AlignmentDir/*.bam

cut -f 1,8- Counts/gene_counts_chr10.txt | tail -n +3 \
| awk '{ sum=0; for(i=2;i<=NF;i++) sum+=$i; if(sum>0) print $0 }' > Counts/count_matrix.txt
```

### Step 4 – Downstream analysis

**Purpose:** Differential expression analysis 
**Tools:** `DESeq2`
**Inputs:**  count matrix and metadata txt files from `data/PreProcessing_Miniaturized/Counts`
**Outputs:** DESeq2_results.csv file with differentially expressed features between BPA and Placebo

The resulting count_matrix from the miniaturized FASTQ files was imported in R and converted into a Deseq2 object using... 

```bash
dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = metadata,
  design = ~ Condition
)
```
This step was followed by DESeq2 analysis using Wald test between BPA and Placebo samples. The resulting output was saved as a .csv file in data/`data/PreProcessing_Miniaturized/Counts`
