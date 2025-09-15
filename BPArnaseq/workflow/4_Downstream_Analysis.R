conda activate R
R
library(DESeq2)

# SET WORKING DIRECTORY
setwd("data/PreProcessing_Miniaturized")

countdata <- read.table("Counts/count_matrix.txt",
                        header = FALSE, row.names = 1)
metadata <- read.table("Counts/metadata.txt", header = TRUE, row.names = 1)
rownames(metadata) <- metadata$SampleName
colnames(countdata) <- metadata$SampleName

dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = metadata,
  design = ~ Condition
)

# Run DESeq analysis
dds <- DESeq(dds)

# Get results
res <- results(dds)
res <- as.data.frame(res)
res$contrast <- resultsNames(dds)[2]
res <- res[order(res$padj), ]

write.csv(res, file = "Counts/DESeq2_results.csv")