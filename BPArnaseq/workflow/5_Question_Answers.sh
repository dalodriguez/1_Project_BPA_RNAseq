# ANSWER TO QUESTION 1
cd ./data/PreProcessing_Miniaturized/
cd Fastq && wc -l *.fastq | sort -nr | awk 'NR==2 {print $2}' | sed 's/_chr10_.*$//'

# ANSWER TO QUESTION 2
gunzip -c SRR28395030_chr10_1.fastq | head -2 | tail -1 | wc -c

# ANSWER TO QUESTION 3
awk '$2 > 200'  Counts/count_matrix.txt | wc -l

#Q4 ANSWER TO QUESTION 4
awk -F, '$1=="\"Erbb3\"" {print $2}' Counts/DESeq2_results.csv | awk '{sum+=$1; count++} END{print sum/count}'

# Q5 ANSWER TO QUESTION 5
tail -n +2 Counts/DESeq2_results.csv | sort -t, -k6,6n | head -n 3 | cut -d, -f1 | tr -d '"'