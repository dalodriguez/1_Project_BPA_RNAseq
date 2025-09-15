
# NOTE: do not run this script, the output files are in data/mm10_UCSC

OUTPUT_DIR="./data/mm10_UCSC"
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

# Download prebuilt iGenome archive
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
tar -xvzf Mus_musculus_UCSC_mm10.tar.gz
rm Mus_musculus_UCSC_mm10.tar.gz
