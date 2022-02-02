#!/bin/bash
#
#SBATCH --job-name=clean_code
#SBATCH --output=logs/%x_%a.out
#SBATCH --error=logs/%x_%a.err  
#SBATCH --array=2-2
#SBATCH --mem=32000
#SBATCH --cpus-per-task=8

# Initialize settings
PROJECT_DIR=  ## your directory
PROJECT_OUT=$PROJECT_DIR"/output_and_analysis"
PROJECT_RAW=$PROJECT_DIR"/raw"
SAMPLE_LINE=$( sed -n ${SLURM_ARRAY_TASK_ID}p manifest.csv )
SAFELIST=$'barcode_safelist.txt'

# Initialize important global variables
FILENAME=
BARCODE=
INDEX1=
INDEX2=

# Read local variables from $SAMPLE_LINE in subshell,
# then pass them along to the global versions
while IFS=, read -r filename barcode index1 index2
do
    { FILENAME=$filename; }
    { BARCODE=$barcode; }
    { INDEX1=$index1; }
    { INDEX2=$index2; }
done <<< "$SAMPLE_LINE"

# and a corresponding file of barcodes.
BASE=`basename $FILENAME`
STEM=$BARCODE"_""${BASE%%.*}"

# Load the necessary modules
module purge
module load novoalign
module load samtools
module load bedtools
module load java

# Prepare output variables
OUT_STEM=$PROJECT_OUT/$STEM
OUT_SAM=$OUT_STEM".sam"
OUT_BAM=$OUT_STEM".bam"
OUT_MAP_SORT_PREFIX=$OUT_STEM"_map_sort"
OUT_BAM_MAP=$OUT_MAP_SORT_PREFIX".bam"
OUT_BC=$OUT_STEM".bc.fastq"
OUT_BED=$OUT_MAP_SORT_PREFIX".bed"
OUT_BEDGRAPH=$OUT_MAP_SORT_PREFIX".bedgraph"

SB_ITR="TAAGTGTATGTAAACTTCCGACTTCAACTGTA"
alt_SB_ITR="AAGTGTATGTAAACTTCCGACTTCAACTGTA"

#PB_LTR="TTTACGCAGACTATCTTTNNNNGGTTAA"  ## changed CTAG to NNNN to allow any basepair here for mutagenesis screen
PB_LTR="TTTACGCAGACTATCTTT"  
Long_PB_LTR="CGTCAAT"$PB_LTR

transposase=PB
genome=hg38


module load cutadapt

# Trim and demultiplex by LTR + GGTTAA insertion site (exact matching)
# use PB_LTR, not Long_PB_LTR for older data (just check the first bases after your barcode in Read1 to know)

cutadapt \
    -g "^"$BARCODE$Long_PB_LTR \
    -o $OUT_STEM".fastq.gz" \
    --minimum-length 1 \
    --discard-untrimmed \
    -e 0.2 \
    --no-indels \
    $PROJECT_RAW/$FILENAME

# Trim any trailing Nextera adapters (allowing mismatches)
cutadapt \
    -a "CTGTCTCTTATACACATCTCCGAGCCCACGAGACTNNNNNNNNNNTCTCGTATGCCGTCTTCTGCTTG" \
    -o $OUT_STEM".temp.fastq.gz" \
    --minimum-length 25 \
    $OUT_STEM".fastq.gz"

# Add SRT barcode to read name with umi-tools, also cuts off this adaptor 
# Use --filter-cell-barcode (after bc-pattern line) to remove reads that are not safelisted if you want (error-detecting bcs)
# Otherwise, it will keep any identified barcode. Barcodes are automatically tagged in the output fastq file
# Barcode pattern C = cell barcode position, and N = umi (required). Our barcode is not really of this form, instead, N should always be G
# This sequence matches the piggyBac TR that we mutagenized  

module unload cutadapt
module load umi-tools

umi_tools extract \
	  --stdin $OUT_STEM".temp.fastq.gz" \
	  --bc-pattern=CCCCNGTTAA \
	  --whitelist=$SAFELIST \
	  --stdout $OUT_STEM".fastq.gz"

module unload umi-tools

# Align the trimmed and SRT_BC annoated reads
novoalign \
    -d /novoalign_indexes/$genome/$genome.nvx \
    -f $OUT_STEM".fastq.gz" \
    -n 40 \
    -o SAM \
    -o SoftClip > $OUT_SAM

# Filter only mapped reads, convert to BAM, and sort
samtools view \
    -bS -h -F 260 \
    $OUT_SAM | \
samtools sort - -o $OUT_BAM_MAP

# Tag reads with barcodes

python3 TagBam.py \
    --tag XP:Z:$BARCODE \
    $OUT_BAM_MAP \
    $OUT_MAP_SORT_PREFIX"_tagged.bam"

# Tag reads wih i7 indexes
python3 TagBam.py \
    --tag XJ:Z:$INDEX1 \
    $OUT_MAP_SORT_PREFIX"_tagged.bam" \
    $OUT_MAP_SORT_PREFIX"_tagged2.bam"

# Tag reads with i5 indexes
python3 TagBam.py \
    --tag XK:Z:$INDEX2 \
    $OUT_MAP_SORT_PREFIX"_tagged2.bam" \
    $OUT_MAP_SORT_PREFIX"_tagged3.bam"


# Tag reads with SRT barcode
python3 TagBamWithSrtBC.py \
	--tag XQ:Z \
	$OUT_MAP_SORT_PREFIX"_tagged3.bam" \
	$OUT_MAP_SORT_PREFIX"_tagged4.bam"

# Tag reads with transposon insertions
python3 AnnotateInsertionSites.py \
    --transposase $transposase \
    -f \
    $OUT_MAP_SORT_PREFIX"_tagged4.bam" \
    /novoalign_indexes/$genome/$genome.2bit \
    $OUT_MAP_SORT_PREFIX"_final.bam"

# Index the final BAM file
samtools index $OUT_MAP_SORT_PREFIX"_final.bam"

# Get an unsorted list of unique insertions
python3 BamToCallingCard.py -b XP XJ XQ -i $OUT_MAP_SORT_PREFIX"_final.bam" -o $OUT_MAP_SORT_PREFIX"_unsorted.ccf" 

# Sort the CCF file 
sort -k1V -k2n -k3n $OUT_MAP_SORT_PREFIX"_unsorted.ccf" > $OUT_MAP_SORT_PREFIX"_final.ccf"

# Clean up
rm $OUT_SAM
rm $OUT_BAM
rm $OUT_STEM".fastq.gz"
rm $OUT_STEM".temp.fastq.gz"
rm $OUT_BAM_MAP
rm $OUT_MAP_SORT_PREFIX"_tagged.bam"
rm $OUT_MAP_SORT_PREFIX"_tagged2.bam"
rm $OUT_MAP_SORT_PREFIX"_tagged3.bam"
rm $OUT_MAP_SORT_PREFIX"_tagged4.bam"
rm $OUT_MAP_SORT_PREFIX"_unsorted.ccf"
