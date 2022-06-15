This tutorial will walk through the commands of the barcoded_calling_cards pipeline, from the fastq input to the ccf output.
We provide all the intermediate files generated throughout the pipeline. 


# Example Analysis

First, follow the **Getting Started** section of the README 

Input: test.fq.gz
This file contains 100,000 reads randomly sampled from a MYOD1 calling cards dataset ( SRR17863637 ) 
For calling card analysis, we only use one read (Read 1) that contains a library barcode, part of the piggyBac transposon, the SRT barcode, the transposon-genome junction, then the genome.   

To run this tuturial, move the test.fq.gz file into your **raw** directory.  
To analyze the test data, change the header of bulkRNACallingCards.sh to read line 6. 

# Executing the main script 

From within your **'CALLINGCARDS'** directory, run the following command to submit the job.  

Using Slurm:   
>sbatch CODE/bulkRNACallingCardsBarcodes.sh  

Using LSF:
>bsub < CODE/bulkRNACallingCardsBarcodes.sh  

# Understanding the code 

The first command uses cutadapt to find reads that have the 'library barcode' specified in the sample manifest, and that also begin with the piggyBAC LTR sequence. 

```
cutadapt \
    -g "^"$BARCODE$Long_PB_LTR \
    -o $OUT_STEM".fastq.gz" \
    --minimum-length 1 \
    --discard-untrimmed \
    -e 0.2 \
    $PROJECT_RAW/$FILENAME
```

The output of this command is named: NN_test_hg38.fastq.gz and contains 98603 reads matching this pattern. 
The next block of code also uses cutadapt and trims any remaining Nextera adapters, and any reads less than 25 bases 

```
cutadapt \
    -a "CTGTCTCTTATACACATCTCCGAGCCCACGAGACTNNNNNNNNNNTCTCGTATGCCGTCTTCTGCTTG" \
    -o $OUT_STEM".temp.fastq.gz" \
    --minimum-length 25 \
    $OUT_STEM".fastq.gz"
```
The output of this command, named NN_test_hg38.temp.fastq.gz, contains 95780 reads. 


