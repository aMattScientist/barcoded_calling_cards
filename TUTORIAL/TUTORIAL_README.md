This tutorial will walk through the commands of the barcoded_calling_cards pipeline, from the fastq input to the ccf output.
We provide all the intermediate files generated throughout the pipeline. 


# Example Analysis

First, follow the **Getting Started** section of the README 

The input for this tutorial is test.fq.gz
This file contains 100,000 reads randomly sampled from a MYOD1 calling cards dataset ( SRR17863637 ) 
For calling card analysis, we only use one read (Read 1) that contains a library barcode, part of the piggyBac transposon, the SRT barcode, the transposon-genome junction, then the genome.  

To run this tuturial, move the test.fq.gz file into your **raw** directory.  
To analyze the test data, change the header of bulkRNACallingCards.sh to read line 6. 

# Executing the main script 

From within your **'CALLINGCARDS'** directory, run the following command to submit the job.  

Using Slurm:   
```
sbatch CODE/bulkRNACallingCardsBarcodes.sh  
```

Using LSF:
```
bsub < CODE/bulkRNACallingCardsBarcodes.sh  
```

# Understanding the code 

The first command uses cutadapt to find reads that have the 'library barcode' specified in the sample manifest, and that also begin with the piggyBac LTR sequence. 

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
The output of this command, named NN_test_hg38.temp.fastq.gz, contains 95780 trimmed reads. 

The next code block uses umi_tools to find the SRT barcode, append the barcode to the readname, then trim the read up to and including the TTAA site. 
We use REGEX to specify the SRT barcode by position: (?P<cell_1>.{0,3})(?P<umi_1>.{4})(?P<discard_2>GTTAA) 
This searches each read for a 4 nt cell-barcode (cell_1), and 1 nt umi then the GTTAA from an inserted SRT

```
umi_tools extract \
    --stdin $OUT_STEM".temp.fastq.gz" \
    --extract-method=regex \
    --bc-pattern='(?P<cell_1>.{0,3})(?P<umi_1>.{4})(?P<discard_1>GGTTAA)' \
    --whitelist=$SAFELIST \
    --stdout $OUT_STEM".temp2.fastq.gz"  
```

The output of this command, named NN_test_hg38.temp2.fastq.gz has 93399 reads with the SRT barcoded appened to the read name (underlined below)  

>@SRR17863637.2609305_<ins>GAAG</ins> NB501801:394:HNC7JBGXG:4:13409:1839:11934 length=75  
>AAAGGCCCACCAGTTTAGTACACAGTGTGACTTTCAG  
>+  
>EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE  
>@SRR17863637.1098745_<ins>TGGT</ins> NB501801:394:HNC7JBGXG:2:13212:12252:18806 length=75  
>ACAGCAATAAGAAACTTCCTGTGACAACATAATAAAT  
>+  
>EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE  

The next block of code aligns the read to the reference genome using STAR. Code for NovoAlign is also available.

```
STAR --runThreadN 1 \
    --genomeDir $STAR_REF \
    --readFilesCommand zcat \
    --readFilesIn $OUT_STEM".temp2.fastq.gz" \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix $OUT_STEM"."
```

The mapping statistics are ouput in: NN_test_hg38.Log.final.out: 

                          Number of input reads |	93399
                      Average input read length |	36
                                    UNIQUE READS:
                   Uniquely mapped reads number |	82479
                        Uniquely mapped reads % |	88.31%

The mapped reads are output as NN_test_hg38.Aligned.sortedByCoord.out.bam

Next, we tag the bam file with the library barcode, INDEX1, and INDEX2 from the manifest. 

```
python $SCRIPTS/TagBam.py \
    --tag XP:Z:$BARCODE \
    $OUT_STEM".Aligned.sortedByCoord.out.bam" \
    $OUT_MAP_SORT_PREFIX"_tagged.bam"
```

After this, we tag the bam file wtih the SRT barcode. 

```
python $SCRIPTS/TagBamWithSrtBC.py \
	--tag XQ:Z \
	$OUT_MAP_SORT_PREFIX"_tagged3.bam" \
	$OUT_MAP_SORT_PREFIX"_tagged4.bam"
```

Next, we calculate the insertion sites of each read and add that info with the tag "XI" 
This command enforces that insertions are only valid if the genomic insertion site is "TTAA" for piggyBac-based experiments 

```
python $SCRIPTS/AnnotateInsertionSites.py \
    --transposase $transposase \
    -f \
    $OUT_MAP_SORT_PREFIX"_tagged4.bam" \
    $twobit_human \
    $OUT_MAP_SORT_PREFIX"_final.bam"
```

Finally, after indexing the bam file, we convert it to the qBED format (https://pubmed.ncbi.nlm.nih.gov/32941613/) 
Here we use the .ccf (calling card format) interchangeably with qBED. 

```
samtools index $OUT_MAP_SORT_PREFIX"_final.bam"
python $SCRIPTS/BamToCallingCard.py -b XP XJ XQ -i $OUT_MAP_SORT_PREFIX"_final.bam" -o $OUT_MAP_SORT_PREFIX"_unsorted.ccf" 
sort -k1V -k2n -k3n $OUT_MAP_SORT_PREFIX"_unsorted.ccf" > $OUT_MAP_SORT_PREFIX"_final.ccf"
```

The final output is a sorted qBED file named NN_test_hg38_map_sort_final.ccf containing 76554 lines.   
Each line corresponds to a genomic insertion of a barcoded SRT.

The resulting ccf file can be used as input for visualizing genomic insertions and for downstream peak calling using the mammalian calling cards toolkit: https://gitlab.com/rob.mitra/mammalian_cc_tools/
