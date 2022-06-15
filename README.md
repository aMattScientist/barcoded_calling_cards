This is a fork of Arnav Moudgil's Calling Card analysis repository (https://github.com/arnavm/calling_cards), with minor modifications to enable analysis of barcoded self-reporting transposon calling cards. Reference: https://www.biorxiv.org/content/10.1101/2021.04.15.439516

## Dependencies

The Python scripts in this folder require `python3`, preferrably a relatively recent version (e.g. ≥ 3.5). In addition, you will need to install the following modules:
- `pysam`
- `numpy`
- `pandas`
- `twobitreader`
- `pybedtools`
- `scipy`

The easiest way to install them is:
- If `python3` is the default:
`pip install pysam numpy pandas twobitreader pybedtools scipy`
- Otherwise:
`pip3 install pysam numpy pandas twobitreader pybedtools scipy`

To figure out which version of `python` is installed by default:
- `python -V`

The main shell script was designed to be executed from a Linux/Unix based high-performance computing environment using either Slurm or LSF job managers.   The software requires several common bioinformatics modules including:   
- `cutadapt`
- `umi-tools`
- `star or (NovoAlign)` 
- `samtools`
- `bedtools`
- `java`

# Getting Started

Within your Linux/Unix based high-performance computing environment, create a new directory. For example, name it CALLINGCARDS.  
Within this directory, make 3 more directories:  output_and_analysis, raw, and CODE.  
Download the contents of this repository to the CODE directory.  

If your cluster uses Slurm, append the SLURM_header to bulkRNACallingCards and save this as a shell script.  
```
cat SLURM_header.txt bulkRNACallingCardsBarcodes > bulkRNACallingCardsBarcodes.sh
```

Alternatively, if your cluster uses LSF, append the LSF_header.  
```
cat LSF_header.txt bulkRNACallingCardsBarcodes > bulkRNACallingCardsBarcodes.sh  
```

Change the relevant path names in bulkRNACallingCardsBarcodes.sh as appropriate    

Move the barcode_safelist.txt and manifest.csv out of the CODE directory into the **'CALLINGCARDS'** directory.  

# Example Analysis

To test the software, download one of our barcoded calling card datasets from the Sequence Read Archive (SRA).   

Navigate into the **'raw'** directory and download the MYOD1 calling card dataset that was jointly prepared with BRB-seq data.   
The SRA number for this dataset is : SRR17863637    
The link for this dataset is:  https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR17863637   
The GEO sample information is: GSM5857636: MYOD1 joint BRB.   

One way to download this data is using the SRA toolkit:  
```
ml sratoolkit   
fastq-dump --gzip SRR17863637  
#Rename the file as MYOD1_jointBRB.fastq.gz  
mv SRR17863637.fastq.gz > MYOD1_jointBRB.fastq.gz  
```

Alternatively, we have provided a sub-sampled dataset in the **TUTORIAL** folder, named: test.fq.gz   
Move this file to your **raw** directory

The SLURM_header and LSF_header are configured to direct bulkRNACallingCards.sh to analyze line 4 of the manifest (MYOD1_jointBRB.fastq.gz) 
To analyze the test data instead, change the headers to read line 6. 

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

# Expected Results

This script should generate a file named NN_MYOD1_jointBRB_hg38_map_sort_final.ccf with approximately 587277 lines. Each line corresponds to a genomic insertion of a barcoded SRT.  

A full test example including all expected intermediate files is provided in the **TUTORIAL** folder.     
The resulting ccf file can be used as input for downstream peak calling using the mammalian calling cards toolkit:   https://gitlab.com/rob.mitra/mammalian_cc_tools/
