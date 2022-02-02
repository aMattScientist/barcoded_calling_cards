#!/bin/bash
#SBATCH -n 1  
#SBATCH -N 1 
#SBATCH --mem=30000
#SBATCH -o call_peaks_NGN.out#standard out goes here
#SBATCH -e call_peaks_NGN.err # standard error goes here
#SBATCH -J PEAK

module load ccf_tools

python $CCF_TOOLS/call_peaks_macs.py -e NGN1_4reps_sorted.ccf -o NGN1_MACS.sig \
	-t /calling_card_ref/human/TTAA_hg38_ccf.txt -a /calling_card_ref/human/refGene.hg38.Sorted.bed \
	-b HEK293T_hypPB_sorted.ccf -pc 0.001 --peak_finder_pvalue 0.01 --window 1000 --step 500 --pseudocounts 0.2
