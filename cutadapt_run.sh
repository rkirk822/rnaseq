#!/bin/bash

# cutadapt_run.sh
# Run through a list of fastq files (or fasta?) and trim them with cutadapt.
# Save results as originalFileName_trimmed.originalExtension.
# Cutadapt's report, which by default is displayed in the terminal, goes to originalFileName_report.txt.
# That report includes the command line parameters, so whatever you use in this script as the adapter
# sequences, quality cutoff, etc, will appear there.
#
# TIME:
# For 5-6 GB fastq files, about 15m.
# For 9-10 GB, closer to 30m.
# For all the RORb files, total ~7.5hrs.

# Parameters to consider:
# 	quality cutoff; Yasu used 20 and I'm seeing online examples using both 20 and 30
#   error tolerance
#	Allow for possibility of >1 adapters on a given read (Anton allowed for up to 3)
# 	Trim flanking Ns (--trim-n); based on what Anton said, no reason not to


###################################################
# # Stuff to set
pathRaw='/Volumes/CodingClub1/RNAseq/NucSeq/raw/'
pathTrimmed='/Volumes/CodingClub1/RNAseq/NucSeq/trimmed/trimmed_m20_q20/'
# pathRaw='/Volumes/CodingClub1/RNAseq/TTX/raw/'
# pathTrimmed='/Volumes/CodingClub1/RNAseq/TTX/trimmed/trimmed_m20_q20/'
extension=.fastq.gz
adapterName=TruSeq_Index
adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
adapterName2=TruSeq_Universal
adapter2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
qualityCutoff=20
minLen=20 # throw away reads shorter than this (set to 0 to not filter this way)
nAdapt=3 # cutadapt will assume there could be as many adapters as this on a given read
fileSuffix='_m20_q20' # to distinguish files trimmed with different parameter values; else make it ''
# fileSuffix='_q20'
cores=8 # number of available cores
###################################################

sampleNames='AACTCG_S2_R1_001 ATGACA_S7_R1_001 CTCCTA_S3_R1_001 GTCGAG_S5_R1_001 ACAACA_S6_R1_001 CCAACA_S10_R1_001 CTCGAA_S12_R1_001 TCTAGG_S11_R1_001 ATCGTC_S9_R1_001 CGTCTA_S1_R1_001 GCACAC_S8_R1_001 TCTTCG_S4_R1_001'
# sampleNames='PVCtlEarly_1 PVCtlEarly_2 PVCtlEarly_3 PVCtlLate_1 PVCtlLate_2 PVCtlLate_3 PVCtlLate_4 PVTTXEarly_1 PVTTXEarly_2 PVTTXEarly_3 PVTTXLate_1 PVTTXLate_2 PVTTXLate_3 PVTTXLate_4'
# sampleNames='EMXCtlEarly_1 EMXCtlEarly_2 EMXCtlEarly_3 EMXCtlEarly_4 EMXCtlLate_1 EMXCtlLate_2 EMXCtlLate_3 EMXCtlLate_4 EMXTTXEarly_1 EMXTTXEarly_2 EMXTTXEarly_3 EMXTTXLate_1 EMXTTXLate_2 EMXTTXLate_3 EMXTTXLate_4'
# sampleNames='BF_RORbHTp200_1 BF_RORbHTp200_2 BF_RORbHTp200_3 BF_RORbHTp200_4 BF_RORbHTp2_1 BF_RORbHTp2_2 BF_RORbHTp2_3 BF_RORbHTp2_4 BF_RORbHTp30_1 BF_RORbHTp30_2 BF_RORbHTp30_3 BF_RORbHTp30_4 BF_RORbHTp7_1 BF_RORbHTp7_2 BF_RORbHTp7_3 BF_RORbHTp7_4 BF_RORbKOp2_1 BF_RORbKOp2_2 BF_RORbKOp2_3 BF_RORbKOp2_4 BF_RORbKOp30_1 BF_RORbKOp30_2 BF_RORbKOp30_3 BF_RORbKOp30_4 BF_RORbKOp7_1 BF_RORbKOp7_2 BF_RORbKOp7_3 BF_RORbKOp7_4'

echo $sampleNames

# Check that you have access to cutadapt
# (see the user guide re adding this to PATH)
cutadapt 2>/dev/null
if [ "$?" -ne "0" ]; then PATH=$PATH:$HOME/.local/bin; fi

for sn in $sampleNames;
do

    # Define names of files to write so we can check if they exist
    outfile=$pathTrimmed$sn'_trimmed'$fileSuffix$extension
    reportfile=$pathTrimmed$sn$fileSuffix'_cutadapt_report.txt'
    echo Trimmed reads will be saved to $outfile
    if [ -f $outfile ]; then echo 'Trimmed file already exists. Skipping file.'; continue; fi
    if [ -f $reportfile ]; then echo 'Report file already exists.  Skipping file.'; continue; fi

    infile=$pathRaw$sn$extension

    cutadapt -a $adapterName=$adapter -a $adapterName2=$adapter2 $infile -o $outfile --trim-n -q $qualityCutoff -m $minLen -n $nAdapt -j $cores > $reportfile

done



