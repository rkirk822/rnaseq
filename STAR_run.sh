#!/bin/bash

# STAR_run
# Map a list of fastq files using STAR.
#
# Need to make this build a command, instead of defining all possible parameters
# as variables.  In bash, do this by defining an array rather than a string:
# args=(-s "$subject" --flag "arg with spaces")
# mail "${args[@]}"
# 
# TIME: 
# Roughly 10m per 5G data file.
# 5h 30m for all 28 RORb data files.
# 4h 45m for the 15 EMX-Cre data files.
# 3hr 30m for the 14 PV data files.
# ~30m for the Nuc-seq data (12 files?), minus one file.

###################################################
# Stuff to set
pathSTAR='/opt/STAR/bin/MacOSX_x86_64/STAR'
pathTrimmed='/Volumes/CodingClub1/RNAseq/RORb/trimmed/trimmed_m20_q20/'
trimmedSuffix='_trimmed_m20_q20.fastq' # used to construct input file names
pathMapped='/Volumes/CodingClub1/RNAseq/RORb/mapped/mapped_sams/'
mappedSuffix='_mapped_' # last part of prefix to output files; include underscore at the end
cores=8 # runThreadN; number of available cores
# THIS WAS GENERATED TO BE OPTIMAL FOR FILES WITH MAX READ LENGTH 75 (or 76?):
genomeDir='/Volumes/CodingClub1/STAR_stuff/indexes/generated_with_gtf_from_mysql/'
# outSAMtype='BAM SortedByCoordinate'

# sampleNames='cb_p036_1 cb_p036_2 cb_p036_3 cb_p036_4 L4_RORb_1 L4_RORb_2 L4_RORb_3 L4_RORb_4'
# sampleNames='cb_p033_4 cb_p033_5 cb_p033_6 cb_p033_7 L4_RORht_1 L4_RORht_2 L4_RORht_3 L4_RORht_4'
# sampleNames='AACTCG_S2_R1_001 ATGACA_S7_R1_001 CTCCTA_S3_R1_001 GTCGAG_S5_R1_001 ACAACA_S6_R1_001 CCAACA_S10_R1_001 CTCGAA_S12_R1_001 TCTAGG_S11_R1_001 ATCGTC_S9_R1_001 CGTCTA_S1_R1_001 GCACAC_S8_R1_001 TCTTCG_S4_R1_001'
# sampleNames='GTCGAG_S5_R1_001'
# sampleNames='PVCtlEarly_1 PVCtlEarly_2 PVCtlEarly_3 PVCtlLate_1 PVCtlLate_2 PVCtlLate_3 PVCtlLate_4 PVTTXEarly_1 PVTTXEarly_2 PVTTXEarly_3 PVTTXLate_1 PVTTXLate_2 PVTTXLate_3 PVTTXLate_4'
# sampleNames='EMXCtlEarly_1 EMXCtlEarly_2 EMXCtlEarly_3 EMXCtlEarly_4 EMXCtlLate_1 EMXCtlLate_2 EMXCtlLate_3 EMXCtlLate_4 EMXTTXEarly_1 EMXTTXEarly_2 EMXTTXEarly_3 EMXTTXLate_1 EMXTTXLate_2 EMXTTXLate_3 EMXTTXLate_4'
# sampleNames=EMXCtlEarly_1
sampleNames='BF_RORbHTp2_1 BF_RORbHTp2_2 BF_RORbHTp2_3 BF_RORbHTp2_4 BF_RORbKOp2_1 BF_RORbKOp2_2 BF_RORbKOp2_3 BF_RORbKOp2_4 BF_RORbHTp7_1 BF_RORbHTp7_2 BF_RORbHTp7_3 BF_RORbHTp7_4 BF_RORbKOp7_1 BF_RORbKOp7_2 BF_RORbKOp7_3 BF_RORbKOp7_4 BF_RORbHTp30_1 BF_RORbHTp30_2 BF_RORbHTp30_3 BF_RORbHTp30_4 BF_RORbKOp30_1 BF_RORbKOp30_2 BF_RORbKOp30_3 BF_RORbKOp30_4 BF_RORbHTp200_1 BF_RORbHTp200_2 BF_RORbHTp200_3 BF_RORbHTp200_4'
# # ENCODE SETTINGS - long RNA
outFilterType='BySJout'
outFilterMultimapNmax=20
alignSJoverhangMin=8  # for paired reads?
alignSJDBoverhangMin=1  # for paired reads?
outFilterMismatchNmax=999
outFilterMismatchNoverLmax=0.04
alignIntronMin=20   # for paired reads?
alignIntronMax=1000000  # for paired reads?
alignMatesGapMax=1000000 # for paired reads?
# # ENCODE SETTINGS - short RNA
# outFilterMismatchNoverLmax=0.03
# outFilterMultimapNmax=20  #same as long RNA
# alignIntronMax=1
# outFilterMatchNmin=16
# outFilterScoreMinOverLread=0
# outFilterMatchNminOverLread=0
###################################################



for sn in $sampleNames;
do

    echo Mapping $sn with STAR

    # Check if you already did this one
    logfile=$pathMapped$sn$mappedSuffix'Log.out'
    if [ -f $logfile ]; then echo 'File ' $logFile ' already exists. Skipping sample.'; continue; fi

# # Long RNA
#     $pathSTAR --runThreadN $cores --genomeDir $genomeDir --readFilesIn $pathTrimmed$sn$trimmedSuffix --outFilterType $outFilterType --outFilterMultimapNmax $outFilterMultimapNmax --alignSJoverhangMin $alignSJoverhangMin --alignSJDBoverhangMin $alignSJDBoverhangMin --outFilterMismatchNmax $outFilterMismatchNmax --outFilterMismatchNoverLmax $outFilterMismatchNoverLmax --alignIntronMin $alignIntronMin --alignIntronMax $alignIntronMax --alignMatesGapMax $alignMatesGapMax --outSAMtype $outSAMtype --outFileNamePrefix $pathMapped$sn$mappedSuffix

# # Long RNA, no outSAMtype because I can't figure out how to explicitly set it to default which is SAM
    $pathSTAR --runThreadN $cores --genomeDir $genomeDir --readFilesIn $pathTrimmed$sn$trimmedSuffix --outFilterType $outFilterType --outFilterMultimapNmax $outFilterMultimapNmax --alignSJoverhangMin $alignSJoverhangMin --alignSJDBoverhangMin $alignSJDBoverhangMin --outFilterMismatchNmax $outFilterMismatchNmax --outFilterMismatchNoverLmax $outFilterMismatchNoverLmax --alignIntronMin $alignIntronMin --alignIntronMax $alignIntronMax --alignMatesGapMax $alignMatesGapMax --outFileNamePrefix $pathMapped$sn$mappedSuffix

# # Short RNA
#     $pathSTAR --runThreadN $cores --genomeDir $genomeDir --readFilesIn $pathTrimmed$sn$trimmedSuffix --outFilterMismatchNoverLmax $outFilterMismatchNoverLmax --outFilterMultimapNmax $outFilterMultimapNmax --alignIntronMax $alignIntronMax --outFilterMatchNmin $outFilterMatchNmin --outFilterScoreMinOverLread $outFilterScoreMinOverLread --outFilterMatchNminOverLread $outFilterMatchNminOverLread --outFileNamePrefix $pathMapped$sn$mappedSuffix

done



