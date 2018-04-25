#!/bin/bash

# STAR_run_paired.sh
#
# Version of STAR_run.sh for paired-end data in separate files.
#
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
pathTrimmed='/Volumes/CodingClub1/RNAseq/smallRNA/20180402_NXT/barcodes_new/Emma/trimmed/'
trimmedSuffix='_trimmed.fastq' # used to construct input file names
pathMapped='/Volumes/CodingClub1/RNAseq/smallRNA/20180402_NXT/barcodes_new/Emma/mapped/'
mappedSuffix='_' # last part of prefix to output files; include underscore at the end
cores=8 # runThreadN; number of available cores
# THIS WAS GENERATED TO BE OPTIMAL FOR FILES WITH MAX READ LENGTH 75 (or 76?):
# printf 'Using indices best for data with max length 50'
# genomeDir='/Volumes/CodingClub1/STAR_stuff/indexes/generated_with_gtf_from_mysql_maxlen50/'
genomeDir='/Volumes/CodingClub1/STAR_stuff/indexes/generated_with_gtf_from_mysql/'
# outSAMtype='BAM SortedByCoordinate'

# sampleNames='GTAGAG_S1 GTCCGC_S2'
sampleNames='GTAGAG_S1 GTCCGC_S2 GTTTCG_S5 CGTACG_S6 GAGTGG_S7 GGTAGC_S8'


# # ENCODE SETTINGS - long RNA
# outFilterType='BySJout'
# outFilterMultimapNmax=20
# alignSJoverhangMin=8  # for paired reads?
# alignSJDBoverhangMin=1  # for paired reads?
# outFilterMismatchNmax=999
# outFilterMismatchNoverLmax=0.04
# alignIntronMin=20   # for paired reads?
# alignIntronMax=1000000  # for paired reads?
# alignMatesGapMax=1000000 # for paired reads?

# # ENCODE SETTINGS - short RNA
outFilterMismatchNoverLmax=0.03
outFilterMultimapNmax=20  #same as long RNA
alignIntronMax=1
# outFilterMatchNmin=16
outFilterScoreMinOverLread=0
outFilterMatchNminOverLread=0

# Number of matched bases you have to have to count as an alignment.
# Using different value than ENCODE for smallRNA, as reads are short.
# They use 16% of read length.  Here that's either 4 or 8, since R1s are 26bp and R2s are mostly 50.
outFilterMatchNmin=8
###################################################

printf '\nSamples:\n%s' "$sampleNames"

printf '\n\nAlignment files will be saved to: %s\n\n' "$pathMapped"

for sn in $sampleNames;
do

    echo Mapping $sn

    infile1=$pathTrimmed$sn'_R1_001'$trimmedSuffix$extension
    infile2=$pathTrimmed$sn'_R2_001'$trimmedSuffix$extension

    # Check if you already did this one.  'Log.out' is what STAR will add to $logfull.
    logfile=$sn$mappedSuffix
    logfull=$pathMapped$logfile
    if [ -f $logfull'Log.out' ]; then echo 'File ' $logFile'Log.out' ' already exists. Skipping sample.'; continue; fi

echo Unzipping $infile1.gz
gunzip $infile1.gz
echo Unzipping $infile2.gz
gunzip $infile2.gz

# # Short RNA
    $pathSTAR --runThreadN $cores --genomeDir $genomeDir --readFilesIn $infile1 $infile2 --outFilterMismatchNoverLmax $outFilterMismatchNoverLmax --outFilterMultimapNmax $outFilterMultimapNmax --alignIntronMax $alignIntronMax --outFilterMatchNmin $outFilterMatchNmin --outFilterScoreMinOverLread $outFilterScoreMinOverLread --outFilterMatchNminOverLread $outFilterMatchNminOverLread --outFileNamePrefix $logfull

done



