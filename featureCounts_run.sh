#!/bin/bash

# featureCounts_run.sh
# Run though a list of BAM files with alignments and get gene read counts with featureCounts.
# 
# TIME: 
# ~1.5m per file, for the ones I noticed.  I think it says the time in the _display.txt files.

# Stuff to set
annotFile=/Volumes/CodingClub1/STAR_stuff/annotations/mm10_refGene.gtf
outPrefix=/Volumes/CodingClub1/RNAseq/TTX/counts/counts_m20_q20/
outSuffix=_fcounts.txt
inPrefix=/Volumes/CodingClub1/RNAseq/TTX/mapped/mapped_m20_q20/
inSuffix=_mapped_Aligned.sortedByCoord.out.bam
# inSuffix=Aligned.sortedByCoord.out.bam
nCores=8
outSuffixDisp=_display.txt # stuff printed to screen will go to a file instead
sampleNames='EMXCtlEarly_1 EMXCtlEarly_2 EMXCtlEarly_3 EMXCtlEarly_4 EMXCtlLate_1 EMXCtlLate_2 EMXCtlLate_3 EMXCtlLate_4 EMXTTXEarly_1 EMXTTXEarly_2 EMXTTXEarly_3 EMXTTXLate_1 EMXTTXLate_2 EMXTTXLate_3 EMXTTXLate_4 PVCtlEarly_1 PVCtlEarly_2 PVCtlEarly_3 PVCtlLate_1 PVCtlLate_2 PVCtlLate_3 PVCtlLate_4 PVTTXEarly_1 PVTTXEarly_2 PVTTXEarly_3 PVTTXLate_1 PVTTXLate_2 PVTTXLate_3 PVTTXLate_4'
# sampleNames='EMXCtlEarly_1 EMXCtlEarly_2'
# sampleNames='BF_RORbHTp200_1 BF_RORbHTp200_2 BF_RORbHTp200_3 BF_RORbHTp200_4 BF_RORbHTp2_1 BF_RORbHTp2_2 BF_RORbHTp2_3 BF_RORbHTp2_4 BF_RORbHTp30_1 BF_RORbHTp30_2 BF_RORbHTp30_3 BF_RORbHTp30_4 BF_RORbHTp7_1 BF_RORbHTp7_2 BF_RORbHTp7_3 BF_RORbHTp7_4 BF_RORbKOp2_1 BF_RORbKOp2_2 BF_RORbKOp2_3 BF_RORbKOp2_4 BF_RORbKOp30_1 BF_RORbKOp30_2 BF_RORbKOp30_3 BF_RORbKOp30_4 BF_RORbKOp7_1 BF_RORbKOp7_2 BF_RORbKOp7_3 BF_RORbKOp7_4'
# sampleNames=BF_RORbHTp2_1


# Check that we have access to featureCounts
featureCounts 2>/dev/null
if [ "$?" -ne "0" ]; then PATH=$PATH:/opt/subread-1.6.0-MacOSX-x86_64/bin; fi


for sn in $sampleNames;
do
	echo Getting counts for $sn
	outputFile=$outPrefix$sn$outSuffix
	inputFile=$inPrefix$sn$inSuffix
	outputFileDisp=$outPrefix$sn$outSuffixDisp

	# Check output file doesn't already exist
	if [ -f $outputFile ]; then echo 'Counts file already exists.  Skipping file.'; continue; fi

	echo Input $inputFile
	echo Output $outputFile
	featureCounts -a $annotFile -o $outputFile $inputFile -T $nCores 2> $outputFileDisp
done
