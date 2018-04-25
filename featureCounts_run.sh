#!/bin/bash

# featureCounts_run.sh
# Run though a list of BAM files with alignments and get gene read counts with featureCounts.
# 
# TIME: 
# ~1.5m per file, for the ones I noticed.  I think it says the time in the _display.txt files.

# # Stuff to set

annotFile=/Volumes/CodingClub1/STAR_stuff/annotations/mm10_refGene.gtf

# echo Using intronic annotation file
# annotFile=/Volumes/CodingClub1/STAR_stuff/annotations/mm10_refSeq_introns_geneids.gtf

outPrefix=/Volumes/CodingClub1/RNAseq/smallRNA/20180402_NXT/barcodes_new/Emma/counts/exonic_multi/
outSuffix=_fcounts.txt
inPrefix=/Volumes/CodingClub1/RNAseq/smallRNA/20180402_NXT/barcodes_new/Emma/mapped/
# inSuffix=_Aligned.sortedByCoord.out.bam
inSuffix=_Aligned.out.sam
nCores=8
outSuffixDisp=_display.txt # stuff printed to screen will go to a file instead

# sampleNames='GTAGAG_S1 GTCCGC_S2'
sampleNames='GTAGAG_S1 GTCCGC_S2 GTTTCG_S5 CGTACG_S6 GAGTGG_S7 GGTAGC_S8'
# sampleNames='cb_p036_1 cb_p036_2 cb_p036_3 cb_p036_4 L4_RORb_1 L4_RORb_2 L4_RORb_3 L4_RORb_4'
# sampleNames='cb_p033_4 cb_p033_5 cb_p033_6 cb_p033_7 L4_RORht_1 L4_RORht_2 L4_RORht_3 L4_RORht_4'
# sampleNames='EMXCtlEarly_1 EMXCtlEarly_2 EMXCtlEarly_3 EMXCtlEarly_4 EMXTTXEarly_1 EMXTTXEarly_2 EMXTTXEarly_3 EMXCtlLate_1 EMXCtlLate_2 EMXCtlLate_3 EMXCtlLate_4 EMXTTXLate_1 EMXTTXLate_2 EMXTTXLate_3 EMXTTXLate_4 PVCtlEarly_1 PVCtlEarly_2 PVCtlEarly_3 PVTTXEarly_1 PVTTXEarly_2 PVTTXEarly_3 PVCtlLate_1 PVCtlLate_2 PVCtlLate_3 PVCtlLate_4 PVTTXLate_1 PVTTXLate_2 PVTTXLate_3 PVTTXLate_4'
# sampleNames='BF_RORbHTp2_1 BF_RORbHTp2_2 BF_RORbHTp2_3 BF_RORbHTp2_4 BF_RORbKOp2_1 BF_RORbKOp2_2 BF_RORbKOp2_3 BF_RORbKOp2_4 BF_RORbHTp7_1 BF_RORbHTp7_2 BF_RORbHTp7_3 BF_RORbHTp7_4 BF_RORbKOp7_1 BF_RORbKOp7_2 BF_RORbKOp7_3 BF_RORbKOp7_4 BF_RORbHTp30_1 BF_RORbHTp30_2 BF_RORbHTp30_3 BF_RORbHTp30_4 BF_RORbKOp30_1 BF_RORbKOp30_2 BF_RORbKOp30_3 BF_RORbKOp30_4 BF_RORbHTp200_1 BF_RORbHTp200_2 BF_RORbHTp200_3 BF_RORbHTp200_4'


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

# featureCounts -a $annotFile -O -o $outputFile $inputFile -T $nCores 2> $outputFileDisp

echo ALLOWING MULTIMAPPERS
featureCounts -a $annotFile -M -o $outputFile $inputFile -T $nCores 2> $outputFileDisp
done
