#!/bin/bash

# cutadapt_run_paired.sh
#
# Version of cutadapt_run.sh for paired-end data in separate files.
#
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
pathRaw='/Volumes/CodingClub1/RNAseq/smallRNA/20180402_NXT/barcodes_new/Emma/raw/'
pathTrimmed='/Volumes/CodingClub1/RNAseq/smallRNA/20180402_NXT/barcodes_new/Emma/trimmed/'
extension=.fastq.gz
adapterNameR=TruSeq_Index_Rev
adapterReverse=GTTCGTCTTCTGCCGTATGCTCTANNNNNNCACTGACCTCAAGTCTGCACACGAGAAGGCTAGA
adapterNameF=TruSeq_Index
adapterForward=AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
qualityCutoff=20
minLen=20
nAdapt=3
outSuffix=''
cores=8 # number of available cores
###################################################

# Filenames must be these appended by _R1_001 or _R2_001 and extension
# sampleNames='GTAGAG_S1 GTCCGC_S2'
sampleNames='GTAGAG_S1 GTCCGC_S2 GTTTCG_S5 CGTACG_S6 GAGTGG_S7 GGTAGC_S8'


# Check that you have access to cutadapt
# (see the user guide re adding this to PATH)
cutadapt 2>/dev/null
if [ "$?" -ne "0" ]; then PATH=$PATH:$HOME/.local/bin; fi

printf '\nSamples:\n%s' "$sampleNames"
printf '\n\nTrimmed reads will be saved to: %s\n\n' "$pathTrimmed"
printf 'Output files:'

for sn in $sampleNames;
do
    # Input filenames
    r1=$sn'_R1_001'
    r2=$sn'_R2_001'
    # Get path in there yes this is terrible coding but I want to be able to display the filenames without the path
    infile1=$pathRaw$r1$extension
    infile2=$pathRaw$r2$extension
    # Output filenames
    r1out=$r1'_trimmed'$outSuffix$extension
    r2out=$r2'_trimmed'$outSuffix$extension
    r1full=$pathTrimmed$r1out
    r2full=$pathTrimmed$r2out
    reportfile=$sn$fileSuffix'_cutadapt_report.txt'
    reportfull=$pathTrimmed$reportfile

    echo $sn: $r1out, $r2out, $reportfile

    # Check if outfiles already exist
    if [ -f $r1full ]; then echo 'Trimmed file for R1 already exists. Skipping sample.'; continue; fi
    if [ -f $r2full ]; then echo 'Trimmed file for R2 already exists. Skipping sample.'; continue; fi
    if [ -f $reportfull ]; then echo 'Report file already exists.  Skipping file.'; continue; fi

    # cutadapt -a $adapterName=$adapter -a $adapterName2=$adapter2 $infile -o $outfile --trim-n -q $qualityCutoff -m $minLen -n $nAdapt -j $cores > $reportfile

    # Suddenly there's no -j option??  Is this because I'm running under ATAC user for the first time and it's using a different version of cutadapt, maybe?  It's using 1.8.1.
    cutadapt -a $adapterNameF=$adapterForward -A $adapterNameR=$adapterReverse -o $r1full -p $r2full $infile1 $infile2 --trim-n -q $qualityCutoff -m $minLen -n $nAdapt > $reportfull

done



