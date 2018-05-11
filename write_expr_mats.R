

# write_expr_mats.R
# 
# Point this script to a list of files, each containing read counts from featureCounts for one sample.  It will write a gene-by-sample matrix to a csv file. Either raw counts or TPM.
#
# Note that the featureCounts "Length" field (which is used to get TPM) appears to be the summed
# lengths of all the transcripts in the annotation file that are tagged as being part of that gene.
# So if you used a file with all exons, it's the gene's exonic length; for introns, it's the gene's
# intronic length.
#
# When this toolbox is a package, turn this script into a function that takes a list of files and an optional tpm argument.  Include example of how to easily define file list in help text.


# Delete when this is turned into a function and in the package
source("/Users/nelsonlab/Documents/Toolboxes/rna-seq/read_fcounts.R")
source("/Users/nelsonlab/Documents/Toolboxes/rna-seq/counts_to_tpm.R")
# source('/Volumes/CodingClub1/RNAseq/code/read_fcounts.R')
# source('/Volumes/CodingClub1/RNAseq/code/counts_to_tpm.R')
# source("/Users/work/Documents/rna-seq/read_fcounts.R")
# source("/Users/work/Documents/rna-seq/counts_to_tpm.R")

###############################################
# Stuff to set:
# Working directory, fileList, resFile, tpm
###############################################
setwd("/Users/work/Documents/counts_nucseq/intronic_counts/")
# setwd("/Volumes/DataStorage2/Emma/From_CodingClub1/RNAseq/TTX/counts/counts_m20_q20_no_comment/")
# cellType = 'EMX'; stage = 'Early'
# comparison = paste(cellType, stage, sep="")
# fileList = list.files(pattern=cellType)
# fileList = fileList[which(regexpr(stage, fileList)>0)]
# comparison = "p2"
# comparison = "Rorb_aging"
comparison = "NucSeq"
# note underscore prevents getting, e.g., "p200" included with "p2"
# fileList = list.files(pattern=paste(comparison, "_", sep=""))
# fileList = c( list.files(pattern="HTp30"), list.files(pattern="HTp200") )
fileList = list.files(pattern="_fcounts.txt")
fileList = fileList[ - which(regexpr('summary', fileList) > 0 ) ]
resFile = paste("/Users/work/Documents/", comparison, "_TPM_intronic.csv", sep="")
# resFile = paste("/Users/nelsonlab/Documents/Results_temporarily_here/RORb_results/", comparison, "_TPM.csv", sep="")
tpm = TRUE # If you want raw counts, set to FALSE
#################################################

# Check if results file already exists
if (file.exists(resFile)) {
    stop("Specified results file already exists.")
}

# Read in the read counts
writeLines("Reading in counts. . .")
fcounts = read_fcounts(fileList)
writeLines("Done.")

# Sample names
samplenames = sub("_fcounts.txt", "", fileList)

# Get TPMs into exprMat
exprMat = fcounts$Count
geneLens = fcounts$Length
if (tpm) { exprMat = counts_to_tpm(exprMat, geneLens) }

# Set matrix row and column names and write to csv
rownames(exprMat) = fcounts$Geneid
colnames(exprMat) = samplenames
write.csv(exprMat, file = resFile)


