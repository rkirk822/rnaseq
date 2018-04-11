

# intronic_read_percent.R
# 
# In a given dataset, get the percentage of intronic reads out of the total read count.
#
# You need the files output by featureCounts that end in "txt.summary", one set produced using
# an annotation file containing only non-intronic regions and one containing only introns.
#
# And you need two text files which, on being read in with read.table, have in the fourth
# column the widths of all the intronic ranges and all the non-intronic ranges, respectively.
# Non-overlapping ranges.  I make these files using write_reduced_ranges.R.


source("/Users/nelsonlab/Documents/Toolboxes/rna-seq/read_fcounts_summary.R")
source("/Users/nelsonlab/Documents/Toolboxes/rna-seq/file_checks.R")
# source('/Users/emmamyers/Documents/Work_temp/read_fcounts_summary.R')

# # int = intronic; non = non-intronic; one annot includes only intronic regions and the other everything else
annotNon = "/Users/nelsonlab/Documents/Results_temporarily_here/NucSeq_results/mm10_refSeq_exons_filtered_sorted_reduced.txt"
annotInt = "/Users/nelsonlab/Documents/Results_temporarily_here/NucSeq_results/mm10_refSeq_introns_filtered_sorted_reduced.txt"
# pathNon = "/Users/nelsonlab/Documents/Smintheus_stuff_copied_here/RORb_stuff/counts_m20_q20/"
# pathInt = "/Users/nelsonlab/Documents/Smintheus_stuff_copied_here/RORb_stuff/intronic_counts/"
# dataname = "Rorb"
# pathInt = "/Users/nelsonlab/Documents/Smintheus_stuff_copied_here/NucSeq_stuff/nucseq_counts/intronic_counts"
# pathNon = "/Users/nelsonlab/Documents/Smintheus_stuff_copied_here/NucSeq_stuff/nucseq_counts/usual_counts/"
# dataname = "Nuc-seq"
pathInt = "/Users/nelsonlab/Documents/Smintheus_stuff_copied_here/NucSeq_stuff/rnaseq_counts/intronic_counts"
pathNon = "/Users/nelsonlab/Documents/Smintheus_stuff_copied_here/NucSeq_stuff/rnaseq_counts/usual_counts/"
dataname = "RNA-seq in same cell types as Nuc-seq"


rangesEx = read.table(annotNon)
totalLenEx = sum(rangesEx[,4])
rangesInt = read.table(annotInt)
totalLenInt = sum(rangesInt[,4])


# Get non-intronic read count
setwd(pathNon)
writeLines(paste("Reading non-intronic counts for dataset", dataname))
filesNon = list.files(pattern="summary")
summaryNon = read_fcounts_summary(filesNon, verbose=FALSE)
names(summaryNon) = gsub("_fcounts.txt.summary", "", names(summaryNon))


# Get intronic read count
setwd(pathInt)
writeLines(paste("Reading intronic counts for dataset", dataname))
filesInt = list.files(pattern="summary")
summaryInt = read_fcounts_summary(filesInt, verbose=FALSE)
names(summaryInt) = gsub("_intronic_fcounts.txt.summary", "", names(summaryInt))

# NORMALIZE TO TOTAL INTRONIC / EXONIC LENGTH
summaryIntNorm = summaryInt / totalLenInt
summaryNonNorm = summaryNon / totalLenEx

# Get percentage of intronic reads
totals = summaryNonNorm + summaryIntNorm
percInt = round( (summaryIntNorm / totals) * 100)
percInt = percInt[ -which(is.nan(rowSums(percInt))) , ]

# Plot - bars are samples; heights are percentages
par(mar=c(6, 6, 2, 1))
percIntAssigned = as.vector(percInt[1,], mode="numeric")
names(percIntAssigned) = names(percInt)
barplot(percIntAssigned, main=paste("Intronic reads in dataset", dataname), las=2, ylim=c(0,100), ylab="Percentage")

# Write to table, with columns:
# Sample, non-intronic count, intronic count, intronic percentage
tablename = paste(dataname, "_intron_counts.csv", sep="")
if (file_checks(tablename, shouldExist=FALSE, verbose=TRUE)) {
    writeLines(paste("Writing values to ", tablename, "...", sep=""))
    # Information into data frame
    countInfo = data.frame(names(totals), as.vector(summaryNon[1,], mode="numeric"), as.vector(summaryInt
    [1,], mode="numeric"), as.vector(summaryNonNorm[1,], mode="numeric"), as.vector(summaryIntNorm[1,], mode="numeric"), as.vector(percIntAssigned, mode="numeric"))
    # Set column names
    colnames(countInfo)=c("Sample", "Non-intronic count", "Intronic count", "Normalized non-intronic", "Normalized intronic", "Percent intronic (based on normalized values)")
    # Write to file
    write.table(countInfo, row.names=FALSE, file=tablename, sep=",")
    writeLines("done.")
}



