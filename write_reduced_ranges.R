


# PUT IN CHECK TO AVOID OVERWRITING FILE and point out if bedfile doesn't exist

# write_reduced_ranges.R
# 
# Given a bedfile, collapse it into non-overlapping ranges and write them to a table.

library(GenomicRanges)

bedfile = "/Users/nelsonlab/Documents/Results_temporarily_here/NucSeq_results/mm10_refSeq_exons_filtered_sorted.bed"
outfile = "mm10_refSeq_exons_filtered_sorted_reduced.txt"

# Get total intron and total exon lengths
writeLines("Reading bedfile info...", sep="")
featureTable =  read.table(bedfile, header=FALSE)
grObject = with(featureTable, GRanges(V1, IRanges(V2,V3))) # GRanges object
writeLines("collapsing to non-overlapping ranges...", sep="")
grObject = reduce(grObject)
writeLines("writing to file...", sep="")
write.table(as.data.frame(grObject)[,1:4], outfile, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
writeLines("done.")
