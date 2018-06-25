#' Reduce genomic ranges to set of non-overlapping ranges
#'
#' Given a bedfile containing genomic ranges, collapse overlapping ranges using the GenomicRanges package.
#' @param inFile String - Name of bedfile containing genomic ranges
#' @param outFile String - Name to give output bedfile with non-overlapping genomic ranges
#' @return reducedRanges
#' @details
#' TIME: 10-20s, mostly to read the input file.
#' @examples
#' bed1 = "/Users/nelsonlab/Documents/Results_temporarily_here/NucSeq_results/what_was_this/mm10_refSeq_exons_filtered_sorted.bed"
#' ranges = genomic_ranges_reduce(bed1, outFile="~/Documents/mm10_refSeq_exons_filtered_sorted_reduced.txt")
#' @author Emma Myers
#' @export

genomic_ranges_reduce = function(inFile, outFile=NULL) {

    # Check arguments
    if ( !file.exists(inFile) ) { stop("Input file does not exist.") }
    if ( !is.null(outFile) && file.exists(outFile) ) { stop("Output file already exists.") }

    # Get total intron and total exon lengths
    writeLines("Reading bedfile...", sep="")
    featureTable =  read.table(inFile, header=FALSE)
    # Define GRanges object
    grObject = with(featureTable, GenomicRanges::GRanges(V1, IRanges::IRanges(V2,V3)))
    writeLines("collapsing to non-overlapping ranges...", sep="")
    grObject = IRanges::reduce(grObject)
    reducedRanges = as.data.frame(grObject)[,1:4]
    writeLines("done.")
    if ( !is.null(outFile) ) {
        writeLines(paste("Reduced ranges will be written to", outFile))
        write.table(reducedRanges, outFile, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
    }

    return(reducedRanges)

}
