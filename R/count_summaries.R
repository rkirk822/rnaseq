#' Get read counts from featureCounts summary files
#'
#' Given a list of featureCounts summary files, return a dataframe of read counts broken down by type (assigned, unassigned due to multimapping, etc)
#' @param inFiles Character - List of featureCounts output files
#' @param verbose Logical - whether to report progress
#' @return Data frame - Columns correspond to input files and rows to read types
#' @details Gets read counts per sample from the "summary.txt" files output by featureCounts.  This means you're getting counts for the features included
#' in the annotation file given to featureCounts.  This is how to get specifically intronic read counts, by giving featureCounts a gtf with only introns
#' and then giving this function the resulting summary files.
#' TIME:
#'
#' @examples
#' summaryFiles = list.files("NucSeq/counts/exonic_counts", pattern="summary")
#' fCounts = count_summaries(paste("NucSeq/counts/exonic_counts/", summaryFiles, sep=""))
#' assignedCounts = as.vector(fCounts[1,], mode="numeric")
#' pdf("NucSeq_exonic_primary_read_counts.pdf")
#' barplot(assignedCounts/1000000, las=2, main="Nuc-seq exonic primary reads", cex.names = 0.5, names.arg=gsub("_fcounts.txt.summary", "", summaryFiles), ylab="Millions")
#' dev.off()
#' @author Emma Myers
#' @export

count_summaries = function(inFiles, verbose = TRUE) {

    # Check arguments
    if ( !all(file.exists(inFiles)) ) {
        writeLines("Missing input files:")
        writeLines(inFiles[which(!file.exists(inFiles))])
        stop("Missing input file(s).  See above.")
    }

    # Read in the first file
    # There's probably a nicer way to do this, but this makes the dataframe end up the way I want it whereas initializing an empty one doesn't
    if (verbose) {writeLines(paste("Reading", inFiles[1]))}
    countSummaries = read.table(inFiles[1], header=TRUE, row.names=1)

    # If there are more files, read in the rest of them and concatenate
    if (length(inFiles) > 1) {
        for (f in inFiles[2:length(inFiles)]) {
            if (verbose) {writeLines(paste("Reading", f))}
            countSummaries = cbind(countSummaries, read.table(f, header=TRUE, row.names=1))
        }
    }

    # Nicer column names than the full filenames that might include paths
    names(countSummaries) = gsub(".summary", "", basename(inFiles))
    # Get rid of "_fcounts.txt" in the names separately, because it's not automatically done by featureCounts and might not be there
    names(countSummaries) = gsub("_fcounts.txt", "", names(countSummaries))

    return(countSummaries)

}

