#' Get read counts from fastqs
#'
#' Given a list of fastq files, return an array of read counts.
#' @param inFiles Character - List of fastqs
#' @param verbose Logical - whether to print progress / results
#' @return Numeric vector? - Read count for each input file
#' @details If you use (unzipped!) fastqs, this function will call wc -l through the command line to get the number of reads in the file,
#' and then divide by 4 to get the read count (fastqs include 4 lines for each read).
#' TIME:
#'
#' @examples
#' readCounts = count_reads(reportFiles)
#' pdf("countsBar.pdf")
#' barplot(readCounts/1000000, las=2, main="Total reads", cex.names = 0.5, names.arg=gsub("_mapped_Aligned.sortedByCoord.out.bam", "", bams), ylab="Millions")
#' dev.off()
#' @author Emma Myers
#' @export

count_reads = function(inFiles, what="all", verbose=FALSE) {

    # Check arguments
    if ( !all(file.exists(inFiles)) ) {
        writeLines("Missing input files:")
        writeLines(inFiles[which(!file.exists(inFiles))])
        stop("Missing input file(s).  See above.")
    }
    # Make sure they're not zipped, which seems to result in an inaccurate count without any kind of error
    if ( any( substr(inFiles, nchar(inFiles)-2, nchar(inFiles)) == ".gz" ) ) {
        stop("At least one input filename ends in \".gz\".  Applying this to zipped files may yield inaccurate counts.")
    }

    # Get counts
    countVec = vector(mode="numeric", length=length(inFiles))
    counter = 1
    for (f in inFiles) {
        if (verbose) { writeLines(paste("\nProcessing file: ", f, "...", sep=""), sep="") }
        tStart = proc.time()[3]
        # Output of wc -l filename.fastq into variable
        wcOut = system2("wc", args = c("-l", f), stdout=TRUE)
        # wcOut is a string consisting of some number of spaces, followed by the line count, followed by the filename.
        # Split wcOut around spaces.  The first element of the result (at [[1]]) contains the split string.
        # The first element of the split string with length > 0 will be the number of lines.
        countIdx = which( nchar( strsplit(wcOut[1], split=" ")[[1]] ) > 0 )[1]
        lineCount = as.numeric(strsplit(wcOut[1], split=" ")[[1]][countIdx])
        # Each read takes up 4 lines
        countVec[counter] = lineCount/4
        counter = counter + 1
        if (verbose) { writeLines(paste("done (", round( (proc.time()[3] - tStart) /60, digits=2), "m).", sep="")) }
    }

    names(countVec) = basename(gsub(".fastq", "", inFiles))
    return(countVec)

}
