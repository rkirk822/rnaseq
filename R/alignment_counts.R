#' Get alignment counts from BAM files
#'
#' Given a list of BAM files, return an array of unique (primary) alignments.
#' @param inFiles Character - BAM file list
#' @param samtoolsPath String - path to samtools directory
#' @param what String - What type of reads to include
#' @param verbose Logical - whether to print progress / results
#' @return Numeric vector? - Read count for each BAM
#' @details Might make sense to give it an option to include, say, multi-mapped alignments.
#' Requires samtools to be installed.
#' TIME:  ~1m / BAM.
#' @examples
#' readCounts = alignment_counts(bams, verbose=TRUE)
#' pdf("countsBar.pdf")
#' barplot(readCounts/1000000, las=2, main="Aligned reads", cex.names = 0.5, names.arg=gsub("_mapped_Aligned.sortedByCoord.out.bam", "", bams), ylab="Millions")
#' dev.off()
#' @author Emma Myers
#' @export

alignment_counts = function(inFiles, samtoolsPath="~/anaconda2/bin/", what="passedQC", verbose=FALSE) {

    # Check arguments
    if ( !all(file.exists(inFiles)) ) {
        writeLines("Missing input files:")
        writeLines(inFiles[which(!file.exists(inFiles))])
        stop("Missing input file(s).  See above.")
    }
    if (! samtoolsPath == "") { samtoolsPath = dir_check(path.expand(samtoolsPath)) }

    # What reads do we want to count
    if (what == "passedQC") {
        if (verbose) { writeLines("\nAll reads that passed quality control will be counted.") }
        args_except_file = c("view", "-f", "0", "-c")
    } else if (what == "primary") {
        if (verbose) { writeLines("\nOnly primary alignments (no multi-mappers or chimeric reads) will be counted.") }
        args_except_file = c("view", "-F", "0x904", "-c")
    }


    counts = array(dim = length(inFiles))
    counter = 1

    for (f in inFiles) {

        # Get primary alignment count
        if (verbose) { writeLines(paste("\nProcessing file:", f, "...", sep=""), sep="") }
        tStart = proc.time()[3]
        counts[counter] = as.numeric(system2( paste(samtoolsPath, "samtools", sep=""), args = c(args_except_file, f), stdout = TRUE ))
        tElapsed = proc.time()[3] - tStart
        if (verbose) { writeLines(paste("done (", round(tElapsed/60, digits=2), "m).", sep="")) }

        # Increment counter
        counter = counter + 1

    }

    return(counts)

}
