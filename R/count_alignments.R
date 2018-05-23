#' Get alignment counts from BAM files
#'
#' Given a list of BAM files, return an array of either unique (primary) alignments or all alignments that passed QC.
#' @param inFiles Character - BAM file list
#' @param samtoolsPath String - path to samtools directory
#' @param what String - What type of reads to include
#' @param verbose Logical - whether to print progress / results
#' @param outFile String - Name of file to write counts to, or NULL (this function can take awhile, and there's no safeguard against overwriting this file here)
#' @return Numeric vector? - Alignment count for each BAM
#' @details Might make sense to give it an option to include, say, multi-mapped alignments.
#' Requires samtools to be installed.
#' For explanations of samtools flags (like "-f 0" and "-F 0x904" below), see this tool:
#' http://broadinstitute.github.io/picard/explain-flags.html
#' Note that "-f" means to include only things with these flags, and "-F" means to exclude anything with these flags.
#' TIME:  ~1m / BAM.
#' @examples
#' readCounts = count_alignments(bamFileNames)
#' pdf("countsBar.pdf")
#' barplot(readCounts/1000000, las=2, main="Uniquely mapped reads", cex.names = 0.5, names.arg=gsub("_mapped_Aligned.sortedByCoord.out.bam", "", bamFileNames), ylab="Millions")
#' dev.off()
#' @author Emma Myers
#' @export

count_alignments = function(inFiles, samtoolsPath="~/anaconda2/bin/", what="primary", verbose=TRUE, outFile="alignmentCounts.txt") {

    # Check arguments
    if ( !all(file.exists(inFiles)) ) {
        writeLines("Missing input files:")
        writeLines(inFiles[which(!file.exists(inFiles))])
        stop("Missing input file(s).  See above.")
    }
    if (! samtoolsPath == "") { samtoolsPath = dir_check(path.expand(samtoolsPath)) }

    # What reads do we want to count
    if (what == "passedQC") {
        if (verbose) { writeLines("\nAll alignments that passed quality control will be counted (including those from a single multi-mapping read!") }
        args_except_file = c("view", "-f", "0", "-c")
    } else if (what == "primary") {
        if (verbose) { writeLines("\nOnly primary alignments (no multi-mappers or chimeric reads) will be counted.") }
        args_except_file = c("view", "-F", "0x904", "-c")
    }


    counts = array(dim = length(inFiles))
    counter = 1

    for (f in inFiles) {

        # Get primary alignment count
        if (verbose) { writeLines(paste("\nProcessing file: ", f, "...", sep=""), sep="") }
        tStart = proc.time()[3]
        counts[counter] = as.numeric(system2( paste(samtoolsPath, "samtools", sep=""), args = c(args_except_file, f), stdout = TRUE ))
        tElapsed = proc.time()[3] - tStart
        if (verbose) { writeLines(paste("done (", round(tElapsed/60, digits=2), "m).", sep="")) }

        # Increment counter
        counter = counter + 1

    }

    # Things that are likely to be in the input file names
    n = gsub("_Aligned", "", basename(inFiles))
    n = gsub(".sortedByCoord", "", n)
    n = gsub(".out.bam", "", n)
    names(counts) = n

    if (!is.null(outFile)) {
        if (verbose) { writeLines(paste( "\nAlignment counts will be written to", outFile )) }
        write.table(counts, file=outFile, col.names=FALSE, quote=FALSE)
    }

    return(counts)

}
