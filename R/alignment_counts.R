#' Get alignment counts from BAM files
#'
#' Given a list of BAM files, return an array of unique (primary) alignments.
#' @param filenames Character - BAM file list
#' @param samtoolsPath String - path to samtools directory
#' @param what String - What type of reads to include
#' @param verbose Logical - whether to print progress / results
#' @return Numeric vector? - Read count for each BAM
#' @details I'm leaving the function name ambiguous because it might make sense to give it an option to include,
#' say, multi-mapped alignments.
#' Requires samtools to be installed.
#' @examples
#' uniqueReads = alignment_counts(c("S1.bam", "S2.bam"), samtoolsPath="/opt/samtools-0.1.18")
#' @author Emma Myers
#' @export

alignment_counts = function(filenames, samtoolsPath="", what="passedQC", verbose=FALSE) {

    # Won't be necessary when this is in package
    # source("/Users/nelsonlab/Documents/Toolboxes/rna-seq/file_checks.R")
    # source("/Volumes/CodingClub1/RNAseq/code/file_checks.R")

    # What reads do we want to count
    if (what == "passedQC") {
        if (verbose) { writeLines("\nAll reads that passed quality control will be counted.") }
        args_except_file = c("view", "-f", "0", "-c")
    } else if (what == "primary") {
        if (verbose) { writeLines("\nOnly primary alignments (no multi-mappers or chimeric reads) will be counted.") }
        args_except_file = c("view", "-F", "0x904", "-c")
    }


    counts = array(dim = length(filenames))
    counter = 1

    for (f in filenames) {

        if (verbose) { writeLines("\n") }

        # Make sure file exists.  Should be checking that it's either .bam or .sam too.
        if ( ! file_checks(f, verbose=TRUE) ) { next }

        # Get primary alignment count
        if (verbose) { writeLines(paste("Processing file:", f, "...", sep=""), sep="") }
        counts[counter] = as.numeric(system2( paste(samtoolsPath, "samtools", sep=""), args = c(args_except_file, f), stdout = TRUE ))
        if (verbose) { writeLines("Done.") }

        # Increment counter
        counter = counter + 1

    }

    return(counts)

}
