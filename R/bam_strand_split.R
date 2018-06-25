#' Split BAM with paired-end data into forward-strand and reverse-strand
#'
#' Given a list of BAM files, write two new BAM files for each, one containing forward-strand reads and the other reverse-strand reads.
#' @param inFiles Character - BAM file list
#' @param samtoolsPath String - path to samtools directory
#' @param outDest String - Save new BAM files here
#' @param verbose Logical - whether to print progress / results
#' @details
#' Requires samtools to be installed.
#' For explanations of samtools flags (like "-f 0x10" and "-f 0x20" below), see this tool:
#' http://broadinstitute.github.io/picard/explain-flags.html
#' Note that "-f" means to include only things with these flags, and "-F" means to exclude anything with these flags.
#' TIME:
#' @examples
#' bamfiles = dir('/Volumes/CodingClub1/RNAseq/smallRNA/20180402_NXT/barcodes_new/Emma/filtered_and_indexed', pattern=".bam", full.names=TRUE)
#' bamfiles = bamfiles[ -which(regexpr(".bai", bamfiles) > 0) ]
#' bamfiles = bamfiles[c(4,5,6,1,2,3)]
#' bam_pair_split(bamfiles, outDest=resPath)
#' @author Emma Myers
#' @export

bam_strand_split = function(inFiles, samtoolsPath="~/anaconda2/bin/", outDest="./", verbose=TRUE) {

    # Check that destination directory exists
    outDest = dir_check(outDest)

    # Check arguments
    if ( !all(file.exists(inFiles)) ) {
        writeLines("Missing input files:")
        writeLines(inFiles[which(!file.exists(inFiles))])
        stop("Missing input file(s).  See above.")
    }
    if (! samtoolsPath == "") { samtoolsPath = dir_check(path.expand(samtoolsPath)) }

    # Check if any output files already exist
    outFilesPlus = paste(outDest, gsub(".bam", "_plus.bam", basename(inFiles)), sep="")
    outFilesMinus = paste(outDest, gsub(".bam", "_minus.bam", basename(inFiles)), sep="")
    if ( any(file.exists(c(outFilesPlus, outFilesMinus))) ) {
        writeLines( outFilesPlus[ which(file.exists(outFilesPlus)) ] )
        writeLines( outFilesMinus[ which(file.exists(outFilesMinus)) ] )
        stop("One or more output files already exists (see filenames printed above).")
    }

    # Go through splitting files
    counter = 0
    for (f in inFiles) {

        counter = counter + 1

        if (verbose) { writeLines( paste("\nProcessing file:", f) ) }

        tStart = proc.time()[3]

        # get forward-strand reads
        if (verbose) { writeLines("Getting forward-strand reads. . .", sep="") }
        system2( paste(samtoolsPath, "samtools", sep=""), args = c("view", "-h", "-f", "0x20", f), stdout = outFilesPlus[counter] )
        if (verbose) { writeLines("done.") }

        # get reverse-strand reads
        if (verbose) { writeLines("Getting reverse-strand reads. . .", sep="") }
        system2( paste(samtoolsPath, "samtools", sep=""), args = c("view", "-h", "-f", "0x10", f), stdout = outFilesMinus[counter] )
        if (verbose) { writeLines("done.") }

        tElapsed = proc.time()[3] - tStart
        if (verbose) { writeLines(paste("done (", round(tElapsed/60, digits=2), "m).", sep="")) }

    }

}
