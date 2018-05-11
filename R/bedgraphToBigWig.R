#' Bedgraph files to bigwig format
#'
#' Take a list of bedgraph files and write a bigwig file for each one.
#' @param filenames Character - List of bedgraph files
#' @param genomeSizeFile String - Filename (with path) to text file containing genome size information
#' @param ucscPath String - Path to UCSC tools executables
#' @param outDest String - Directory where bigwig files should be written
#' @param outSuffix String - will be appended to original filename (and ".bedgraph" is replaced by ".bigwig")
#' @details Use UCSC function bedgraphToBigWig to convert bedgraph files to bigwig format.
#' Requires UCSC tools to be installed.
#' Time:  ~40s - 1m per bedgraph, where the bedgraphs are generally ~400-800 MB.
#' Example at the command line:
#' system2('/opt/UCSCtools/wigToBigWig', args=c('samplename.gencov.bedgraph', '/path/to/genomeSizeFile.txt', samplename.bigwig))
#' @examples
#' bedgraphToBigWig(filenames, genomeSizeFile="/path/to/genomeSizeFile.txt")
#' @author Emma Myers
#' @export

bedgraphToBigWig = function(filenames, genomeSizeFile, ucscPath="/opt/UCSCtools/", outDest="./", outSuffix="") {

    # Check arguments
    ucscPath = dir_check(ucscPath)
    outDest = dir_check(outDest)
    if ( ! file.exists(genomeSizeFile) ) { stop("Specified genomeSizeFile does not exist.") }


    for (f in filenames) {

        writeLines("\n")
        # Make sure file exists, but let's not be picky about the extension as bedgraphs might not have .bedgraph
        if ( ! file_checks(f, verbose=TRUE) ) { next }

        writeLines(paste("\nProcessing file:", f))

        # Define arguments to the genomeCoverageBed command
        fOut = paste(outDest, sub(tools::file_ext(f), paste(outSuffix, "bigwig", sep=""), basename(f)), sep="")
        arguments = c(f, genomeSizeFile, fOut)

        writeLines(paste(arguments, collapse=" "))

        #################
        # Get output file
        #################
        if ( file_checks(fOut, shouldExist=FALSE, verbose=TRUE) ) {
            writeLines(paste("Creating file ", fOut, "...", sep=""), sep="")
            tStart = proc.time()[3]
            system2(paste(ucscPath, "bedgraphToBigWig", sep=""), args = arguments, stdout = fOut)
            tElapsed = proc.time()[3] - tStart
            writeLines(paste("done (", round(tElapsed/60, digits=2), "m).", sep=""))
            writeLines("Done with file.")
        }

    }

}

