#' Write expression values to csv
#'
#' Read counts from featureCounts output and write them to a csv file.  Gene-by-sample matrix.
#' @param inFiles Character - List of featureCounts output files (not the ones that end in "summary" or "report")
#' @param outPrefix String - Output filename will be this followed by ".csv" (or "_TPM.csv", if TPM=TRUE).  Include path if you don't want current directory.
#' @param tpm Logical - If true, convert values to transcripts per million
#' @param outDest String - Path to directory where output file should be written
#' @param verbose Logical - If true, reports name of each file as it's read
#' @details Note that the featureCounts "Length" field (which is used to get TPM) appears to be the summed
#' lengths of all the transcripts in the annotation file that are tagged as being part of that gene.
#' So if you used a file with all exons, it's the gene's exonic length; for introns, it's the gene's intronic length.
#' @examples
#' setwd("counts_nucseq/intronic_counts/")
#' fileList = list.files(pattern="_fcounts.txt") # note underscore prevents getting, e.g., "p200" included with "p2"
#' fileList = fileList[ - which(regexpr('summary', fileList) > 0 ) ]
#' write_expr(fileList, outPrefix="Nucseq_intronic", tpm=TRUE, outDest="tpm_nucseq")
#' @author Emma Myers
#' @export

write_expr = function(inFiles, outPrefix, tpm = FALSE, verbose=TRUE) {

    # Check arguments
    if ( !all(file.exists(inFiles)) ) {
        writeLines("Missing input files:")
        writeLines(inFiles[which(!file.exists(inFiles))])
        stop("Missing input file(s).  See above.")
    }

    # Define full output filename
    outFile = outPrefix
    if (tpm) { outFile = paste(outFile, "_TPM", sep="") }
    outFile = paste(outFile, ".csv", sep="")
    if (verbose) { writeLines(paste("Expression matrix will be saved to:", outFile, sep="\n")) }

    # Check if output file already exists
    if ( file.exists(outFile) ) {
        stop("File already exists.")
    }

    # Read in the read counts
    if (verbose) { writeLines("Reading in counts. . .") }
    fcounts = read_fcounts(inFiles, verbose=verbose)
    if (verbose) { writeLines("Done.") }

    # Get expression values into exprMat
    exprMat = fcounts$Count
    geneLens = fcounts$Length

    # Convert to TPM if requested
    if (tpm) { exprMat = counts_to_tpm(exprMat, geneLens) }

    # Write csv
    Geneid_Count = data.frame(fcounts$Geneid, fcounts$Count)
    write.table( Geneid_Count, file=outFile, row.names=FALSE, col.names=c("Geneid", colnames(fcounts$Count)), sep="," )

}






