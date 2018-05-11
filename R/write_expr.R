#' Write expression values to csv
#'
#' Read counts from featureCounts output and write them to a csv file.  Gene-by-sample matrix.
#' @param filenames Character - List of featureCounts output files (not the ones that end in "summary" or "report")
#' @param outPrefix String - Output filename (sans path) will be this followed by ".csv" or "_TPM.csv", if TPM=TRUE
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

write_expr = function(filenames, outPrefix, tpm = FALSE, outDest="./", verbose = TRUE) {

    # Check arguments
    outDest = dir_check(outDest)

    # Define output file name
    fOut = paste(outDest, outPrefix, sep="")
    if (tpm) { fOut = paste(fOut, "_TPM", sep="") }
    fOut = paste(fOut, ".csv", sep="")
    writeLines(fOut)

    # Check if output file already exists
    if ( !file_checks(fOut, shouldExist=FALSE, verbose=TRUE) ) {
        stop("File already exists.")
    }

    # Check for nonexistent files
    msgFiles = filenames[which(!file.exists(filenames))]
    if (length(msgFiles) != 0) {
        writeLines('The following files do not exist and will be ignored:')
        writeLines(msgFiles)
        filenames = filenames[which(file.exists(filenames))]
    }

    # Read in the read counts
    writeLines("Reading in counts. . .")
    fcounts = read_fcounts(fileList, verbose=verbose)
    writeLines("Done.")

    # Get expression values into exprMat
    exprMat = fcounts$Count
    geneLens = fcounts$Length

    # Convert to TPM if requested
    if (tpm) { exprMat = counts_to_tpm(exprMat, geneLens) }

    # Set matrix row and column names and write to csv
    rownames(exprMat) = fcounts$Geneid
    samplenames = sub("_fcounts.txt", "", fileList)
    colnames(exprMat) = samplenames
    write.csv(exprMat, file = fOut)

}






