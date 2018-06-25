#' Get counts or TPM matrix
#'
#' Read counts from featureCounts output into a gene-by-sample matrix.  Optionally convert to TPM; optionally write to a csv file.
#' @param inFiles Character - List of featureCounts output files (not the ones that end in "summary" or "report")
#' @param tpm Logical - If true, convert values to transcripts per million
#' @param outFileNamePrefix String - If given, expression matrix will be written to csv with filename outFileNamePrefix.csv (or "_TPM.csv", if TPM=TRUE).
#' @param verbose Logical - If true, announce output filename (if there is one).
#' @return exprMat
#' @details Note that the featureCounts "Length" field (which is used to get TPM) appears to be the summed
#' lengths of all the transcripts in the annotation file that are tagged as being part of that gene.
#' So if you used a file with all exons, it's the gene's exonic length; for introns, it's the gene's intronic length.
#' @examples
#' This way of creating the input file list lets you put the columns of the matrix in the order you want
#' comparisons=c("HTp2_","HTp7_","HTp30_", "KOp2_", "KOp7_", "KOp30_") # underscore at the end prevents confusing p2 with p200
#' countFiles=vector(mode="character")
#' for (c in comparisons) { countFiles=c(countFiles, dir("RORb/counts/counts_m20_q20", pattern=c, full.names=TRUE)) }
#' countFiles = countFiles[ -which( regexpr("summary", countFiles) > 0 ) ]
#' countFiles = countFiles[ -which( regexpr("display", countFiles) > 0 ) ]
#' rorbTPM = get_expr(countFiles, tpm=TRUE, outFileNamePrefix="~/Documents/RORb", verbose=FALSE)
#' @author Emma Myers
#' @export

get_expr = function(inFiles, tpm = FALSE, outFileNamePrefix=NULL, verbose=TRUE) {

    # Check arguments
    if ( !all(file.exists(inFiles)) ) {
        writeLines("Missing input files:")
        writeLines(inFiles[which(!file.exists(inFiles))])
        stop("Missing input file(s).  See above.")
    }

    # Read in the read counts
    if (verbose) { writeLines("Reading in counts. . .") }
    fcounts = count_features(inFiles, verbose=verbose)
    if (verbose) { writeLines("Done.") }

    # Get expression values into exprMat (distinct from Count, which we may or may not be transforming)
    exprMat = fcounts$Count
    geneLens = fcounts$Length

    # Convert to TPM if requested
    if (tpm) { exprMat = counts_to_tpm(exprMat, geneLens) }

    # Write csv if requested
    if ( !is.null(outFileNamePrefix) ) {

         # Define full output filename
        outFile = outFileNamePrefix
        if (tpm) { outFile = paste(outFile, "_TPM", sep="") }
        outFile = paste(outFile, ".csv", sep="")

        # Check if output file already exists before writing
        if ( file.exists(outFile) ) {
            stop( paste("File already exists: ", outFile) )
        } else {
            if (verbose) { writeLines(paste( "Writing expression matrix to", outFile) ) }
            write.csv(exprMat, file=outFile)
        }
    }


    return(exprMat)

}






