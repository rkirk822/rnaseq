#' Read in raw counts for genomic features
#'
#' Read counts from featureCounts output into a dataframe.
#' @param inFiles Character - List of featureCounts output files (not the ones that end in "summary" or "report")
#' @param verbose Logical - If true, reports name of each file as it's read
#' @details  Output dataframe has columns corresponding to featureCounts output files, and rows corresponding to features (usually genes).
#' @examples
#' fcs = dir(paste(projectPath, "counts/exonic_counts/", sep=""), pattern="_fcounts.txt")
#' fcs = fcs[ which(regexpr("summary", fcs) < 0) ]
#' counts_nucseq_ex = count_features(paste(projectPath, "counts/exonic_counts/",fcs,sep=""))
#' dim(counts_nucseq_ex)
#' [1] 24746     7
#' @author Emma Myers
#' @export

count_features = function(inFiles, verbose = TRUE) {

    # Check arguments
    if ( !all(file.exists(inFiles)) ) {
        writeLines("Missing input files:")
        writeLines(inFiles[which(!file.exists(inFiles))])
        stop("Missing input file(s).  See above.")
    }

    # Read in the first file and rename last column
    if (verbose) {writeLines(paste('Reading', inFiles[1]))}
    fcounts = read.table(inFiles[1], header=TRUE, stringsAsFactors=FALSE)
    names(fcounts)[7] = 'Count'
    # Count field has to be dataframe instead of matrix, or having only one file will fuck it up
    fcounts$Count = as.data.frame(fcounts$Count)

    # Read in the rest of the files
    if ( length(inFiles) > 1 ) {
        for (f in inFiles[2:length(inFiles)]) {
            # Read into data frame and make sure all fields except Count match up
            if (verbose) {writeLines(paste('Reading', f))}
            nextCounts = read.table(f, header=TRUE, stringsAsFactors=FALSE)
            names(nextCounts)[7] = 'Count'
            if (any(fcounts$Geneid != nextCounts$Geneid)) {'Mismatched Geneid fields; skipping'; continue}
            if (any(fcounts$Chr != nextCounts$Chr)) {'Mismatched Chr fields; skipping'; continue}
            if (any(fcounts$Start != nextCounts$Start)) {'Mismatched Start fields; skipping'; continue}
            if (any(fcounts$End != nextCounts$End)) {'Mismatched End fields; skipping'; continue}
            if (any(fcounts$Strand != nextCounts$Strand)) {'Mismatched Strand fields; skipping'; continue}
            if (any(fcounts$Length != nextCounts$Length)) {'Mismatched Length fields; skipping'; continue}
            # Concatenate Count field
            fcounts$Count = cbind(fcounts$Count, nextCounts$Count)
        }
    }

    # Give the expression matrix row names
    rownames(fcounts$Count) = fcounts$Geneid

    # Use filenames (sans path) as column names
    # (and if featureCounts output has been consistently named we can clean it up a little)
    colnames(fcounts$Count) = gsub("_fcounts.txt", "", basename(inFiles))
    return(fcounts)

}
