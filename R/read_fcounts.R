#' Read in raw counts
#'
#' Read counts from featureCounts output into a dataframe.
#' @param filenames Character - List of featureCounts output files (not the ones that end in "summary" or "report")
#' @param verbose Logical - If true, reports name of each file as it's read
#' @details  Output dataframe has columns corresponding to featureCounts output files, and rows corresponding to features (usually genes).
#' The column names of the output's Count field are the full filenames, so you might want to do something like:
#' counts = read_fcounts(fileList)
#' colnames(counts$Count) = gsub("_fcounts.txt", "", colnames(counts$Count))
#' @examples
#' example here
#' @author Emma Myers
#' @export

read_fcounts = function(filenames, verbose = TRUE) {

    # Check for nonexistent files
    msgFiles = filenames[which(!file.exists(filenames))]
    if (length(msgFiles) != 0) {
        writeLines('The following files do not exist and will be ignored:')
        writeLines(msgFiles)
        filenames = filenames[which(file.exists(filenames))]
    }

    # Read in the first file and rename last column
    if (verbose) {writeLines(paste('Reading', filenames[1]))}
    fcounts = read.table(filenames[1], header=TRUE)
    names(fcounts)[7] = 'Count'

    # Read in the rest of the files
    for (f in filenames[2:length(filenames)]) {
        # Read into data frame and make sure all fields except Count match up
        if (verbose) {writeLines(paste('Reading', f))}
        nextCounts = read.table(f, header=TRUE)
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

    colnames(fcounts$Count) = filenames
    rownames(fcounts$Count) = fcounts$Geneid
    return(fcounts)

}
