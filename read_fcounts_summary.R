

read_fcounts_summary = function(filenames, verbose = TRUE) {

# read_fcounts_summary.R
# Read summary output of featureCounts for specified samples into a single data frame.
# Put a few lines up here in the help text for making a quick bar plot of, say assigned reads
# per sample.  Maybe.

    # Check for nonexistent files
    msgFiles = filenames[which(!file.exists(filenames))]
    if (length(msgFiles) != 0) {
        writeLines("The following files do not exist and will be ignored:")
        writeLines(msgFiles)
        filenames = filenames[which(file.exists(filenames))]
    }

    # Read in the first file and rename the one column
    if (verbose) {writeLines(paste("Reading", filenames[1]))}
    countSummaries = read.table(filenames[1], header=TRUE, row.names=1)
    names(countSummaries) = filenames[1]
    
    # If there are more files
    if (length(filenames) > 1) {
        # Read in the rest of them
        for (f in filenames[2:length(filenames)]) {
            # Read into data frame and rename the column
            if (verbose) {writeLines(paste("Reading", f))}
            nextSummary = read.table(f, header=TRUE, row.names=1)
            names(nextSummary) = f
            # Concatenate Count field
            countSummaries = cbind(countSummaries, nextSummary)
        }
    }
    
    return(countSummaries)
    
}
