
# read_counts.R
# 
# Given a list of text files with the output of cutadapt (that gets displayed to screen unless directed to a file), return a list of read counts.
#
#
# USAGE:
# sampleCounts = read_counts( list.files(pattern='cutadapt_report.txt' ) )
# names(sampleCounts) = sapply( names(sampleCounts), function(x) {substr(x,1,2)} )
# barplot(sampleCounts, las = 2)

read_counts = function(filenames) {

    source('/Users/nelsonlab/Documents/Toolboxes/rna-seq/file_checks.R')

    countVec = vector()
    samplenames = filenames
    
    for (f in filenames) {
    
        # Make sure file exists
        if ( ! file_checks(f)) { samplenames = samplenames[-which(samplenames==f)]; next }
        lines = readLines(con = f)
        # identify the line with the information we want
        lineWithCount = lines[which(regexpr("Total reads processed:", lines)>0)]
        # Split line around spaces, and get the last piece, which contains the total read count
        splitLine = strsplit(lineWithCount, ' ')
        countStr = splitLine[[1]][length(splitLine[[1]])]
        # Remove commas, coerce to integer, and put into list
        countVec = c(countVec, as.numeric(gsub(',', '', countStr)))
        
    
    }
    
    names(countVec) = samplenames
    return(countVec)
    
}
