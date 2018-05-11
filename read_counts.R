
# read_counts.R
# 
# Given a list of text files with the output of cutadapt (that gets displayed to screen unless directed to a file), return a list of read counts.
#
#
# USAGE:
# sampleCounts = read_counts( list.files(pattern='cutadapt_report.txt' ) )

# # Getting samples in the order you want for a bar plot
# names(sampleCounts) = gsub("_m20_q20_cutadapt_report.txt", "", names(sampleCounts))
# sampleCountsSorted=sampleCounts[which( regexpr('cb_', names(sampleCounts)) > 0)]
# sampleCountsSorted=c(sampleCountsSorted, sampleCounts[which( regexpr('L4_', names(sampleCounts)) > 0)])
# sampleCountsSorted=c(sampleCountsSorted, sampleCounts[which( regexpr('L6_', names(sampleCounts)) > 0)])
# sampleCountsSorted=c(sampleCountsSorted, sampleCounts[which( regexpr('BF_', names(sampleCounts)) > 0)])
# sampleCountsSorted=c(sampleCountsSorted, sampleCounts[which( regexpr('PV', names(sampleCounts)) > 0)])
# sampleCountsSorted=c(sampleCountsSorted, sampleCounts[which( regexpr('EMX', names(sampleCounts)) > 0)])
#
# # or, to use just the first two characters of each sample name:
# names(sampleCounts) = sapply( names(sampleCounts), function(x) {substr(x,1,2)} )
# sampleCountsSorted=sampleCounts[which(names(sampleCounts)=='cb')]
# sampleCountsSorted=c(sampleCountsSorted, sampleCounts[which(names(sampleCounts)=='L4')])
# sampleCountsSorted=c(sampleCountsSorted, sampleCounts[which(names(sampleCounts)=='L6')])
# sampleCountsSorted=c(sampleCountsSorted, sampleCounts[which(names(sampleCounts)=='BF')])
# sampleCountsSorted=c(sampleCountsSorted, sampleCounts[which(names(sampleCounts)=='PV')])
# sampleCountsSorted=c(sampleCountsSorted, sampleCounts[which(names(sampleCounts)=='EM')])
# sampleCountsSorted
#
# # Bar plot
# par(mar=c(3, 6, 2, 1))
# barplot(sampleCountsSorted, las = 2, cex.names=0.5)

read_counts = function(filenames) {

    # source("/Volumes/CodingClub1/RNAseq/code/file_checks.R")
    source("/Users/nelsonlab/Documents/Toolboxes/rna-seq/file_checks.R")
    
    countVec = vector()
    samplenames = filenames
    
    for (f in filenames) {
        
        # Make sure file exists
        if ( ! file_checks(f, verbose=TRUE)) { samplenames = samplenames[-which(samplenames==f)]; next }
        lines = readLines(con = f)
        # identify the line with the information we want
        lineWithCount = lines[which(regexpr("Total reads processed:", lines)>0)]
        print(lineWithCount)
        # Make sure this file actually has that line
        # if (isEmpty(lineWithCount)) {
        #     warning(paste("No line containing reads reported in", f))
        #     countVec = c(countVec, NA)
        #     next
        # }
        # Split line around spaces, and get the last piece, which contains the total read count
        splitLine = strsplit(lineWithCount, ' ')
        countStr = splitLine[[1]][length(splitLine[[1]])]
        # Remove commas, coerce to integer, and put into list
        countVec = c(countVec, as.numeric(gsub(',', '', countStr)))
    
    }
    
    names(countVec) = samplenames
    return(countVec)
    
}
