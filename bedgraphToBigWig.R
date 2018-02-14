

# bedgraphToBigWig.R
#
# Take a list of bedgraph files and make bigwigs, using UCSC's bedgraphToBigWig function.
#
#
# USAGE:
# > bedgraphToBigWig(filenames, ucscPath="/opt/UCSC/", genomeSizeFile="/path/to/genomeSizeFile.txt")
#
# TIME:
# On yasuimac: ~1m per file.
#
#
########################################################
# Here are the equivalent commands at the command line:
########################################################
# $ system2('/opt/UCSC/wigToBigWig', args=c('samplename.gencov.bedgraph', '/path/to/genomeSizeFile.txt', samplename.bigwig))


bedgraphToBigWig = function(filenames, ucscPath, genomeSizeFile, bedgraphDest="./") {

    # Won't be necessary when this is in package
    source("/Users/nelsonlab/Documents/Toolboxes/code/file_checks.R")
    source("/Users/nelsonlab/Documents/Toolboxes/code/dir_check.R")
    library(tools)

    # Check arguments
    ucscPath = dir_check(ucscPath)
    bwDest = dir_check(bwDest)
    if ( ! file.exists(genomeSizeFile) ) { stop("Specified genomeSizeFile does not exist.") }
    
    
    for (f in filenames) {
        
        writeLines(paste("\nProcessing file:", f))
        
        # Make sure file exists, but let's not be picky about the extension as bedgraphs might not have .bedgraph
        if ( ! file_checks(f) ) { next }
        
        # Define arguments to the genomeCoverageBed command
        fOut = paste(bwDest, sub(file_ext(f), "bigwig", f), sep="")
        arguments = c(f, genomeSizeFile, fOut)
        
        #################
        # Get output file
        #################
        writeLines(paste("Creating file ", fOut, "...", sep=""), sep="")
        if ( file_checks(fOut, shouldExist=FALSE) ) {
            tStart = proc.time()[3]
            system2(paste(ucscPath, "bedgraphToBigWig", sep=""), args = arguments, stdout = fOut)
            tElapsed = proc.time()[3] - tStart
            writeLines(paste("done (", round(tElapsed/60, digits=2), "m).", sep=""))
        }
        
        writeLines("Done with file.")
        
    }
    
}

