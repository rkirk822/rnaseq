

# plotHeatmap.R
#
# Very simple; just run through matrix files created by computeMatrix, creating an eps with histogram and heatmap for each one.
#
# I'm not yet building in flexibility with arguments to plotHeatmap, but that should happen.
#
# Note that you can't give plotHeatmap a missingDataColor value using #ffffff format when calling
# it through R.  So far as I can tell the # is getting interpreted as a comment, and it thinks you're
# not giving it anything.
#
#
# USAGE:
# >
#
# TIME:
#
#
#
########################################################
# Here are the equivalent commands at the command line.
# If you want to play with the parameters while looking at
# just one file, this might be easiest.
########################################################
# $ plotHeatmap -m samplename_genebody.mat.gz -out samplename_genebody.eps --sortRegions descend --colorMap hot_r -min 0 -max 1 --missingDataColor "#ffffff"

plotHeatmap = function(matrixFiles, deeptoolsPath="", epsDest="./", outSuffix="") {
    
    # Won't be necessary when this is in package
    source("/Volumes/CodingClub1/RNAseq/code/file_checks.R")
    source("/Volumes/CodingClub1/RNAseq/code/dir_check.R")
    library(tools)
    
    # Check arguments
    if (! deeptoolsPath == "") {deeptoolsPath = dir_check(deeptoolsPath)}
    epsDest = dir_check(epsDest)
    
    # For each matrix file
    for (m in matrixFiles) {
        
        writeLines(paste("\nProcessing file:", m))
        
        # Make sure file exists and ends in .mat.gz
        if ( ! file_checks(m, extension=".mat.gz") ) { next }
        
        # Define arguments to the plotHeatmap command
        fOut = paste(epsDest, sub(".mat.gz", paste(outSuffix, ".eps", sep=""), m), sep="")
        arguments = c("-m", m, "-out", fOut, "--sortRegions", "descend", "--colorMap", "hot_r", "-min", "0", "-max", "1", "--missingDataColor", "white")
        
        #################
        # Get output file
        #################
        writeLines(paste("Creating file ", fOut, "...", sep=""), sep="")
        if ( file_checks(fOut, shouldExist=FALSE) ) {
            tStart = proc.time()[3]
            system2(paste(deeptoolsPath, "plotHeatmap", sep=""), args = arguments, stdout = fOut)
            tElapsed = proc.time()[3] - tStart
            writeLines(paste("done (", round(tElapsed/60, digits=2), "m).", sep=""))
        }
        
        writeLines(paste("Done with file."))
        
    }
    
}
