

# plotHeatmap.R
#
# Very simple; just run through matrix files created by computeMatrix, creating an eps with histogram and heatmap for each one.
#
# IMPROVE: If given sample and region labels, some kind of check to see that they correspond to what's
# in the matfile somehow.  And check the length of the character vectors to make sure it matches the
# length of matrixFiles.  Also you have to be in the directory where the input files are right now.
# 
# Note that you can't give plotHeatmap a missingDataColor value using #ffffff format when calling
# it through R.  So far as I can tell the # is getting interpreted as a comment, and it thinks you're
# not giving it anything.
# 
# The latest version of deeptools by default uses interpolation method "bilinear", so I've set that to
# be the default here.  However, it fails for our data (at least the TTX and the small RNA), whereas
# "nearest" works.  
#
# Though it's told the color for missing data is white, it seems to be doing black anyway.  I believe this was
# also true when doing it directly in the command line, in which case it's not about going through R.
# But check that.
#
# outMatrix, if true, tells it to output the values in the matrix into an additional output file.
# I think this might be mainly useful when you have a bunch of input files and this is combining them?
#
# USAGE:
# > mats = list.files(pattern='.mat.gz')
# > plotHeatmap(mats, samplesLabels=sub('.*_', '', sub('_Aligned.*', '', matfiles)))
#
# TIME:
# ~10s per matrix file
#
#
########################################################
# Here are the equivalent commands at the command line.
# If you want to play with the parameters while looking at
# just one file, this might be easiest.
########################################################
# $ plotHeatmap -m samplename_genebody.mat.gz -out samplename_genebody.eps --sortRegions descend --colorMap hot_r -min 0 -max 1 --missingDataColor "#ffffff"

plotHeatmap = function(matrixFiles, deeptoolsPath="", outDest="./", outSuffix="", outMatrix=FALSE, 
                       samplesLabels=NULL, regionsLabels=NULL, minVal=NULL, maxVal=NULL, interpolationMethod=bilinear) {
    
    # Won't be necessary when this is in package
    # source("/Volumes/CodingClub1/RNAseq/code/file_checks.R")
    # source("/Volumes/CodingClub1/RNAseq/code/dir_check.R")
    # source("/Users/nelsonlab/Documents/Toolboxes/rna-seq/file_checks.R")
    # source("/Users/nelsonlab/Documents/Toolboxes/rna-seq/dir_check.R")
    source("/Users/work/Documents/rna-seq/file_checks.R")
    source("/Users/work/Documents/rna-seq/dir_check.R")
    library(tools) # for I forget what
    library(glue)  # for collapse()
    
    # Check arguments
    if (! deeptoolsPath == "") {deeptoolsPath = dir_check(deeptoolsPath)}
    outDest = dir_check(outDest)
    
    counter = 1
    
    # For each matrix file
    for (m in matrixFiles) {
        
        writeLines("\n")
        
        # Make sure file exists and ends in .mat.gz
        if ( ! file_checks(m, extension=".mat.gz", verbose=TRUE) ) { next }
        
        writeLines(paste("Processing file:", m))
        
        # Define arguments to the plotHeatmap command
        fOut = paste(outDest, sub(".mat.gz", paste(outSuffix, ".eps", sep=""), m), sep="")
        arguments = c("-m", m, "-out", fOut, "--sortRegions", "descend", "--colorMap", "hot_r", 
                      "--missingDataColor", "white", 
                      "--interpolationMethod", interpolationMethod)
        if (outMatrix) {
            arguments = c(arguments, "--outFileNameMatrix", sub("mat.gz", "txt", m))
        }
        
        if ( !is.null(samplesLabels) ) {
            arguments = c(arguments, "--samplesLabel", samplesLabels[counter])
        }
        
        if ( !is.null(regionsLabels) ) {
            arguments = c(arguments, "--regionsLabel", collapse(paste("\'", regionsLabels, "\'", sep=""), sep=" "))
        }
        
        if ( !is.null(minVal) ) {
            arguments = c(arguments, "-min", minVal)
        }
        
        if ( !is.null(maxVal) ) {
            arguments = c(arguments, "-max", maxVal)
        }
        
        #################
        # Get output file
        #################
        if ( file_checks(fOut, shouldExist=FALSE, verbose=TRUE) ) {
            writeLines(paste("Creating file ", fOut, "...", sep=""), sep="")
            tStart = proc.time()[3]
            system2(paste(deeptoolsPath, "plotHeatmap", sep=""), args = arguments, stdout = fOut)
            tElapsed = proc.time()[3] - tStart
            writeLines(paste("done (", round(tElapsed/60, digits=2), "m).", sep=""))
            writeLines(paste("Done with file."))
        }
        
        counter = counter + 1
        
    }
    
}
