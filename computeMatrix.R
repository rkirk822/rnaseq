


# computeMatrix.R
#
# Take a list of bigwigs and use the deeptools computeMatrix function to get a matrix file for each,
# which can then be used by plotHeatmap to create an eps.
#
# I'm not yet building in flexibility with arguments to computeMatrix, but that should happen.
#
#
# USAGE:
# > computeMatrix("samplename.gencov.bigwig", regionsFile="All_mm10_wholeGenes.bed", outSuffix="genebody")
#
# TIME:
# ~14m for a 90MB bigwig, to 5MB matrix.
#
#
########################################################
# Here are the equivalent commands at the command line.
# If you want to play with the parameters while looking at
# just one file, this might be easiest.
########################################################
#
# $ computeMatrix scale-regions -S test.bw -R All_mm10_wholeGenes.bed -b 1000 -a 1000 --regionBodyLength 1000 --skipZeros -o test.mat.gz -p 2


computeMatrix = function(bigwigs, regionsFile, deeptoolsPath="", matDest="./", outSuffix="") {

    # Won't be necessary when this is in package
    source("/Volumes/CodingClub1/RNAseq/code/file_checks.R")
    source("/Volumes/CodingClub1/RNAseq/code/dir_check.R")
    library(tools)

    # Check arguments
    if (! deeptoolsPath == "") {deeptoolsPath = dir_check(deeptoolsPath)}
    matDest = dir_check(matDest)
    if ( ! file.exists(regionsFile) ) { stop("Specified regionsFile does not exist.") }
    
    # For each bigwig
    for (b in bigwigs) {
        
        writeLines(paste("\nProcessing file:", b))
        
        # Make sure file exists, but let's not be picky about the extension as bigwigs might not have .bigwig
        if ( ! file_checks(b) ) { next }
        
        # Define arguments to the computeMatrix command
        fOut = paste(matDest, sub(file_ext(b), paste(outSuffix, ".mat.gz", sep=""), b), sep="")
        arguments = c("scale-regions", "-S", b, "-R", regionsFile, "-b", "1000", "-a", "1000", "--regionBodyLength", "1000", "--skipZeros", "-o", fOut, "-p", "2")
        
        #################
        # Get output file
        #################
        writeLines(paste("Creating file ", fOut, "...", sep=""), sep="")
        errFile = sub(".gz", ".err.txt", fOut)
        if ( file_checks(fOut, shouldExist=FALSE) ) {
            tStart = proc.time()[3]
            system2(paste(deeptoolsPath, "computeMatrix", sep=""), args = arguments, stdout = fOut, stderr = errFile)
            tElapsed = proc.time()[3] - tStart
            writeLines(paste("done (", round(tElapsed/60, digits=2), "m).", sep=""))
        }
        
        writeLines(paste("Done with file.  See", errFile, "for messages displayed by computeMatrix."))
        
    }
    
}

