#' Read density per region
#'
#' Write matrix containing read density per region from bigwig files
#' @param inFiles Character - List of BAM files
#' @param deeptoolsPath String - Path to deeptools directory
#' @param genomeSizeFile String - Filename (with path) to text file containing genome size information
#' @param outDest String - Directory where matrix files should be written
#' @param outSuffix String - will be appended to original filename (and ".bedgraph" is replaced by ".bigwig")
#' @details Take a list of bigwigs and use the deeptools computeMatrix function to get a matrix file for each,
#' which can then be used by plotHeatmap to create an eps.
#' TIME:  Varies a lot, maybe based on how many cores you're using?  15-35m per bigwig is typical.
#' COMMAND LINE EXAMPLE (if you want to play with the parameters while looking at just one file, this might be easiest):
#' computeMatrix scale-regions -S samplename.bw -R All_mm10_wholeGenes.bed -b 1000 -a 1000 --skipZeros -o samplename_regions.mat.gz -p 2
#' IMPROVE:
#' Have it write the deeptools command to the top of the file where text printed to screen gets sent.
#' Make it so you can do reference-point instead of scale-regions.
#' Could be two different functions if the options are really different.  a and b have different defaults.
#' Don't bother making it possible to group bigwigs.  They're still separate plots, just in the same output file.
#' @examples
#' computeMatrix("samplename.bigwig", regionsFiles=c("All_mm10_wholeGenes.bed", "Shuffled_RORht_mm10.bed"), outSuffix="genebody")
#' @author Emma Myers
#' @export

computeMatrix = function(inFiles, regionsFiles, deeptoolsPath="", outDest="./", outSuffix="", textMatrix=FALSE, startLabel=NULL, endLabel=NULL,
                         skipZeros=FALSE, missingDataAsZero=FALSE, maxThreshold=NULL, minThreshold=NULL,
                         upstream=0, downstream=0, regionBodyLength=1000, binSize=10,
                         nProcessors=1) {

    # Check arguments
    if ( !all(file.exists(inFiles)) ) {
        writeLines("Missing input files:")
        writeLines(inFiles[which(!file.exists(inFiles))])
        stop("Missing input file(s).  See above.")
    }
    if (! deeptoolsPath == "") {deeptoolsPath = dir_check(deeptoolsPath)}
    outDest = dir_check(outDest)
    if ( ! all(file.exists(regionsFiles)) ) { stop("One or more of the specified regions files do not exist.") }

    # For each bigwig
    for (f in inFiles) {

        writeLines(paste("\nProcessing file:", f))

        # Define arguments to the computeMatrix command
        fOut = paste(outDest, sub(file_ext(f), paste(outSuffix, "mat.gz", sep=""), basename(f)), sep="")
        arguments = c("scale-regions", "-S", b, "-R", paste(regionsFiles, collapse=" "),
                      "-b", upstream, "-a", downstream,
                      "--regionBodyLength", regionBodyLength, "--binSize", binSize,
                      "-o", fOut, "-p", nProcessors)
        if (skipZeros) { arguments = c(arguments, "--skipZeros") }
        if (missingDataAsZero) { arguments = c(arguments, "--missingDataAsZero") }
        if ( !is.null(maxThreshold) ) { arguments = c(arguments, "--maxThreshold", maxThreshold) }
        if ( !is.null(minThreshold) ) { arguments = c(arguments, "--minThreshold", minThreshold) }
        if (textMatrix) { arguments = c(arguments, "--outFileNameMatrix", sub(".gz", ".txt", fOut)) }
        if ( !is.null(startLabel) ) { arguments = c(arguments, "--startLabel", paste("'", startLabel, "'", sep="")) }
        if ( !is.null(endLabel) ) { arguments = c(arguments, "--endLabel", paste("'", endLabel, "'", sep="")) }

        writeLines( paste(arguments, collapse=" ") )

        #################
        # Get output file
        #################
        fDisp = sub(".gz", ".err.txt", fOut)
        if ( file_checks(fOut, shouldExist=FALSE, verbose=TRUE) ) {
            writeLines(paste("Creating file ", fOut, "...", sep=""), sep="")
            tStart = proc.time()[3]
            system2(paste(deeptoolsPath, "computeMatrix", sep=""), args = arguments, stdout = fOut, stderr = fDisp)
            tElapsed = proc.time()[3] - tStart
            writeLines(paste("done (", round(tElapsed/60, digits=2), "m).", sep=""))
            writeLines(paste("Done with file.  See", fDist, "for messages displayed by computeMatrix."))
        }

    }

}

