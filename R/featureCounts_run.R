#' Read counts from BAM or SAM files
#'
#' Get raw counts of reads mapped to specified genomic features using featureCount tool from the subread package.
#' @param inFiles Character - BAM or SAM file list
#' @param annotFile String - Name (with path) of gtf file with annotations to define features
#' @param fcPath String - Path to directory with featureCounts executable
#' @param outDest String - Directory where output files should be saved
#' @param outSuffix String - will be appended to original filename (and followed by "_fcounts.txt")
#' @param runThreadN Numeric - How many cores to use
#' @param pairedEnd Logical - Whether data is paired-end
#' @param isStrandSpecific Logical - Strand-specific counting (reads will only be counted if they match the strand of the feature)
#' @param annotFormat String - "GTF" or "SAF".  Format of annotation file.
#' @param multimappers Logical - Whether to count multi-mapping reads.  All reported alignments will be counted, using the 'NH' tag in the BAMs/SAMs.
#' @param dispToText Logical - Whether to send messages that featureCounts normally displays to the screen, to a text file instead
#' @details Take a list of SAM or BAM files, and get counts of reads mapped to genomic features in specified annotations file.
#' Command issued will be written to a text file, particularly so you can confirm the annotation file used (unless its filename gets changed).
#' TIME:  1-3m per BAM.  Paired-end data is significantly longer.  3-9m on small RNA data when features are ATAC peaks, so probably longer for larger files when features are all exons or all introns
#' Example at the command line (if you want to play with the parameters while looking at just one file, this might be easiest):
#' annotFile=/Volumes/CodingClub1/RNAseq/Metadata/mm10_refGene.gtf # for intronic read counts, use mm10_refSeq_introns_geneids.gtf
#' featureCounts -a $annotFile -o $outputFile $inputFile -T 8
#' @author Emma Myers
#' @export

featureCounts_run = function(inFiles, annotFile, fcPath="/opt/subread-1.6.0-MacOSX-x86_64/bin/", outDest="./", outSuffix="", runThreadN=1,
                             pairedEnd=FALSE, isStrandSpecific=FALSE, annotFormat="GTF", multimappers=FALSE, dispToText=FALSE) {

    # Check arguments
    if ( !all(file.exists(inFiles)) ) {
        writeLines("Missing input files:")
        writeLines(inFiles[which(!file.exists(inFiles))])
        stop("Missing input file(s).  See above.")
    }
    if ( !file.exists(annotFile) ) { stop("Specified annotations file does not exist.") }
    if ( !(fcPath == "") ) { fcPath = dir_check(fcPath) }
    outDest = dir_check(outDest)

    for (f in inFiles) {

        writeLines(paste("\nProcessing file:", f))

        fOut = paste(outDest, gsub(".bam", "", basename(f)), outSuffix, "_fcounts.txt", sep="")
        # Check if there's already an output file with this name
        if ( file_checks(fOut, shouldExist=FALSE, verbose=TRUE) ) {

            # Define arguments to the featureCounts command
            arguments = c(f, "-a", annotFile, "-o", fOut, "-T", runThreadN, "-F", annotFormat)
            if ( multimappers ) { arguments = c(arguments, "-M") }
            if ( pairedEnd ) { arguments = c(arguments, "-p") }
            if ( isStrandSpecific ) { arguments = c(arguments, "-s") }

            # Get the counts
            tStart = proc.time()[3]
            if (dispToText) {
                system2( paste(fcPath, "featureCounts", sep=""), args = arguments, stderr = gsub("_fcounts.txt", "_display.txt", fOut) )
            } else { system2( paste(fcPath, "featureCounts", sep=""), args = arguments) }

            # Save information about command to a textfile
            fRecordCon = file(gsub("_fcounts.txt", "_record.txt", fOut), "w")
            writeLines( paste( paste(fcPath, "featureCounts", sep=""), paste(arguments, collapse=" ") ), con = fRecordCon )
            close(fRecordCon)

            tElapsed = proc.time()[3] - tStart
            writeLines(paste("done (", round(tElapsed/60, digits=2), "m).", sep=""))

        }

    }

}

