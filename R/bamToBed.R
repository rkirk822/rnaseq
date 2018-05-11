#' Convert BAM files to bed format
#'
#' Take a list of BAM files and write a bed file for each one.
#' @param filenames Character - List of BAM files
#' @param bedtoolsPath String - Path to bedtools directory
#' @param outDest String - Directory where bed files should be written
#' @param outSuffix String - will be appended to original filename (and ".bam" is replaced by ".bed")
#' @details Use bedtools function bamToBed convert BAM files to bed format.  Requires bedtools to be installed.
#' Mitochondrial, unknown, and "random" reads should have been removed already (see sam_filter.R).
#' If you give the output file a suffix with "outSuffix", you need to start it with a "_" or "." or
#' whatever if you don't want it just smooshed onto the end of the input file name.
#' Time:  For NucSeq data, ranged from ~2.5m for cb_p036_1 to ~14.5m for L4_RORb_1.
#' Example at the command line:
#' /opt/bedtools2/bin/bamToBed -i samplename.bam > samplename.bed
#' @examples
#' samplenames = read.csv('samplenames.txt', comment.char='#')
#' filenames = paste(samplenames[,1], '_mapped_Aligned.out.bam', sep='')
#' bamToBed(filenames, bedtoolsPath = '/opt/bedtools2/bin/')
#' @author Emma Myers
#' @export

bamToBed = function(filenames, bedtoolsPath = "", outDest = "./", outSuffix="" ) {

    # Check arguments
    if (! bedtoolsPath == "") { bedtoolsPath = dir_check(bedtoolsPath) }
    outDest = dir_check(outDest)

    for (f in filenames) {

        writeLines("\n")
        # Make sure file exists and has extension .bam
        if ( ! file_checks(f, extension = ".bam", verbose=TRUE) ) { next }

        writeLines(paste("Processing file:", f))

        #################
        # Get bed file
        #################
        fBed = paste(outDest, sub(".bam", paste(outSuffix, ".bed", sep=""), basename(f)), sep="")
        writeLines(fBed)
        if ( file_checks(fBed, shouldExist=FALSE, verbose=TRUE) ) {
            writeLines(paste("Creating file ", fBed, "...", sep=""), sep="")
            tStart = proc.time()[3]
            system2(paste(bedtoolsPath, "bamToBed", sep=""), args = c("-i", f, "-split"), stdout = fBed)
            tElapsed = proc.time()[3] - tStart
            writeLines(paste("done (", round(tElapsed/60, digits=2), "m).", sep=""))
            writeLines("Done with file.")
        }

    }

}

