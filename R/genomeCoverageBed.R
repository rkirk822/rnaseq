#' BAM or bed files to bedgraph with histogram of coverage values
#'
#' For each specified bed or bam file, apply the bedtools genomeCoverageBed function to get
#' a bedgraph file with histogram of coverage values.
#' @param filenames Character - List of BAM files
#' @param genomeSizeFile String - Filename (with path) to text file containing genome size information
#' @param bedtoolsPath String - Path to bedtools executables
#' @param samtoolsPath String - Path to samtools exectuables
#' @param outDest String - Directory where bedgraph files should be written
#' @param outSuffix String - will be appended to original filename
#' @param norm Logical - whether to normalize read counts to RPM (requires BAM input)
#' @param strand String - For paired-end data; only do forward- or only do reverse-strand reads ("+" or "-"): "plus" or "minus" will be appended to filename
#' @details For each specified bed or bam file, apply the bedtools genomeCoverageBed function to get
#' a bedgraph file with histogram of coverage values.
#' Requires bedtools and samtools to be installed.
#' Requires tools package.
#' Note that it's a bit faster with bed than bam input files.  But if you do "norm", you have to give it a BAM file, not a bed file.
#' TIME:  ~2.5m for 600 MB bam file; ~4m for an 815 MB bam file.  Slightly shorter starting from bed file (2 and 3m).
#' IMPROVE:
#' Make it time calculating scaling factors too.
#' And give total time for each file and for all files.
#' Command line examples:
#' $ /opt/bedtools2/bin/genomeCoverageBed -g /path/to/genome_size.tab.txt -i samplename.bed -bg -split > samplename.gencov.bedgraph
#' Or, for bam input:
#' $ /opt/bedtools2/bin/genomeCoverageBed -g /path/to/genome_size.tab.txt -ibam samplename.bam -bg -split > samplename.gencov.bedgraph
#' To normalize by reads per million:
#' $ uniqueReads=$(samtools view -F 0x904 -c samplename.bam)  # those flags in the samtools view command get the number of uniquely mapped reads in BAM
#' $ /opt/bedtools2/bin/genomeCoverageBed -g /path/to/genome_size.tab.txt -ibam samplename.bam -bg -scale 1000000/uniqueReads > samplename.gencov.bedgraph
#' @examples
#' genomeCoverageBed('samplename.filtered.bed', '/path/to/genome_size.tab.txt', outSuffix = 'unscaled', bedtoolsPath='/opt/bedtools2/bin/')
#' @author Emma Myers
#' @export


genomeCoverageBed = function(filenames, genomeSizeFile, bedtoolsPath="/opt/bedtools2/bin/", samtoolsPath="~/anaconda2/bin/",
                             outDest="./", outSuffix="", norm=FALSE, strand=NULL) {

    # Check arguments
    if (! bedtoolsPath == "") { bedtoolsPath = dir_check(bedtoolsPath) }
    if (! samtoolsPath == "") { samtoolsPath = dir_check(path.expand(samtoolsPath)) }
    outDest = dir_check(outDest)
    if ( ! file.exists(genomeSizeFile) ) { stop("Specified genomeSizeFile does not exist.") }
    if ( norm & !all(tools::file_ext(filenames)=='bam')) { stop("I can only normalize coverage if given BAM files.") }

    # If we're NOT scaling, put that in the filename
    if (!norm) { outSuffix = paste("_unscaled", outSuffix, sep="") }
    # And if we're only doing one strand of paired-end data, put that in the filename
    if (!is.null(strand)) {
        if ( strand == "+" ) { outSuffix = paste("_plus", outSuffix, sep="") }
        if ( strand == "-" ) { outSuffix = paste("_minus", outSuffix, sep="") }
    }
    for (f in filenames) {

        writeLines("\n")
        # Make sure file exists.  Should be checking that it's either .bam or .bed too.
        if ( ! file_checks(f, verbose=TRUE) ) { next }

        writeLines(paste("Processing file:", f))

        # Check if output file already exists
        fOut = paste(outDest, sub(paste(".", tools::file_ext(f), sep=""), paste(outSuffix, "_gencov.bedgraph", sep=""), basename(f)), sep="")
        if ( file_checks(fOut, shouldExist=FALSE, verbose=TRUE) ) {

            # If we're normalizing coverage values, calculate scaling factor
            scalingFactor = 1
            if (norm) {
                writeLines("Calculating scaling factor for normalization...", sep="")
                readCount = count_alignments(f, samtoolsPath, what="primary")
                writeLines(paste(readCount, "uniquely mapped reads in sample..."), sep="")
                scalingFactor = 10^6 / readCount
                writeLines("done.")
            }

            # Define arguments to the genomeCoverageBed command
            arguments = c(f, "-g", genomeSizeFile, "-bg", "-scale", scalingFactor, "-split")
            if ( !is.null(strand) ) { arguments = c(arguments, "-strand", strand) }
            if (tools::file_ext(f) == "bam") { arguments = c("-ibam", arguments) } else { arguments = c("-i", arguments) }

            # Get output file
            writeLines(paste("Creating file ", fOut, "...", sep=""), sep="")
            tStart = proc.time()[3]
            system2(paste(bedtoolsPath, "genomeCoverageBed", sep=""), args = arguments, stdout = fOut)
            tElapsed = proc.time()[3] - tStart
            writeLines(paste("done (", round(tElapsed/60, digits=2), "m).", sep=""))
            writeLines("Done with file.")
        }

    }

}

