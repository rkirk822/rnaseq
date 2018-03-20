

# UNFINISHED
# Make it time calculating scaling factors too.
# And give total time for each file and for all files.
# But it's working, as is.

# genomeCoverageBed.R
#
# For each specified bed or bam file, apply the bedtools genomeCoverageBed function to get
# a bedgraph file with histogram of coverage values.  Right now the only flexibility I'm
# giving you is to choose bam or bed as your input, but I can incorporate other arguments to
# genomeCoverageBed.
#
# Note that it's a bit faster with bed than bam input files.
#
# Mitochondrial, unknown, and "random" reads should have been removed already (see sam_filter.R).
#
# If you do "norm", samtools has to be installed, and you have to give it a BAM file, not a bed file.
#
#
# USAGE:
# > genomeCoverageBed('samplename_mapped_Aligned.out.filtered.bed', '/path/to/genome_size.tab.txt', outSuffix = 'unscaled', bedtoolsPath='/opt/bedtools2/bin/')
#
# TIME:
# On yasuimac: 2.5m for 600 MB bam file
#              4m for an 815 MB bam file
#              Slightly shorter starting from bed file (2 and 3m).
#
########################################################
# Here are the equivalent commands at the command line:
########################################################
# $ /opt/bedtools2/bin/genomeCoverageBed -g /path/to/genome_size.tab.txt -i samplename.bed -bg > samplename.gencov.bedgraph
# Or, for bam input:
# $ /opt/bedtools2/bin/genomeCoverageBed -g /path/to/genome_size.tab.txt -ibam samplename_mapped_Aligned_out.bam -bg > samplename.gencov.bedgraph
#
# To normalize by reads per million:
# $ uniqueReads=$(samtools view -F 0x904 -c samplename_mapped_Aligned_out.bam)
# $ /opt/bedtools2/bin/genomeCoverageBed -g /path/to/genome_size.tab.txt -ibam samplename_mapped_Aligned_out.bam -bg -scale 1000000/uniqueReads > samplename.gencov.bedgraph

genomeCoverageBed = function(filenames, genomeSizeFile, bedtoolsPath="", bedgraphDest="./", outSuffix="", norm=FALSE, samtoolsPath="") {

    # Won't be necessary when this is in package
    library(tools)
    source("/Volumes/CodingClub1/RNAseq/code/file_checks.R")
    source("/Volumes/CodingClub1/RNAseq/code/dir_check.R")
    source("/Volumes/CodingClub1/RNAseq/code/alignment_counts.R")

    # Check arguments
    if (! bedtoolsPath == "") { bedtoolsPath = dir_check(bedtoolsPath) }
    bedgraphDest = dir_check(bedgraphDest)
    if ( ! file.exists(genomeSizeFile) ) { stop("Specified genomeSizeFile does not exist.") }
    if ( norm & !all(file_ext(filenames)=='bam')) { stop("I can only normalize coverage if given BAM files.") }
    
    
    for (f in filenames) {
        
        writeLines("\n")
        # Make sure file exists.  Should be checking that it's either .bam or .bed too.
        if ( ! file_checks(f, verbose=TRUE) ) { next }
        
        writeLines(paste("Processing file:", f))
        
        # Check if output file already exists
        fOut = paste(bedgraphDest, sub(paste(".", file_ext(f), sep=""), paste(outSuffix, ".gencov.bedgraph", sep=""), f), sep="")
        if ( file_checks(fOut, shouldExist=FALSE, verbose=TRUE) ) {
            
            # If we're normalizing coverage values, calculate scaling factor
            if (norm) {
                writeLines("Calculating scaling factor for normalization...", sep="")
                readCount = alignment_counts(f, samtoolsPath)
                writeLines(paste(readCount, "reads in sample..."), sep="")
                scalingFactor = 10^6 / readCount
                writeLines("done.")
            }
            
            # Define arguments to the genomeCoverageBed command
            arguments = c(f, "-g", genomeSizeFile, "-bg", "-scale", scalingFactor, "-split")
            if (file_ext(f) == "bam") { arguments = c("-ibam", arguments) } else { arguments = c("-i", arguments) }
            
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

