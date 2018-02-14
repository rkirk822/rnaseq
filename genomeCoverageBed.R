

# genomeCoverageBed.R
#
# For each specified bed or bam file, apply the bedtools genomeCoverageBed function to get
# a bedgraph file with histogram of coverage values.  Right now the only flexibility I'm
# giving you is choose bam or bed as your input, but I can incorporate other arguments to
# genomeCoverageBed.
#
# Note that it's a bit faster with bed than bam input files.
#
# Mitochondrial, unknown, and "random" reads should have been removed already (see sam_filter.R).
#
#
# USAGE:
# > genomeCoverageBed('samplename_mapped_Aligned.out.filtered.bed', bedtoolsPath='/opt/bedtools2/bin/', genomeSizeFile = '/path/to/genome_size.tab.txt')
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


genomeCoverageBed = function(filenames, genomeSizeFile, bedtoolsPath="", bedgraphDest='./') {

    # Won't be necessary when this is in package
    source("/Volumes/CodingClub1/RNAseq/code/file_checks.R")
    source("/Volumes/CodingClub1/RNAseq/code/dir_check.R")

    # Check arguments
    if (! bedtoolsPath == "") { bedtoolsPath = dir_check(bedtoolsPath) }
    bedgraphDest = dir_check(bedgraphDest)
    if ( ! file.exists(genomeSizeFile) ) { stop("Specified genomeSizeFile does not exist.") }
    
    
    for (f in filenames) {
        
        writeLines(paste("\nProcessing file:", f))
        
        # Make sure file exists and has extension .bam or .bed
        # Note the 4 is because both of those are length 4
        inFormat = substr(f, (nchar(f)-4)+1, nchar(f))
        if ( ! file_checks(f, extension = inFormat) ) { next }
        
        # Define arguments to the genomeCoverageBed command
        arguments = c(f, "-g", genomeSizeFile, "-bg")
        if (inFormat == ".bam") { arguments = c("-ibam", arguments) } else { arguments = c("-i", arguments) }
        
        #################
        # Get output file
        #################
        fOut = paste(bedgraphDest, sub(inFormat, ".gencov.bedgraph", f), sep="")
        writeLines(paste("Creating file ", fOut, "...", sep=""), sep="")
        if ( file_checks(fOut, shouldExist=FALSE) ) {
            tStart = proc.time()[3]
            system2(paste(bedtoolsPath, "genomeCoverageBed", sep=""), args = arguments, stdout = fOut)
            tElapsed = proc.time()[3] - tStart
            writeLines(paste("done (", round(tElapsed/60, digits=2), "m).", sep=""))
        }
        
        writeLines("Done with file.")
        
    }

}

