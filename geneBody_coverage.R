

# geneBody_coverage.R
#
#
#
#
# USAGE:
# > bamfiles = list.files(pattern = '.bam')
# > bamfiles = bamfiles[which(regexpr('.bai', bamfiles) < 0)]
# > geneBody_coverage(bamfiles, regionsFile='/path/to/bedfile/All_mm10_wholeGenes.bed', outPrefixMarker='_mapped', dispMessages=FALSE)
#
# TIME:
#
#
#
########################################################
# Here are the equivalent commands at the command line:
########################################################
#
# $ geneBody_coverage.py -r /path/to/bedfile/All_mm10_wholeGenes.bed -i /path/to/bamfile/samplename_mapped_Aligned.out.filtered.bam -o /path/to/destination/samplename
#
# For geneBody_coverage2.py, replace the bam file with a bigwig.


geneBody_coverage = function(filenames, regionsFile, outDest = "./", outPrefixMarker = NA, dispMessages = TRUE) {

    # Won't be necessary when this is in package
    # source("/Volumes/CodingClub1/RNAseq/code/file_checks.R")
    # source("/Volumes/CodingClub1/RNAseq/code/dir_check.R")
    source("/Users/nelsonlab/Documents/Toolboxes/rna-seq/file_checks.R")
    source("/Users/nelsonlab/Documents/Toolboxes/rna-seq/dir_check.R")
    library(tools)

    # Check arguments
    outDest = dir_check(outDest)
    
    for (f in filenames) {
        
        writeLines("\n")
        # Make sure file exists and has extension .bam or .bigwig
        inFormat = file_ext(f)
        if ( ! file_checks(f, extension = inFormat, verbose=TRUE) ) { next }
        
        writeLines(paste("Processing file:", f))
        
        # Are we calling geneBody_coverage or geneBody_coverage2
        if (inFormat == "bam") {
            funcName = "geneBody_coverage.py"
        } else if (inFormat == "bigwig") {
            funcName = "geneBody_coverage2.py"
        } else {
            stop("Input file must be either bam or bigwig format.")
        }
        
        # Define arguments to the geneBody_coverage (or geneBody_coverage2) command
        outPrefix = file_path_sans_ext(f)
        if ( ! is.na(outPrefixMarker) ) {
            outPrefix = strsplit(outPrefix, outPrefixMarker)[[1]][1]
        }
        arguments = c("-r", regionsFile, "-i", f, "-o", paste(outDest, outPrefix, sep="") )
        
        ##################
        # Process the file
        ##################
        outText = paste(outDest, outPrefix, ".geneBodyCoverage.txt", sep="")
        outR = paste(outDest, outPrefix, ".geneBodyCoverage.r", sep="")
        if (funcName == "geneBody_coverage.py") {
            outPdf = paste(outDest, outPrefix, ".geneBodyCoverage.curves.pdf")
        } else { outPdf = paste(outDest, outPrefix, ".geneBodyCoverage.pdf") }
        outErr = paste(outDest, outPrefix, ".err.txt", sep="")
        if ( file_checks(outText, shouldExist=FALSE, verbose=TRUE) && file_checks(outR, shouldExist=FALSE, verbose=TRUE) && file_checks(outPdf, shouldExist=FALSE, verbose=TRUE) ) {
            writeLines(paste("Creating output files with prefix ", paste(outDest, outPrefix, "...", sep=""), sep=""))
            tStart = proc.time()[3]
            if (dispMessages) {
                system2(funcName, args = arguments)
            } else {
                system2(funcName, args = arguments, stderr = outErr )
            }
            tElapsed = proc.time()[3] - tStart
            writeLines(paste("done (", round(tElapsed/60, digits=2), "m).", sep=""))
            writeLines("Done with file.", sep="")
            if ( ! dispMessages) {
                writeLines(paste("See ", outErr, " for messages displayed by ", funcName, ".", sep=""))
            } else { writeLines("\n") }
        }

    }

}

