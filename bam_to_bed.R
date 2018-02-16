

# bam_to_bed.R
#
# Take bam file and get bed file.
# Mitochondrial, unknown, and "random" reads should have been removed already (see sam_filter.R).
#
#
# USAGE:
# > samplenames = read.csv('samplenames.txt', comment.char='#')
# > filenames = paste(samplenames[,1], '_mapped_Aligned.out.bam', sep='')
# > bam_to_bed(filenames, bedtoolsPath = '/opt/bedtools2/bin/')
#
# TIME:
# For NucSeq data, ranged from ~2.5m for cb_p036_1 to ~14.5m for L4_RORb_1.
#
#
########################################################
# Here are the equivalent commands at the command line:
########################################################
#
# $ /opt/bedtools2/bin/bamToBed -i samplename_mapped_Aligned.out.bam > samplename.bed


bam_to_bed = function(filenames, bedtoolsPath = "", bedDest = "./" ) {

    # Won't be necessary when this is in package
    source("/Volumes/CodingClub1/RNAseq/code/file_checks.R")
    source("/Volumes/CodingClub1/RNAseq/code/dir_check.R")

    # Check arguments
    if (! bedtoolsPath == "") { bedtoolsPath = dir_check(bedtoolsPath) }
    bedDest = dir_check(bedDest)
    
    for (f in filenames) {
        
        writeLines(paste("\nProcessing file:", f))
    
        # Make sure file exists and has extension .bam
        if ( ! file_checks(f, extension = ".bam") ) { next }
        
        #################
        # Get bed file
        #################
        fBed = paste(bedDest, sub(".bam", ".bed", f), sep="")
        writeLines(paste("Creating file ", fBed, "...", sep=""), sep="")
        if ( file_checks(fBed, shouldExist=FALSE) ) {
            tStart = proc.time()[3]
            system2(paste(bedtoolsPath, "bamToBed", sep=""), args = c("-i", f), stdout = fBed)
            tElapsed = proc.time()[3] - tStart
            writeLines(paste("done (", round(tElapsed/60, digits=2), "m).", sep=""))
        }
        writeLines("Done with file.")

    }

}

