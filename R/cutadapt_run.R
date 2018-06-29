#' Trim reads in fastq files
#'
#' Trim reads using cutadapt.  Written using cutadapt v1.16.
#' @param readFilesIn Character - files with sequences to trim; can be gzipped (if paired-end data, these are the R1 files)
#' @param adapters List - Each element is an adapter sequence; element names are adapter names
#' @param cutadaptPath String - Path to directory with cutadapt executable
#' @param outDest String - Directory where output files should be saved
#' @param outSuffix String - will be appended to original filename (and followed by "_trimmed")
#' @param readFilesInR2 Character - The R2 files, if paired-end data
#' @param qualityCutoff Numeric - If one value, trim bases with lower quality score than this from 3' end; if two csv values, trim from 3' and 5' ends respectively.  Pre-adpater removal.
#' @param minLen Numeric - Reads shorter than this will be tossed
#' @param nAdapt Numeric - Cutadapt will assume there could be as many adapters as this on a given read
#' @param trimn Logical - Whether to trim flanking Ns (unknown bases)
#' @details Cutadapt's report, normally displayed in the terminal, goes to originalFileName_report.txt.  Keep these.  They're good if you need to quickly look back
#' at an early stage of processing, and the count_reads function reads them to get a vector of total read counts so you can quickly plot counts per sample.
#' That report includes the command line parameters, so whatever you use in this script as the adapter sequences, quality cutoff, etc, will appear there.
#'
#' TIME: 15-30m for 5-10 GB files.  So, generally a bunch of hours for a whole dataset, though small RNA is faster.
#'
#' Example at the command line (if you want to play with the parameters while looking at just one file, this might be easiest):
#' ~/.local/bin/cutadapt -a TruSeq_Index=AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG BF_RORbHTp2_1.fastq -o BF_RORbHTp2_1_trimmed.fastq --trim-n -q 20,20 -m 20 -n 3 > BF_RORbHTp2_1_report.txt
#' Example using paired-end data.  Just (1) put in -A for each adapter to trim from R2s, (2) put both input filenames, and (3) following the R1 output file name, put -p <R2_filename>.
#' adapterForward=AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
#' adapterRev=GTTCGTCTTCTGCCGTATGCTCTANNNNNNCACTGACCTCAAGTCTGCACACGAGAAGGCTAGA
#' ~/.local/bin/cutadapt -a TruSeq_Index=$adapterForward -A TruSeq_Index_Rev=$adapterRev CGTACG_S6_R1_001.fastq CGTACG_S6_R2_001.fastq -o CGTACG_S6_R1_001_trimmed.fastq -p CGTACG_S6_R1_001_trimmed.fastq --trim-n -q 20,20 -m 20 -n 3 > CGTACG_S6_report.txt
#' @examples
#' Single-end data:
#' fastqs = dir( paste(projectPath,"raw/",sep=""), pattern=".fastq" )
#' cutadapt_run(paste(projectPath, "raw/", fastqs, sep=""), outDest=paste(projectPath, "trimmed/", sep=""), qualityCutoff=c(20,20), minLen=20, nAdapt=3, trimn=TRUE)
#' Paired-end data:
#' fastqs=dir(paste(projectPath,'raw',sep=''), pattern='fastq', full.names=TRUE)
#' r1s=fastqs[which(regexpr("R1",fastqs)>0)]
#' r2s=fastqs[which(regexpr("R2",fastqs)>0)]
#' cutadapt_run(r1s, outDest=paste(projectPath, "trimmed/", sep=""), qualityCutoff=c(20,20), minLen=20, nAdapt=3, trimn=TRUE, readFilesInR2=r2s)
#' @author Emma Myers
#' @export

cutadapt_run = function(readFilesIn,
                        adapters = list("TruSeq_Universal" = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
                                        "TruSeq_Index" = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"),
                        cutadaptPath="~/.local/bin/", outDest="./", outSuffix="", qualityCutoff=c(0,0), minLen=0, nAdapt=1, trimn=FALSE,
                        readFilesInR2=NULL, adaptersRev=list("TruSeq_Index_Rev" = "GTTCGTCTTCTGCCGTATGCTCTANNNNNNCACTGACCTCAAGTCTGCACACGAGAAGGCTAGA") ) {

    # Check arguments
    if ( !all(file.exists(readFilesIn)) ) {
        writeLines("Missing input files:")
        writeLines(readFilesIn[which(!file.exists(readFilesIn))])
        stop("Missing input file(s).  See above.")
    }
    if (! cutadaptPath == "") { cutadaptPath = dir_check(path.expand(cutadaptPath)) }
    outDest = dir_check(outDest)

    # Check if paired-end data
    if ( !is.null(readFilesInR2) ) {
        if ( !all(file.exists(readFilesInR2)) ) {
            writeLines("Missing R2 files:")
            writeLines(readFilesIn[which(!file.exists(readFilesInR2))])
            stop("Missing R2 file(s).  See above.")
        }
        if ( !(length(readFilesIn)==length(readFilesInR2)) ) {
            stop("Different number of R1 and R2 files.")
        }
        # Make sure there are R1 files in readFilesIn and R2 files in readFilesInR2
        # I don't actually know if small letters are ever used but let's be safe and convert everything to upper-case
        fn1 = toupper( tools::file_path_sans_ext(basename(readFilesIn)) )
        fn2 = toupper( tools::file_path_sans_ext(basename(readFilesInR2)) )
        if ( !( all( regexpr("R1", fn1) > 0 ) && all( regexpr("R2", fn2) > 0 ) ) ) {
            stop("For paired-end data, readFilesIn must all have R1 or r1 in the file name and readFilesInR2 must all have R2 or r2.")
        }
        # Try taking out R1/R2 and make sure filenames match
        # This isn't perfect but it should catch most mismatches
        if ( !( all( gsub("R1", "", fn1) == gsub("R2", "", fn2) ) ) ) {
            stop("Mismatched filenames in R1 and R2 files (removing R1 and R2 from filenames does not yield pairs of identical filenames).  Note both files of a pair must be either zipped or unzipped.")
        }
    }

    # Turn adapter list into the form we need
    adapterArgs = ""
    for (a in 1:length(adapters) ) {
        adapterArgs = c(adapterArgs, "-a", paste(names(adapters)[a], "=", adapters[a], sep=""))
    }
    adapterArgs = adapterArgs[-1]
    # Same for adapters to be trimmed off of R2s, if paired-end
    adapterRevArgs = ""
    if ( !is.null(readFilesInR2) ) {
        for (a in 1:length(adaptersRev) ) {
            adapterRevArgs = c(adapterRevArgs, "-A", paste(names(adaptersRev)[a], "=", adaptersRev[a], sep=""))
        }
    }

    # Turn qualityCutoff into the form we need
    qualityCutoff = paste(qualityCutoff, collapse=",")

    counter = 1

    for (f in readFilesIn) {

        writeLines(paste("\nProcessing file:", f))

        # Create output filenames (actual trimmed file and "report" that cutadapt usually prints to screen)
        # if input file is gzipped, avoid treating .gz as the file extension
        if ( tools::file_ext(f) == "gz" ) { fileExt = paste(tools::file_ext( gsub(".gz", "", f) ), ".gz", sep="") } else { fileExt = tools::file_ext(f) }
        fOut = paste(outDest, gsub(paste(".", fileExt, sep=""), "", basename(f)), outSuffix, "_trimmed.", fileExt, sep="")
        fReport = paste(outDest, gsub(paste(".", fileExt, sep=""), "", basename(f)), outSuffix, "_report.txt", sep="")

        # Check if output file already exists
        if ( file_checks(fOut, shouldExist=FALSE, verbose=TRUE) ) {

            writeLines("Trimmed file will be saved to:", fOut)

            # Define arguments to the cutadapt command, starting with input file
            arguments = f
            # if paired-end data, include R2 filename
            if ( !is.null(readFilesInR2) ) { arguments = c(arguments, readFilesInR2[counter]) }
            arguments = c(arguments, "-o", fOut)
            if ( !is.null(readFilesInR2) ) {
                f2 = readFilesInR2[counter]
                if ( tools::file_ext(f2) == "gz" ) { fileExt2 = paste(tools::file_ext( gsub(".gz", "", f2) ), ".gz", sep="") } else { fileExt2 = tools::file_ext(f2) }
                fOut2 = paste(outDest, gsub(paste(".", fileExt2, sep=""), "", basename(f2)), outSuffix, "_trimmed.", fileExt2, sep="")
                arguments = c(arguments, "-p", fOut2)
            }
            arguments = c(arguments, "-q", qualityCutoff, "-m", minLen, "-n", nAdapt, adapterArgs, adapterRevArgs)
            if ( trimn ) { arguments = c(arguments, "--trim-n") }

            # Get the counts and display time taken
            writeLines("Trimming. . .", sep="")
            tStart = proc.time()[3]
            system2( paste(cutadaptPath, "cutadapt", sep=""), args = arguments, stdout = fReport )
            writeLines(paste("done (", round( (proc.time()[3] - tStart) /60, digits=2), "m).", sep=""))

        } else {  writeLines("Skipping file.") }

        counter = counter + 1

    }

}



