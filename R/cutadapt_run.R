#' Trim reads in fastq files
#'
#' Trim reads using cutadapt.  Written using cutadapt v1.16.
#' @param readFilesIn Character - Fastq file list (can be gzipped)
#' @param adapters List - Each element is an adapter sequence; element names are adapter names
#' @param cutadaptPath String - Path to directory with cutadapt executable
#' @param outDest String - Directory where output files should be saved
#' @param outSuffix String - will be appended to original filename (and followed by "_trimmed")
#' @param qualityCutoff Numeric - If one value, trim bases with lower quality score than this from 3' end; if two csv values, trim from 3' and 5' ends respectively.  Pre-adpater removal.
#' @param minLen Numeric - Reads shorter than this will be tossed
#' @param nAdapt Numeric - Cutadapt will assume there could be as many adapters as this on a given read
#' @param trimn Logical - Whether to trim flanking Ns (unknown bases)
#' @details Cutadapt's report, normally displayed in the terminal, goes to originalFileName_report.txt.  Keep these.  They're good if you need to quickly look back
#' at an early stage of processing, and the count_reads function reads them to get a vector of total read counts so you can quickly plot counts per sample.
#' That report includes the command line parameters, so whatever you use in this script as the adapter sequences, quality cutoff, etc, will appear there.
#' TIME:
#' For 5-6 GB fastq files, about 15m.
#' For 9-10 GB, closer to 30m.
#' For all the RORb files, total ~7.5hrs.
#' Example at the command line (if you want to play with the parameters while looking at just one file, this might be easiest):
#' TESTED THIS
#' ~/.local/bin/cutadapt -a TruSeq_Index=AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG BF_RORbHTp2_1.fastq -o BF_RORbHTp2_1_trimmed.fastq --trim-n -q 20,20 -m 20 -n 3 > BF_RORbHTp2_1_report.txt
#' @examples
#' TESTED THIS
#' fqs = list.files( paste(projectPath,"raw/",sep=""), pattern=".fastq" )
#' cutadapt_run(paste(projectPath, "raw/", fqs, sep=""), outDest=paste(projectPath, "trimmed/", sep=""), qualityCutoff=c(20,20), minLen=20, nAdapt=3, trimn=TRUE)
#' @author Emma Myers
#' @export

cutadapt_run = function(readFilesIn,
                        adapters = list("TruSeq_Universal" = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
                                        "TruSeq_Index" = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"),
                        cutadaptPath="~/.local/bin/", outDest="./", outSuffix="",
                        qualityCutoff=c(0,0), minLen=0, nAdapt=1, trimn=FALSE) {

    # Check arguments
    if ( !all(file.exists(readFilesIn)) ) {
        writeLines("Missing input files:")
        writeLines(readFilesIn[which(!file.exists(readFilesIn))])
        stop("Missing input file(s).  See above.")
    }
    if (! cutadaptPath == "") { cutadaptPath = dir_check(path.expand(cutadaptPath)) }
    outDest = dir_check(outDest)

    # Turn adapter list into the form we need
    adapterArgs = ""
    for (a in 1:length(adapters) ) {
        adapterArgs = c(adapterArgs, "-a", paste(names(adapters)[a], "=", adapters[a], sep=""))
    }
    adapterArgs = adapterArgs[-1]

    # Turn qualityCutoff into the form we need
    qualityCutoff = paste(qualityCutoff, collapse=",")

    for (f in readFilesIn) {

        writeLines(paste("\nProcessing file:", f))

        # Create output filenames
        fOut = paste(outDest, gsub(paste(".", tools::file_ext(f), sep=""), "", basename(f)), outSuffix, "_trimmed.", tools::file_ext(f), sep="")
        fReport = paste(outDest, gsub(paste(".", tools::file_ext(f), sep=""), "", basename(f)), outSuffix, "_report.txt", sep="")

        # Check if output file already exists
        if ( file_checks(fOut, shouldExist=FALSE, verbose=TRUE) ) {

            writeLines("Trimmed file will be saved to:", fOut)

            # Define arguments to the cutadapt command
            arguments = c(f, "-o", fOut, "-q", qualityCutoff, "-m", minLen, "-n", nAdapt, adapterArgs)
            if ( trimn ) { arguments = c(arguments, "--trim-n") }

            # Get the counts and display time taken
            writeLines("Trimming. . .", sep="")
            tStart = proc.time()[3]
            system2( paste(cutadaptPath, "cutadapt", sep=""), args = arguments, stdout = fReport )
            writeLines(paste("done (", round( (proc.time()[3] - tStart) /60, digits=2), "m).", sep=""))

        } else {  writeLines("Skipping file.") }

    }

}



