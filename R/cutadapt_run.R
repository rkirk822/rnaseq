#' Trim reads in fastq files
#'
#' Trim reads using cutadapt.
#' @param readFilesIn Character - Fastq file list (can be gzipped)
#' @param adapters List - Each element is an adapter sequence; element names are adapter names
#' @param cutadaptPath String - Path to directory with cutadapt executable
#' @param outDest String - Directory where output files should be saved
#' @param outSuffix String - will be appended to original filename (and followed by "_trimmed")
#' @param qualityCutoff Numeric -
#' @param minLen Numeric - Reads shorter than this will be tossed
#' @param nAdapt Numeric - Cutadapt will assume there could be as many adapters as this on a given read
#' @param trimn Logical - Whether to trim flanking Ns (unknown bases)
#' @details Cutadapt's report, which by default is displayed in the terminal, goes to originalFileName_report.txt.
#' That report includes the command line parameters, so whatever you use in this script as the adapter
#' sequences, quality cutoff, etc, will appear there.
#' Requires cutadapt to be installed.
#' TIME:
#' For 5-6 GB fastq files, about 15m.
#' For 9-10 GB, closer to 30m.
#' For all the RORb files, total ~7.5hrs.
#' Example at the command line (if you want to play with the parameters while looking at just one file, this might be easiest):
#' cutadapt -a $adapterName=$adapter -a $adapterName2=$adapter2 $infile -o $outfile --trim-n -q $qualityCutoff -m $minLen -n $nAdapt > $reportfile
#' @examples
#' get one in here
#' @author Emma Myers
#' @export

cutadapt_run = function(readFilesIn,
                        adapters = list("TruSeq_Universal" = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
                                        "TruSeq_Index" = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"),
                        cutadaptPath="~/.local/bin/", outDest="./", outSuffix="",
                        qualityCutoff=0, minLen=0, nAdapt=1, trimn=FALSE) {

    # Check arguments
    if (! cutadaptPath == "") { cutadaptPath = dir_check(path.expand(cutadaptPath)) }
    outDest = dir_check(outDest)

    # Turn adapter list into the form we need
    adapterArgs = ""
    for (a in 1:length(adapters) ) {
        adapterArgs = c(adapterArgs, "-a", paste(names(adapters)[a], "=", adapters[a], sep=""))
    }
    adapterArgs = adapterArgs[-1]

    for (f in readFilesIn) {

        writeLines("\n")
        # Make sure file exists
        if ( ! file_checks(f, verbose=TRUE) ) { next }

        writeLines(paste("Processing file:", f))

        # If input file includes path, strip it before creating output filename

        fOut = paste(outDest, gsub(paste(".", tools::file_ext(f), sep=""), "", basename(f)), outSuffix, "_trimmed.", tools::file_ext(f), sep="")
        fReport = paste(outDest, gsub(paste(".", tools::file_ext(f), sep=""), "", basename(f)), outSuffix, "_report.txt", sep="")

        # Check if output file already exists
        if ( file_checks(fOut, shouldExist=FALSE, verbose=TRUE) ) {

            # Define arguments to the cutadapt command
            arguments = c(f, "-o", fOut, "-q", qualityCutoff, "-m", minLen, "-n", nAdapt, adapterArgs)
            if ( trimn ) { arguments = c(arguments, "--trim-n") }

            # Get the counts
            system2( paste(cutadaptPath, "cutadapt", sep=""), args = arguments, stdout = fReport )

        }

    }

}



