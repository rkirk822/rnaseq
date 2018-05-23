#' Convert bedgraph values between positive and negative
#'
#' Get raw counts of reads mapped to specified genomic features using featureCount tool from the subread package.
#' @param inFiles Character - Bedgraph file list
#' @param outDest String - Directory where output files should be saved
#' @param outSuffix String - will be appended to original filename; default "_flipped"
#' @details Take a list of bedgraph files and convert positive values to negative and vice-versa.  This is a hack for getting
#' value from two different files into the same deeptools heatmap, but with a different color scheme for each
#' (using computeMatrix() and plotHeatmap()).  Handy for being able to distinguish forward- and reverse-strand reads.
#' TIME:  3-5 s per bedgraph, at least for ~15-20MB ones.
#' @examples
#' bgs = list.files(pattern = "minus.bedgraph")
#' bedgraph_flip(bgs)
#' @author Emma Myers
#' @export

bedgraph_flip = function(inFiles, outDest="./", outSuffix="_flipped") {

    # Check arguments
    outDest = dir_check(outDest)
    if ( !all(file.exists(inFiles)) ) {
        writeLines("Missing input files:")
        writeLines(inFiles[which(!file.exists(inFiles))])
        stop("Missing input file(s).  See above.")
    }
    if ( isEmpty(inFiles) ) { stop("No input files given.") }

    # Don't want it writing table with scientific notation
    options(scipen=999)

    for (f in inFiles) {

        writeLines(paste("Processing file:", f))

        fOut = paste(outDest, gsub(".bedgraph", paste(outSuffix, ".bedgraph", sep=""), basename(f)), sep="")

        # Check if there's already an output file with this name
        if ( file_checks(fOut, shouldExist=FALSE, verbose=TRUE) ) {
            writeLines(paste("File will be saved to:", fOut))
          # Read, flip, and write
            bedgraphTable = read.table(f)
            bedgraphTable[,4] = bedgraphTable[,4] * -1
            # 1st and 2nd columns only are delimited by a tab; the rest by spaces
            colsAfterFirst = apply(bedgraphTable[,2:dim(bedgraphTable)[2]], 1, function(x) {paste(x, collapse=" ")})
            names(colsAfterFirst) = NULL
            bedgraphTable[,2] = colsAfterFirst
            bedgraphTable = bedgraphTable[, c(1,2)]
            write.table(bedgraphTable, file=fOut, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
        }

    }

}

