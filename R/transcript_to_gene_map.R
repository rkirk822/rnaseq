#' Create file associating transcript ids with gene ids
#'
#' Given a GTF file with both gene ids and transcript ids, write a csv with gene ids in the first column and transcript ids in the second.
#' @param inFile String - Name of GTF with transcript ids in gene_id field
#' @param outFile String - Name of csv with gene ids in the first column and transcript ids in the second
#' @return gene_transcript
#' @details Note that gene ids (in the first column of the output file) will be repeated when associated with multiple transcript ids.
#' TIME: A minute or two.
#' @examples
#' ids = transcript_to_gene_map("/Volumes/CodingClub1/RNAseq/Metadata/mm10_refGene.gtf", "/Volumes/CodingClub1/RNAseq/Metadata/RefSeq_id_map.csv")
#' @author Emma Myers
#' @export

transcript_to_gene_map = function(inFile, outFile) {

    # Check arguments
    if ( !file.exists(inFile) ) { stop("Input file does not exist.") }
    if ( file.exists(outFile) ) { stop("Output file already exists.") }

    # Load gtf file I usually use, that may or may not include introns
    writeLines("Reading gtf...", sep="")
    fullGtf = read.table(inFile, sep="\t")
    writeLines("done.")

    # Get the key-value field or whatever it's called, that includes gene ids and transcript ids
    fullCol9 = fullGtf[,9]
    fullCol9C = as.character(fullCol9)

    # Define a function that returns the gene and transcript id for each entry in the GTF
    # ASSUMING GENE ID IS IN FIRST BIT AND TRANSCRIPT ID IN SECOND
    get_ids = function(x) {
        splitStr = strsplit(x, ";")
        geneId = substr(splitStr[[1]][1], nchar("gene_id ")+1, nchar(splitStr[[1]][1]))
        trId = substr(splitStr[[1]][2], nchar(" transcript_id ")+1, nchar(splitStr[[1]][2]))
        return(c(geneId, trId))
    }

    # Apply that to the
    writeLines("Getting gene and transcript ids...", sep="")
    gene_transcript = t( sapply(fullCol9C, get_ids) )
    rownames(gene_transcript) = NULL
    gene_transcript = unique(gene_transcript)
    writeLines("done.")


    # Write to file
    writeLines("Writing gene and transcript ids to file...")
    write.table(gene_transcript, file=outFile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
    writeLines(paste("...done.  Filename: ", outFile, "\n"), sep="")

    return(gene_transcript)

}
