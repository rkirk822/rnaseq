#' Replace transcript ids with gene ids in GTF file
#'
#' Given a GTF file with transcript ids in the "gene_id" field and a file mapping these to gene ids, write a new GTF with gene ids in that field.
#' @param gtf String - Name of GTF with transcript ids in gene_id field
#' @param mapFile String - Name of csv with gene ids in the first column and transcript ids in the second
#' @details  Use the transcript_to_gene_map function to create the mapFile input.
#' TIME:  Most of it's the "Replacing transcript ids with gene ids in key-value field" part, which is anywhere from 15-30m.
#' @examples
#' transcript_to_gene_fix("/Volumes/CodingClub1/RNAseq/Metadata/mm10_refSeq_introns.gtf", "/Volumes/CodingClub1/RNAseq/Metadata/RefSeq_id_map.csv")
#' @author Emma Myers
#' @export

transcript_to_gene_fix = function(gtf, mapFile, outSuffix="_geneids") {

    # Check arguments (and define output file name)
    if ( !file.exists(gtf) ) { stop("GTF input file does not exist.") }
    if ( !file.exists(mapFile) ) { stop("Map file does not exist.") }
    outFile = paste( tools::file_path_sans_ext(gtf), outSuffix, ".gtf", sep="" )
    if ( file.exists(outFile) ) { stop("Output file already exists.") }

    gene_transcript = read.csv(mapFile, header=FALSE, stringsAsFactors=FALSE)

    writeLines("Reading gtf with transcript ids as gene ids...", sep="")
    fullGtf = read.table(gtf, sep="\t", quote="")
    writeLines("done.")

    # Get the key-value field or whatever it's called, that includes gene ids and transcript ids
    oldCol9 = fullGtf[,9]
    oldCol9C = as.character(oldCol9)

    # Will use sapply to apply this function to each row of fullCol9C
    # t will be the table with real gene ids in the first column and transcript ids in the second
    replace_ids = function(x, t) {
        # Find the actual gene id
        splitStr = strsplit(x, ";")
        trName = substr(splitStr[[1]][2], nchar(" transcript_id ")+1, nchar(splitStr[[1]][2]))
        trNameSplit = strsplit(trName, "_intron")
        trId = substr(trNameSplit[[1]][1], 2, nchar(trNameSplit[[1]][1])) # removing escaped quote
        gId = t[,1][which(t[,2]==trId)]
        # If this transcript id isn't in the big gtf, don't replace it with an empty character string
        if (length(gId)==0) { gId = trId }
        gId = paste("\"", gId, "\"", sep="")
        # Put it in the right place in the original input
        return(gsub(paste("gene_id", trName), paste("gene_id", gId), x))
    }

    # The long part (could be 15m, could be 30)
    writeLines("Replacing transcript ids with gene ids in key-value field...", sep="")
    newCol9C = sapply(oldCol9C, replace_ids, t=gene_transcript)
    writeLines("done.")

    # Now we need to insert that back into fullGtf as the new ninth column
    fullGtfNew = fullGtf
    fullGtfNew[,9] = newCol9C

    # And write it to file
    writeLines("Writing new gtf to file...", sep="")
    write.table(fullGtfNew, file=outFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    writeLines(paste("done.  File is:", outFile))

}
