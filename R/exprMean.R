#' Mean expression heatmaps
#'
#' Create single column heatmap of mean expression values (across samples)
#' Needs plotly
#' @param exprDataFrame Data frame - Gene x sample expression values (counts, tpm, whatever)
#' @param genes Character - Gene symbols that appear as row names in exprDataFrame
#' @param samples Character - Sample names that appear as column names in exprDataFrame
#' @param L2 Logical - Whether to take log2 of values
#' @param yticklabSize Numeric - Font size for gene symbols
#' @param figHeightPerGene Numeric -
#' @param figWidth Numeric -
#' @param colorsPlot Color ramp -
#' @param ncolors Numeric -
#' @param plotTitle String -
#' @param minVal Numeric -
#' @param maxVal Numeric -
#' @param fileOut String - If given, save png to with this filename
#' @return Plotly object
#' @details Given a gene-by-sample dataframe with expression values, and (optionally) a list of the genes and
#' samples you want included, make a one-column heatmap of mean expression using plotly.
#' Note that, probably because I'm not set up to do the paying thing, some stuff doesn't get
#' incorporated into the output png.
#' @author Emma Myers
#' @export

exprMean = function(exprDataFrame, genes=NULL, samples=NULL, L2=FALSE,
                yticklabSize=8, figHeightPerGene=20, figWidth = 200,
                colorsPlot = colorRamp(c("turquoise1", "magenta")),
                ncolors=5, plotTitle=NULL, minVal=NULL, maxVal=NULL, fileOut=NULL) {

    ### Check inputs ######################################################

    # If genes and samples unspecified, use them all
    if (is.null(genes)) { genes = rownames(exprDataFrame) }
    if (is.null(samples)) { samples = colnames(exprDataFrame) }

    # Check for missing genes and samples
    msgGeneIdx = which( !is.element(genes, rownames(exprDataFrame)) )
    msgSampleIdx = which( !is.element(samples, colnames(exprDataFrame)) )
    if ( any(c(msgGeneIdx, msgSampleIdx)) ) {
        writeLines("Missing genes:")
        writeLines(genes[msgGeneIdx])
        writeLines("Missing samples:")
        writeLines(samples[msgSampleIdx])
        stop("At least one gene or sample requested is not in exprDataFrame (see above).")
    }


    ### Get submatrix #######################################################
    exprSubset = exprDataFrame[genes, samples] # hang onto these values
    exprForPlot = exprSubset # we're gonna transform these for the plot


    ### Do some transformations ##############################################
    # Take log2 if requested
    if ( L2 ) {
        exprTemp = exprForPlot
        exprForPlot = log2(exprForPlot)
        # zeros where expression was zero
        exprForPlot[exprTemp == 0] = 0
    }

    # Vertically flip for heatmap
    exprForPlot = apply(exprForPlot, 2, rev)

    ### Dimensions and colorscale stuff ######################################
    figHeightThis = figHeightPerGene*dim(exprForPlot)[1]
    # Min value for getting the color scale should be slightly smaller than actual min in data
    if ( is.null(minVal) ) {
        if (min(exprForPlot) >= 0) { minFactor = 0.95 } else {minFactor = 1.05}
        minVal = min(exprForPlot) * minFactor
    }
    # and vice versa
    if ( is.null(maxVal) ) {
        if (max(exprForPlot) >= 0) { maxFactor = 1.05 } else {maxFactor = 0.95}
        maxVal = max(exprForPlot) * maxFactor
    }

    ### Making the plot #######################################################
    # Make plotly object
    meanPlotlyObj = plotly::plot_ly(z = cbind(rowMeans(exprForPlot), rowMeans(exprForPlot)),
        x = rep('.', times=2), y = rownames(exprForPlot),
        type='heatmap', colors = colorsPlot, height=figHeightThis, width = figWidth)
    # Set some layout stuff
    meanPlotlyObj = plotly::layout(meanPlotlyObj,
                           yaxis = list(tickfont = list(size = yticklabSize, ticklen = 0)),
                           xaxis = list(ticklen = 0),
                           title = plotTitle)



    # Save if given fileOut name
    if ( !is.null(fileOut) ) {
        fileOutFull = fileOut
        if ( !(substr(fileOutFull, nchar(fileOutFull)-3, nchar(fileOutFull)) == ".png") ) { fileOutFull = paste(fileOutFull, ".png", sep="") }
        if ( !is.null(fileOut) ) {
            if ( file_checks(fileOutFull, shouldExist=FALSE, verbose=TRUE) ) {
                writeLines(paste("Saving image to", fileOutFull))
                plotly::plotly_IMAGE(meanPlotlyObj, format="png", out_file=fileOutFull)
            }
        }
    }

    return(meanPlotlyObj)

# # These don't retain stuff you do with layout():
# htmlwidgets::saveWidget(as_widget(exprPlotlyObj), "/Users/nelsonlab/Documents/exprTest.html")
# export(p=exprPlotlyObj, file="/Users/nelsonlab/Documents/exprTest.pdf")

}
