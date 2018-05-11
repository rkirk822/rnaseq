#' Log2 fold change heatmap
#'
#' Create single column heatmap of log2 fold change across two conditions
#' Needs plotly
#' @param exprDataFrame Data frame - Gene x sample expression values (counts, tpm, whatever)
#' @param sampleGroups List - Two elements, each containing indices of samples in exprDataFrame
#' @param genes Character? -
#' @param ylabSize Numeric - Font size for gene symbols
#' @param figHeightPerGene Numeric -
#' @param figWidth Numeric -
#' @param colorsPlot Color ramp -
#' @param ncolors Numeric -
#' @param plotTitle String -
#' @param minVal Numeric -
#' @param maxVal Numeric -
#' @param fileOut String - If given, save png to with this filename
#' @param addVal Numeric - Very small value added to expression values to avoid taking log2 of 0
#' @return Plotly object
#' @details  Given a gene-by-sample dataframe with expression values, a list with two vectors indicating
#' samples in each of two groups, and (optionally) a list of the genes you want included,
#' make a one-column heatmap of log fold change across the group means for each gene, using plotly.
#' Note that, probably because I'm not set up to do the paying thing, some stuff doesn't get
#' incorporated into the output png.
#' Fold change will be group2 / group1 (so if you put a control condition in group 1, the sign
#' will be positive if expression increased in the treatment condition, not that that changes
#' anything in the heatmap since we'll show abs(LFC)).
#' The value of addVal will be added to all expression levels to avoid dividing by zero or taking
#' the log of zero.  Default is 0.1 because I looked at one sample's distribution of TPM (excluding
#' zero values) to see what constitutes a tiny amount, and 0.1 was under the 3rd percentile.
#' @examples
#' exprData = read.table("Rorb_p2_TPM.csv", header=TRUE, row.names=1, sep=",")
#' colnames(exprDataFrame) = gsub("BF_RORb", "", colnames(exprDataFrame))
#' geneList = c("Rorb", "Plxnd1", "Has2", "Sparcl1", "Pde1a", "Has3")
#' groups = vector("list", 2)
#' groups[[1]] = colnames(exprDataFrame)[1:4]
#' groups[[2]] = colnames(exprDataFrame)[5:8]
#' exprLFC(exprData, genes=geneList, sampleGroups=groups, fileOut="exprLFC.png")
#' @author Emma Myers
#' @export


exprLFC = function(exprDataFrame, sampleGroups, genes=NULL,
                ylabSize=8, figHeightPerGene=20, figWidth = 200,
                colorsPlot = colorRamp(c("blue", "green")), ncolors=5,
                plotTitle="Log2 FC", minVal=NULL, maxVal=NULL, fileOut=NULL,
                addVal=0.1) {

    ### Check inputs ######################################################

    # If genes and samples unspecified, use them all
    if (is.null(genes)) { genes = rownames(exprDataFrame) }

    # Check for missing genes and samples
    msgGeneIdx = which( !is.element(genes, rownames(exprDataFrame)) )
    msgSampleIdx = which( !(sapply(sampleGroups, function(x, y) {is.element(x,y)}, y=colnames(exprDataFrame))) )
    if ( any(c(msgGeneIdx, msgSampleIdx)) ) {
        writeLines("Missing genes:")
        writeLines(genes[msgGeneIdx])
        writeLines("Missing samples:")
        writeLines(samples[msgSampleIdx])
        stop("At least one gene or sample requested is not in exprDataFrame (see above).")
    }


    ### Get LFCs #############################################################
    # avoid dividing by or taking log of zero
    exprDataFrame = exprDataFrame + addVal
    groupMeans = vector("list", 2)
    groupMeans[[1]] = rowMeans(exprDataFrame[genes, sampleGroups[[1]]])
    groupMeans[[2]] = rowMeans(exprDataFrame[genes, sampleGroups[[2]]])
    lfc = log2( groupMeans[[2]] / groupMeans[[1]] )

    ### Do some transformations ###############################################

    # Get absolute value and vertically flip for heatmap
    lfc = rev(abs(lfc))

    ### Dimensions and colorscale stuff ######################################
    figHeightThis = figHeightPerGene*length(lfc)
    # Min value for getting the color scale should be slightly smaller than actual min in data
    if ( is.null(minVal) ) {
        if (min(lfc) >= 0) { minFactor = 0.95 } else {minFactor = 1.05}
        minVal = min(lfc) * minFactor
    }
    # and vice versa
    if ( is.null(maxVal) ) {
        if (max(lfc) >= 0) { maxFactor = 1.05 } else {maxFactor = 0.95}
        maxVal = max(lfc) * maxFactor
    }

    ### Making the plot #######################################################
    # Make plotly object
    lfcPlotlyObj = plot_ly(z = cbind(lfc, lfc), x = rep('.', times=2), y = names(lfc),
                    type='heatmap', colors = colorsPlot, zmin = minVal, zmax = maxVal,
                    height=figHeightThis, width = figWidth)
    # Set some layout stuff
    lfcPlotlyObj = layout(lfcPlotlyObj,
                        yaxis = list(tickfont = list(size = ylabSize, ticklen = 0),
                        xaxis = list(ticklen = 0),
                        title = plotTitle))


    # Save if given fileOut name
    if ( !is.null(fileOut) ) {
        fileOutFull = fileOut
        if ( !(substr(fileOutFull, nchar(fileOutFull)-3, nchar(fileOutFull)) == ".png") ) { fileOutFull = paste(fileOutFull, ".png", sep="") }
        if ( !is.null(fileOut) ) {
            if ( file_checks(fileOutFull, shouldExist=FALSE, verbose=TRUE) ) {
                writeLines(paste("Saving image to", fileOutFull))
                plotly_IMAGE(lfcPlotlyObj, format="png", out_file=fileOutFull)
            }
        }
    }

    return(lfcPlotlyObj)

# # These don't retain stuff you do with layout():
# htmlwidgets::saveWidget(as_widget(exprPlotlyObj), "/Users/nelsonlab/Documents/exprTest.html")
# export(p=exprPlotlyObj, file="/Users/nelsonlab/Documents/exprTest.pdf")

}
