

# exprHeatmap.R
#
# Given a gene-by-sample dataframe with expression values, and (optionally) a list of the genes and
# samples you want included, make an expression heatmap using plotly.
#
# Note that, probably because I'm not set up to do the paying thing, some stuff doesn't get
# incorporated into the output png.

# USAGE
# exprData = read.table("Rorb_p2_TPM.csv", header=TRUE, row.names=1, sep=",")
# colnames(exprDataFrame) = gsub("BF_RORb", "", colnames(exprDataFrame))
# geneList = c("Rorb", "Plxnd1", "Has2", "Sparcl1", "Pde1a", "Has3")
# sampleList = c("HTp2_1", "HTp2_2", "KOp2_1", "KOp2_2")
# exprHeatmap(exprData, genes=geneList, samples=sampleList, fileOut="expr.png")

exprHeatmap = function(exprDataFrame, genes=NULL, samples=NULL, L2=TRUE, scaleGenes=TRUE, scaleByGroup=NULL,  
                ylabSize=8, figHeightPerGene=20, figWidth = 300, colorsPlot = colorRamp(c("yellow", "red")), ncolors=5,
                plotTitle="Expression heatmap", minVal=NULL, maxVal=NULL, fileOut=NULL) {
                    
    # Unnecessary once packaged
    library(plotly)
    # source("/Users/nelsonlab/Documents/Toolboxes/rna-seq/file_checks.R")
    source("/Users/work/Documents/rna-seq/file_checks.R")

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
    if ( scaleGenes ) {
        # If we're scaling the genes but weren't given control samples, scale to 0-to-1 range
        if (is.null(scaleByGroup) ) {
            normFun = function(v) {vnorm = (v-min(v)) / (max(v)-min(v)); return(vnorm)}
            exprForPlot = t(apply(exprForPlot, 1, normFun))
        } else {
            exprForPlot = exprForPlot - rowMeans(exprForPlot[,scaleByGroup])
        }
    } else {
        if (!scaleGenes) { stop("You gave me a group of samples to use in scaling expression, but scaleGenes is FALSE.") }
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
    exprPlotlyObj = plot_ly(z = exprForPlot, x = colnames(exprForPlot), y = rownames(exprForPlot),
                    type='heatmap', colors = colorsPlot,
                    zmin = minVal, zmax = maxVal,
                    height=figHeightThis, width = figWidth)

    # Set some layout stuff
    exprPlotlyObj = layout(exprPlotlyObj,
                    yaxis = list(tickfont = list(size = ylabSize, tickvals = 1:dim(exprForPlot)[2]), ticklen = 0),
                    xaxis = list(ticklen = 0),
                    title = plotTitle)
    

    # Save if given fileOut name
    if ( !is.null(fileOut) ) {
        fileOutFull = fileOut
        if ( !(substr(fileOutFull, nchar(fileOutFull)-3, nchar(fileOutFull)) == ".png") ) { fileOutFull = paste(fileOutFull, ".png", sep="") }
        if ( !is.null(fileOut) ) {
            if ( file_checks(fileOutFull, shouldExist=FALSE, verbose=TRUE) ) {
                writeLines(paste("Saving image to", fileOutFull))
                plotly_IMAGE(exprPlotlyObj, format="png", out_file=fileOutFull)
            }
        }
    }
    
    return(exprPlotlyObj)

# # These don't retain stuff you do with layout():
# htmlwidgets::saveWidget(as_widget(exprPlotlyObj), "/Users/nelsonlab/Documents/exprTest.html")
# export(p=exprPlotlyObj, file="/Users/nelsonlab/Documents/exprTest.pdf")

}
