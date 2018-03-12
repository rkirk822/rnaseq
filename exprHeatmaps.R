
# UNFINISHED
# Works, but need to clean it up.
# Separate the stuff for selecting the genes from the stuff for making the plots.  I want one function that takes a gene 
# list, identfies them in the given files, and makes these figures for them.

# First dealing with the fact that fewer genes end up in the heatmaps than have fold changes above the minimum according to topTable.

# Given files with TPM/sample and p-values for differential expression across conditions, create and print to pdf expression 
# heatmaps of the top-ranked genes by p-value. Each comparison has two expression heatmaps: one for the top-ranked genes 
# that decrease in the treatment condition, and one for those that increase.  Genes are ordered by effect size (i.e. log fold
# change), decreasing.  Each comparison also gets a one-column heatmap of effect size for each of the top genes, and another 
# one of their mean expression across all samples for the comparison.
# P-values are based on limma-voom.
# Expression values are log2(TPM).
# Color scales are not equal across plots.
# Formatted similarly to Fig. 3 in Okaty et al 2009, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2749660/#R70.

library(plotly)

# # Define some variables
# projectPath =  '/Users/nelsonlab/Documents/Results_temporarily_here/Aging/'
projectPath =  "/Users/nelsonlab/Documents/Results_temporarily_here/RORb_results/"
# projectPath = '/Users/emmamyers/Documents/Work_temp/Aging/'
# treatment = 'p200'
treatment = "KO"
N = 100000      # Set to > genes than are in the file with limma pvals to not limit how many genes are in each heatmap
ncolors = 5  # number of colors to have in palette
ylabSize = 8 # font size for gene symbols
figHeight = 2000
figWidth = 300
lfcMin = 1.5  # Make 0 to not filter by LFC
alpha = 0.01 # Make 1 to not filter by p-value

getFileNames = function(path, comparison) {
  # comparison = paste(cellType, stage, sep = '_')
  fileTPM = paste(projectPath, comparison, '_TPM.csv', sep = '')
  fileLimma = paste(projectPath, comparison, '_limma_ranked_genes.csv', sep = '')
  filenames = list(fileTPM=fileTPM, fileLimma=fileLimma)
}

getExprMat = function(filename) {
  # Input: Csv file containing expression matrix
  # Output: Expression matrix, with gene symbols as row names
  exprMat = read.csv(filename)
  rownames(exprMat) = exprMat[,1]  # Move the gene symbols into row names
  exprMat = exprMat[,-1] # Get rid of first column, that has gene symbols
  return(exprMat)
  }

getTopIdx = function(exprMat, filename, n, lfcMin, treatment) {
  # Inputs:  Expression matrix, file with pvals from limma, number of genes
  # Outputs: Data frame with indexes (in expr matrix) of the 100 top-ranked genes down-regulated in treated samples; same for up-regulated, and the effect size (log fold change) for each.  All ordered by effect size.
  df = read.csv(filename)
  dfDn = df[which(df$Direction == paste("dn", treatment, sep="")),]
  dfUp = df[which(df$Direction == paste("up", treatment, sep="")),]
  dfDnTop = dfDn[1:min(n, dim(dfDn)[1]),]
  dfUpTop = dfUp[1:min(n, dim(dfUp)[1]),]
  lfcRanksDn = order(abs(dfDnTop$LFC), decreasing=TRUE)
  lfcRanksUp = order(abs(dfUpTop$LFC), decreasing=TRUE)
  dfDnTopSorted = dfDnTop[lfcRanksDn,]
  dfUpTopSorted = dfUpTop[lfcRanksUp,]
  # Impose minimum LFC
  dfDnTopSorted = dfDnTopSorted[which(abs(dfDnTopSorted$LFC) >= lfcMin),]
  dfUpTopSorted = dfUpTopSorted[which(abs(dfUpTopSorted$LFC) >= lfcMin),]
  # Limit to significantly DEX genes
  dfDnTopSorted = dfDnTopSorted[which(abs(dfDnTopSorted$p.value) < alpha),]
  dfUpTopSorted = dfUpTopSorted[which(abs(dfUpTopSorted$p.value) < alpha),]
  # Find this limma_ranked_genes files stuff in the expression matrix; those are the indexes to output
  topIdx = list(dnIdx = match(dfDnTopSorted$Gene.symbol, rownames(exprMat)),
                      upIdx = match(dfUpTopSorted$Gene.symbol, rownames(exprMat)),
                      dnLFC = dfDnTopSorted$LFC, upLFC = dfUpTopSorted$LFC)
  return(topIdx)
}

# For use in prepExprMat.  Or actually maybe right before plotting instead.
normFun = function(v) {vnorm = (v-min(v)) / (max(v)-min(v)); return(vnorm)}

prepExprMat = function(exprMat, idx) {
  # Inputs: Expression matrix, row indexes
  # Outputs:  Matrix of log2 expression, with 0 values where expression was 0, vertically flipped for heatmap
  exprMatL2 = log2(exprMat[idx,])
  exprMatL2[exprMat[idx,] == 0] = 0
  exprMatL2Rev = apply(exprMatL2, 2, rev)
  return(exprMatL2Rev)
}

makeHeatmap = function(exprForPlot, xlabs, figHeight, figWidth, ylabSize, comparison, direction) {
  # Figure out if minimum value is negative, since that affects how we want to set the lower limit of the color scale
  minFactor = 0.95; if (min(exprForPlot) < 0) {minFactor = 1.05}
  minVal = min(exprForPlot)*minFactor
  # Define plotly object
  p = plot_ly(z = exprForPlot, x = xlabs, y = rownames(exprForPlot), type='heatmap', colors = colorRamp(c('yellow', 'red')), zmin = minVal, zmax = max(exprForPlot*minFactor), height=figHeight, width = figWidth)
  # p = layout(p, yaxis = list(tickfont = list(size = ylabSize), ticklen = 0), 
  #      xaxis = list(ticklen = 0), 
  #      title = paste(comparison, ',', direction))
  # # lenmode and len aren't actually affecting the color bar height
  # colorbar(p, limits = c(minVal, max(exprForPlot*minFactor)), lenmode = 'fraction',  len = 0.25, thickness = 30)
  return(p)
}



# This is the only code that should change for each comparison
# comparison = "Aging"
comparison = "Rorb_p2"
# comparison = "EMX_Late"


# Get TPMs and indexes in TPM matrix of top-ranking genes (up and down, separately)
filenames = getFileNames(projectPath, comparison)
tpms = getExprMat(filenames$fileTPM)
topIdx = getTopIdx(tpms, filenames$fileLimma, N, lfcMin, treatment)

# Get values to display for up- and down-regulated genes
# Normalize within-gene to 0-to-1 range
writeLines("Getting expression matrices...", sep="")
exprForPlotDn_unscaled = prepExprMat(tpms, topIdx$dnIdx)
exprForPlotUp_unscaled = prepExprMat(tpms, topIdx$upIdx)
exprForPlotDn = t(apply(exprForPlotDn_unscaled, 1, normFun))
exprForPlotUp = t(apply(exprForPlotUp_unscaled, 1, normFun))
writeLines("Done.")

# Get x-labels
# xlabs = sub("BF_RORbHT", "", colnames(tpms))
xlabs = sub("BF_RORb", "", colnames(tpms))
# xlabs = sub(cellType, '', sub(paste(stage, '_', sep=''), '', colnames(tpms)))

# My makeHeatmap() function won't do the layout and colorbar stuff so there's some more hackery going on here

# Heatmap for down-regulated genes
writeLines("Defining expression heatmap for down-regulated genes...", sep="")
pDn = makeHeatmap(exprForPlotDn, xlabs, figHeight, figWidth, ylabSize, comparison, paste('down-regulated in', treatment, 'samples'))
minFactor = 0.95; if (min(exprForPlotDn) < 0) {minFactor = 1.05}
minVal = min(exprForPlotDn)*minFactor
layout(pDn, yaxis = list(tickfont = list(size = ylabSize), ticklen = 0), 
       xaxis = list(ticklen = 0), 
       title = paste('down-regulated in', treatment, 'samples'))
# lenmode and len aren't actually affecting the color bar height
# colorbar(pDn, limits = c(minVal, max(exprForPlotDn*minFactor)), lenmode = 'fraction',  len = 0.25, thickness = 30)
writeLines("Done.")

# LFC
writeLines("Defining LFC heatmap for down-regulated genes...", sep="")
c = character(2); c[1] = '.'; c[2] = '.'
pLfcDn = plot_ly(z = cbind(rev(abs(topIdx$dnLFC)), rev(abs(topIdx$dnLFC))), type = 'heatmap', 
        x = rep('.', times=2), y = rownames(exprForPlotDn),
        colors = colorRamp(c('blue', 'green')), 
        height = figHeight, width = 200)
layout(pLfcDn, 
       yaxis = list(tickfont = list(size = ylabSize), ticklen = 0), 
       xaxis = list(ticklen = 0), 
       title = paste('LFC, down in', treatment, 'samples'))
writeLines("Done.")

# Mean log2(TPM)
writeLines("Defining mean expression heatmap for down-regulated genes...", sep="")
c = character(2); c[1] = '.'; c[2] = '.'
meanColRamp = c('turquoise1', 'magenta')
# meanColRamp = c('skyblue', 'magenta')
# meanColRamp = c('turquoise', 'violet')
pMeanDn = plot_ly(z = cbind(rowMeans(exprForPlotDn_unscaled), rowMeans(exprForPlotDn_unscaled)), type = 'heatmap', 
        x = rep('.', times=2), y = rownames(exprForPlotDn),
        colors = colorRamp(rev(meanColRamp)), 
        height = figHeight, width = 200)
layout(pMeanDn, 
       yaxis = list(tickfont = list(size = ylabSize), ticklen = 0), 
       xaxis = list(ticklen = 0), 
       title = paste('Mean log(TPM), down in', treatment, 'samples'))
writeLines("Done.")

############################################################
# Heatmap for up-regulated genes
############################################################

writeLines("Defining expression heatmap for up-regulated genes...", sep="")
pUp = makeHeatmap(exprForPlotUp, xlabs, figHeight, figWidth, ylabSize, comparison, paste('up-regulated in', treatment, 'samples'))
minFactor = 0.95; if (min(exprForPlotUp) < 0) {minFactor = 1.05}
minVal = min(exprForPlotUp)*minFactor
layout(pUp, yaxis = list(tickfont = list(size = ylabSize), ticklen = 0), 
       xaxis = list(ticklen = 0), 
       title = paste('down-regulated in', treatment, 'samples'))
# lenmode and len aren't actually affecting the color bar height
# colorbar(pUp, limits = c(minVal, max(exprForPlotUp*minFactor)), lenmode = 'fraction',  len = 0.25, thickness = 30)
writeLines("Done.")

# LFC
writeLines("Defining LFC heatmap for up-regulated genes...", sep="")
c = character(2); c[1] = '.'; c[2] = '.'
pLfcUp = plot_ly(z = cbind(rev(abs(topIdx$upLFC)), rev(abs(topIdx$upLFC))), type = 'heatmap', 
        x = rep('.', times=2), y = rownames(exprForPlotUp),
        colors = colorRamp(c('blue', 'green')), 
        height = figHeight, width = 200)
layout(pLfcUp, 
       yaxis = list(tickfont = list(size = ylabSize), ticklen = 0), 
       xaxis = list(ticklen = 0), 
       title = paste('LFC, up in', treatment, 'samples'))
writeLines("Done.")

# Mean log2(TPM)
writeLines("Defining mean expression heatmap for up-regulated genes...", sep="")
c = character(2); c[1] = '.'; c[2] = '.'
# meanColFun = colorRampPalette(c('turquoise1', 'magenta'))
# meanCols = meanColFun(N)
meanColRamp = c('turquoise1', 'magenta')
# meanColRamp = c('skyblue', 'magenta')
# meanColRamp = c('turquoise', 'violet')
pMeanUp = plot_ly(z = cbind(rowMeans(exprForPlotUp_unscaled), rowMeans(exprForPlotUp_unscaled)), type = 'heatmap', 
        x = rep('.', times=2), y = rownames(exprForPlotUp),
        colors = colorRamp(rev(meanColRamp)), 
        height = figHeight, width = 200)
layout(pMeanUp, 
       yaxis = list(tickfont = list(size = ylabSize), ticklen = 0), 
       xaxis = list(ticklen = 0), 
       title = paste('Mean log(TPM), up in', treatment, 'samples'))
writeLines("Done.")



# stop("Not creating pdfs")

# # Doesn't include the stuff that was in layout and colorbar commands; doing that stuff in Illustrator
writeLines("Exporting pdfs for down-regulated genes...", sep="")
export(p=pDn, file=paste(projectPath, comparison, '_dn.pdf', sep = ''))
export(p=pLfcDn, file=paste(projectPath, comparison, '_dn_lfc.pdf', sep = ''))
export(p=pMeanDn, file=paste(projectPath, comparison, '_dn_mean.pdf', sep = ''))
writeLines("Done.")

writeLines("Exporting pdfs for up-regulated genes...", sep="")
export(p=pUp, file=paste(projectPath, comparison, '_up.pdf', sep = ''))
export(p=pLfcUp, file=paste(projectPath, comparison, '_up_lfc.pdf', sep = ''))
export(p=pMeanUp, file=paste(projectPath, comparison, '_up_mean.pdf', sep = ''))
writeLines("Done.")


