
# UNFINISHED
#
# Seriously needs to be cleaned up.
#
# And put in a thing to notice if either up or down has no significant genes
# before it hits an error trying to get the expression matrix.
#
# And make sure not overwriting previous files.

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
# projectPath =  "/Users/nelsonlab/Documents/Results_temporarily_here/RORb_results/"
dataPath =  "/Users/nelsonlab/Documents/Results_temporarily_here/TTX_results/"
resPath = "/Users/nelsonlab/Documents/Results_temporarily_here/TTX_results/"
# projectPath = "/Users/emmamyers/Documents/Work_temp/"
# treatment = 'p200'
# treatment = "KO"
treatment = "TTX"
N = 100      # Set to > genes than are in the file with limma pvals to not limit how many genes are in each heatmap
ncolors = 5  # number of colors to have in palette
ylabSize = 8 # font size for gene symbols
figHeightPerGene = 20
figWidth = 300
lfcMin = 0  # Make 0 to not filter by LFC; I've been doing 1.5
alpha = 1 # Make 1 to not filter by p-value; I've been doing 0.01
scaleGenes = FALSE # If true, will standardize values within-gene to range from 0 to 1

getFileNames = function(path, comparison) {
  fileTPM = paste(path, comparison, "_TPM.csv", sep = '')
  fileLimma = paste(path, comparison, "_limma_ranked_genes.csv", sep = '')
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
  
  # Impose minimum LFC
  df = df[which(abs(df$LFC) >= lfcMin),]
  # Separate up from down
  dfDn = df[which(df$Direction == paste("dn", treatment, sep="")),]
  dfUp = df[which(df$Direction == paste("up", treatment, sep="")),]
  # Limit to significantly DEX genes
  dfDnSig = dfDn[which(abs(dfDn$p.value) < alpha),]
  dfUpSig = dfUp[which(abs(dfUp$p.value) < alpha),]
  # Limit to the top n by p-value, at most
  dfDnTop = dfDnSig[1:min(n, dim(dfDnSig)[1]),]
  dfUpTop = dfUpSig[1:min(n, dim(dfUpSig)[1]),]
  # # Don't worry about sorting by effect size now
  # lfcRanksDn = order(abs(dfDnTop$LFC), decreasing=TRUE)
  # lfcRanksUp = order(abs(dfUpTop$LFC), decreasing=TRUE)
  # dfDnTopSorted = dfDnTop[lfcRanksDn,]
  # dfUpTopSorted = dfUpTop[lfcRanksUp,]

  # Find this limma_ranked_genes files stuff in the expression matrix; those are the indexes to output
  topIdx = list(dnIdx = match(dfDnTop$Gene.symbol, rownames(exprMat)),
                      upIdx = match(dfUpTop$Gene.symbol, rownames(exprMat)))
  return(topIdx)
}

getLFCs = function(filename, geneSyms) {
  df = read.csv(filename)
  lfcs = df$LFC[match(geneSyms, df$Gene.symbol)]
  return(lfcs)
}

getPvals = function(filename, geneSyms) {
    df = read.csv(filename)
    pvals = df$p.value[match(geneSyms, df$Gene.symbol)]
    return(pvals)
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

makeHeatmap = function(exprForPlot, xlabs, figHeightPerGene, figWidth, ylabSize, comparison, direction) {
  # Figure out if minimum value is negative, since that affects how we want to set the lower limit of the color scale
  minFactor = 0.95; if (min(exprForPlot) < 0) {minFactor = 1.05}
  minVal = min(exprForPlot)*minFactor
  # Define plotly object
  figHeightThis = figHeightPerGene*dim(exprForPlot)[1]
  p = plot_ly(z = exprForPlot, x = xlabs, y = rownames(exprForPlot), type='heatmap', colors = colorRamp(c('yellow', 'red')), zmin = minVal, zmax = max(exprForPlot*minFactor), height=figHeightThis, width = figWidth)
  # p = layout(p, yaxis = list(tickfont = list(size = ylabSize), ticklen = 0), 
  #      xaxis = list(ticklen = 0), 
  #      title = paste(comparison, ',', direction))
  # # lenmode and len aren't actually affecting the color bar height
  # colorbar(p, limits = c(minVal, max(exprForPlot*minFactor)), lenmode = 'fraction',  len = 0.25, thickness = 30)
  return(p)
}



# This is the only code that should change for each comparison
# comparison = "Aging"
# comparison = "Rorb_p30"
comparison = "PVLate"


# Get TPMs and indexes in TPM matrix of top-ranking genes (up and down, separately)
writeLines("Getting expression, LFCs and p-values of genes to plot. . .", sep="")
filenames = getFileNames(dataPath, comparison)
tpms = getExprMat(filenames$fileTPM)
geneIdx = getTopIdx(tpms, filenames$fileLimma, N, lfcMin, treatment)
# load("Rorb_sharedIdx.Rdata")
# geneIdx = sharedIdx
geneLFCs = list()
geneLFCs$up = getLFCs(filenames$fileLimma, rownames(tpms)[geneIdx$upIdx])
geneLFCs$dn = getLFCs(filenames$fileLimma, rownames(tpms)[geneIdx$dnIdx])
genePvals = list()
genePvals$up = getPvals(filenames$fileLimma, rownames(tpms)[geneIdx$upIdx])
genePvals$dn = getPvals(filenames$fileLimma, rownames(tpms)[geneIdx$dnIdx])
writeLines("done.")


# Get values to display for up- and down-regulated genes
writeLines("Getting expression matrices...", sep="")
exprForPlotDn_unscaled = prepExprMat(tpms, geneIdx$dnIdx)
exprForPlotUp_unscaled = prepExprMat(tpms, geneIdx$upIdx)
writeLines("Done.")

# Normalize within-gene to 0-to-1 range, if requested
if (scaleGenes) {
    exprForPlotDn = t(apply(exprForPlotDn_unscaled, 1, normFun))
    exprForPlotUp = t(apply(exprForPlotUp_unscaled, 1, normFun))
} else {
    exprForPlotDn = exprForPlotDn_unscaled
    exprForPlotUp = exprForPlotUp_unscaled
}


# Write gene symbols to file
geneFn = paste(resPath, comparison, "_genes_in_plots.csv", sep="")
writeLines(paste("Writing selected gene symbols to ", geneFn, "...", sep=""), sep="")
upSyms = paste("up,", rev(rownames(exprForPlotUp)), ",", geneLFCs$up, ",", genePvals$up)
dnSyms = paste("dn,", rev(rownames(exprForPlotDn)), ",", geneLFCs$dn, ",", genePvals$dn)
con = file(geneFn)
writeLines(c(dnSyms, upSyms), con, sep=",\n")
close(con)
writeLines("done.")

# Get x-labels
# xlabs = sub("BF_RORbHT", "", colnames(tpms))
# xlabs = sub("BF_RORb", "", colnames(tpms))
xlabs = sub(comparison, '', colnames(tpms))

# My makeHeatmap() function won't do the layout and colorbar stuff so there's some more hackery going on here

# Heatmap for down-regulated genes
writeLines("Defining expression heatmap for down-regulated genes...", sep="")
pDn = makeHeatmap(exprForPlotDn, xlabs, figHeightPerGene, figWidth, ylabSize, comparison, paste('down-regulated in', treatment, 'samples'))
minFactor = 0.95; if (min(exprForPlotDn) < 0) {minFactor = 1.05}
minVal = min(exprForPlotDn)*minFactor
layout(pDn, yaxis = list(tickfont = list(size = ylabSize, tickvals = 1:dim(exprForPlotDn)[2]), ticklen = 0),
       xaxis = list(ticklen = 0),
       title = paste('down-regulated in', treatment, 'samples'))
# lenmode and len aren't actually affecting the color bar height
# colorbar(pDn, limits = c(minVal, max(exprForPlotDn*minFactor)), lenmode = 'fraction',  len = 0.25, thickness = 30)
writeLines("Done.")

# LFC
writeLines("Defining LFC heatmap for down-regulated genes...", sep="")
c = character(2); c[1] = '.'; c[2] = '.'
pLfcDn = plot_ly(z = cbind(rev(abs(geneLFCs$dn)), rev(abs(geneLFCs$dn))), type = 'heatmap', 
        x = rep('.', times=2), y = rownames(exprForPlotDn),
        colors = colorRamp(c('blue', 'green')), 
        height = figHeightPerGene*length(geneLFCs$dn), width = 200)
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
        height = figHeightPerGene*dim(exprForPlotDn_unscaled)[1], width = 200)
layout(pMeanDn, 
       yaxis = list(tickfont = list(size = ylabSize), ticklen = 0), 
       xaxis = list(ticklen = 0), 
       title = paste('Mean log(TPM), down in', treatment, 'samples'))
writeLines("Done.")

############################################################
# Heatmap for up-regulated genes
############################################################

writeLines("Defining expression heatmap for up-regulated genes...", sep="")
pUp = makeHeatmap(exprForPlotUp, xlabs, figHeightPerGene, figWidth, ylabSize, comparison, paste('up-regulated in', treatment, 'samples'))
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
pLfcUp = plot_ly(z = cbind(rev(abs(geneLFCs$up)), rev(abs(geneLFCs$up))), type = 'heatmap', 
        x = rep('.', times=2), y = rownames(exprForPlotUp),
        colors = colorRamp(c('blue', 'green')), 
        height = figHeightPerGene*length(geneLFCs$up), width = 200)
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
        height = figHeightPerGene*dim(exprForPlotUp_unscaled)[1], width = 200)
layout(pMeanUp, 
       yaxis = list(tickfont = list(size = ylabSize), ticklen = 0), 
       xaxis = list(ticklen = 0), 
       title = paste('Mean log(TPM), up in', treatment, 'samples'))
writeLines("Done.")



# # Doesn't include the stuff that was in layout and colorbar commands; doing that stuff in Illustrator
writeLines("Exporting pdfs for down-regulated genes...", sep="")
export(p=pDn, file=paste(resPath, comparison, '_dn.pdf', sep = ''))
export(p=pLfcDn, file=paste(resPath, comparison, '_dn_lfc.pdf', sep = ''))
export(p=pMeanDn, file=paste(resPath, comparison, '_dn_mean.pdf', sep = ''))
writeLines("Done.")

writeLines("Exporting pdfs for up-regulated genes...", sep="")
export(p=pUp, file=paste(resPath, comparison, '_up.pdf', sep = ''))
export(p=pLfcUp, file=paste(resPath, comparison, '_up_lfc.pdf', sep = ''))
export(p=pMeanUp, file=paste(resPath, comparison, '_up_mean.pdf', sep = ''))
writeLines("Done.")


