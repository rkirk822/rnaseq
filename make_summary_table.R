
#  make_summary_table.R
#
# This is jury-rigged for the TTX thing; would be nice to generalize but not now.
# Just do one table for each comparison in a separate file.  Then put them together manually.
# Also removal of the bad sample is hard-coded in.
#
# Columns should be:
# Gene symbols, ctl mean, ctl stderr, TTX mean, TTX stderr, effect size, p-value.

# Functions we'll be calling
source('/Volumes/CodingClub1/RNAseq/code/read_fcounts.R')
source('/Volumes/CodingClub1/RNAseq/code/counts_to_tpm.R')


setwd('/Volumes/CodingClub1/RNAseq/TTX/counts/counts_m20_q20_no_comment/')


# Get the gene ids and read counts (gene lengths for TPM are in here too)
# In this jury-rigging doing just one comparison per table
cellType = 'EMX'
stage = 'Late'
resfile = paste('/Users/nelsonlab/Documents/Toolboxes/limma_voom/TTX_results/', cellType, '_', stage, '_tpm_summary.csv', sep='')
fileList = list.files(pattern=cellType)
fileList = fileList[which(regexpr(stage, fileList)>0)]
if (any(fileList=='EMXTTXEarly_1_fcounts.txt')) {
    fileList = fileList[-which(fileList=='EMXTTXEarly_1_fcounts.txt')]
}
writeLines('Reading in counts. . .')
fcounts = read_fcounts(fileList)
writeLines('Done.')

# LIMMA-RELATED
# What genes did we keep in the limma analysis (because they're expressed somewhere)
dexTable = read.csv(paste('/Volumes/CodingClub1/RNAseq/TTX/limma_voom/', cellType, '_', stage, '_limma_ranked_genes.csv', sep=''), header=TRUE)
# dexTable = read.csv(paste('/Users/nelsonlab/Documents/Toolboxes/limma_voom/TTX_results/', cellType, '_', stage, '_limma_ranked_genes.csv', sep=''), header=TRUE)
emptyIdx = which(!is.element(fcounts$Geneid, dexTable$Gene.symbol))

# Want to have empty rows / NAs for genes that weren't kept.
geneSyms = fcounts$Geneid
countMat = fcounts$Count; countMat[emptyIdx,] = NA
geneLens = fcounts$Length; geneLens[emptyIdx] = NA
# Grab the pvals while we're at it (match returns index where each element of geneSyms appears in table)
pvals_limma=dexTable$p.value[match(geneSyms, dexTable$Gene.symbol)]

# Convert counts to TPM
exprMat = counts_to_tpm(countMat, geneLens)

# Get mean and standard error for ctl samples
ctlIdx = which(regexpr('Ctl', fileList)>0)
ctlMeans = rowMeans(exprMat[,ctlIdx])
ctlSdevs = apply(exprMat[,ctlIdx], 1, sd)
ctlSEs = ctlSdevs/sqrt(length(ctlIdx))

# Get mean and standard error for TTX samples
ttxIdx = which(regexpr('TTX', fileList)>0)
ttxMeans = rowMeans(exprMat[,ttxIdx])
ttxSdevs = apply(exprMat[,ttxIdx], 1, sd)
ttxSEs = ttxSdevs/sqrt(length(ttxIdx))

# Get log fold changes of mean values (control relative to TTX)
lfcs = log2(ctlMeans/ttxMeans)

# do this separately in case we want to change / not do it
ctlMeans = round(ctlMeans, digits=2)
ctlSEs = round(ctlSEs, digits=2)
ttxMeans = round(ttxMeans, digits=2)
ttxSEs = round(ttxSEs, digits=2)
lfcs = round(lfcs, digits=2)

# Make dataframe and write to csv file
df = data.frame(geneSyms, ctlMeans, ctlSEs, ttxMeans, ttxSEs, lfcs, pvals_limma)
colnames = c('Gene symbol', 'Mean (ctl)', 'SE (ctl)', 'Mean (ttx)', 'SE (ttx)', 'LFC (ctl mean:ttx mean)', 'p-value (limma)')
writeLines('Writing table. . .')
write.table(df, quote=FALSE, row.names=FALSE, col.names=colnames, file=resfile, sep=',')
writeLines('Done.')



