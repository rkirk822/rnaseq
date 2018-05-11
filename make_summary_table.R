

# Works but needs to be cleaned up and made into a function.
# And column names of csv, near the bottom, still use "ttx" as treatment
# condition name and have to be changed manually!  Fix that.

#  make_summary_table.R
#
# Just do one table for each comparison in a separate file.  Then put them together manually.
#
# Columns should be:
# Gene symbols, ctl mean, ctl stderr, treatment mean, treatment stderr, effect size, p-value.

# Functions we'll be calling
# source('/Users/emmamyers/Documents/Work_temp/code/read_fcounts.R')
# source('/Users/emmamyers/Documents/Work_temp/code/counts_to_tpm.R')
# source('/Volumes/CodingClub1/RNAseq/code/read_fcounts.R')
# source('/Volumes/CodingClub1/RNAseq/code/counts_to_tpm.R')
source('/Users/nelsonlab/Documents/Toolboxes/rna-seq/read_fcounts.R')
source('/Users/nelsonlab/Documents/Toolboxes/rna-seq/counts_to_tpm.R')

# setwd('/Users/emmamyers/Documents/Work_temp/Aging/counts_no_comment/')
# setwd('/Volumes/CodingClub1/RNAseq/TTX/counts/counts_m20_q20_no_comment/')
setwd('/Users/nelsonlab/Documents/Smintheus_stuff_copied_here/TTX_stuff/counts_m20_q20_no_comment/')

# Get the gene ids and read counts (gene lengths for TPM are in here too)
# In this jury-rigging doing just one comparison per table
cellType = 'PV'
stage = 'Late'
resfile = paste('/Users/nelsonlab/Documents/Results_temporarily_here/TTX_results/', cellType, '_', stage, '_tpm_summary.csv', sep='')
fileList = list.files(pattern=cellType)
fileList = fileList[which(regexpr(stage, fileList)>0)]
# # No longer removing bad TTX sample, but that's what this was
#  if (any(fileList=='EMXTTXEarly_1_fcounts.txt')) {
#     fileList = fileList[-which(fileList=='EMXTTXEarly_1_fcounts.txt')]
# }
# resfile = '/Users/emmamyers/Documents/Work_temp/Aging/Rorb_aging_new_tpm_summary.csv'
# fileList = list.files(pattern='fcounts')
writeLines('Reading in counts. . .')
fcounts = read_fcounts(fileList)
writeLines('Done.')

# LIMMA-RELATED
# What genes did we keep in the limma analysis (because they're expressed somewhere)
# dexTable = read.csv(paste('/Volumes/CodingClub1/RNAseq/TTX/limma_voom/', cellType, '_', stage, '_limma_ranked_genes.csv', sep=''), header=TRUE)
dexTable = read.csv(paste('/Users/nelsonlab/Documents/Results_temporarily_here/TTX_results/', cellType, stage, '_limma_ranked_genes.csv', sep=''), header=TRUE)
# dexTable = read.csv('/Users/emmamyers/Documents/Work_temp/Aging/Rorb_aging_new_limma_ranked_genes.csv')
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
# ctlIdx = which(regexpr('p30', fileList)>0)
ctlIdx = which(regexpr('Ctl', fileList)>0)
ctlMeans = rowMeans(exprMat[,ctlIdx])
ctlSdevs = apply(exprMat[,ctlIdx], 1, sd)
ctlSEs = ctlSdevs/sqrt(length(ctlIdx))

# Get mean and standard error for treatment samples
# tmIdx = which(regexpr('p200', fileList)>0)
tmIdx = which(regexpr('TTX', fileList)>0)
tmMeans = rowMeans(exprMat[,ttxIdx])
tmSdevs = apply(exprMat[,ttxIdx], 1, sd)
tmSEs = ttxSdevs/sqrt(length(ttxIdx))

# Get log fold changes of mean values (control relative to TTX)
lfcs = log2(ctlMeans/tmMeans)

# do this separately in case we want to change / not do it
ctlMeans = round(ctlMeans, digits=2)
ctlSEs = round(ctlSEs, digits=2)
tmMeans = round(tmMeans, digits=2)
tmSEs = round(tmSEs, digits=2)
lfcs = round(lfcs, digits=2)

# Make dataframe and write to csv file
df = data.frame(geneSyms, ctlMeans, ctlSEs, tmMeans, tmSEs, lfcs, pvals_limma)
colnames = c('Gene symbol', 'Mean (ctl)', 'SE (ctl)', 'Mean (ttx)', 'SE (ttx)', 'LFC (ctl mean:ttx mean)', 'p-value (limma)')
writeLines('Writing table. . .')
write.table(df, quote=FALSE, row.names=FALSE, col.names=colnames, file=resfile, sep=',')
writeLines('Done.')



