


# limma_voom.R
#
# Use limma-voom to rank genes by changed expression across two conditions,
# using one as the reference in the linear models.  Save R objects to an Rdata file, write
# some results from the linear modelling to a text file, and write a ranked gene list to a csv
# file (ranked by BH-adjusted p-value).
#
# TIME:  Under a minute, for one RORb age (comparing 2 genotypes with 4 samples each).
#
# Important notes!!!!!!!!!!!!!!!!!!!
#
# 1) THIS SCRIPT DOES NOT CURRENTLY CHECK BEFORE SAVING RESULTS FILES.
# IF YOU PICK A RESULTS FILE NAME FOR A FILE THAT ALREADY EXISTS, IT WILL BE OVERWRITTEN.
#
# 2) There are some hacks in here for getting the RORb gene ranks from all 3 ages into one big table.
# I really should clean that up, but meanwhile:
# FOR GETTING ALL RORB AGES IN ONE FILE
# To have it loop through all ages and put ranked gene lists into one csv file:
# - comment out the definition of the variable "day",
# - uncomment the "if (day>2)" statement,
# - make the filename for the ranked lists file be appropriate,
# - down where genes_by_rank and colnames are defined, chose the ones that include "day" as first column,
# - and do this:
# > days=c(2,7,30)
# > for (day in days) { source('/path/to/script/dex_script.R') }
#
# 3) In the fcounts output files, the commented-out line at the top has to be deleted,
# and the headers left in.  Here's a few lines of R code for removing the comment:
# setwd('/path/to/featureCounts/output/')
# fileList=list.files()
# # remove files likely to be in that directory that we don't want
# fileList=fileList[which(regexpr('summary',fileList)<0)]
# fileList=fileList[which(regexpr('display',fileList)<0)]
# fileList=fileList[which(regexpr('pdf',fileList)<0)]
# for (f in fileList) {
#     t = read.table(f)
#     write.table(t, paste('/destination/of/fixed/files/',f,sep=''),quote=FALSE, row.names=FALSE,col.names=FALSE,sep='\t')
# }
#

# # Loading edgeR also loads limma.
library('edgeR')

########################################################################
# This part you have to really look at and make sure it's what you want.
# Define lfcMin, fn_efit, fn_ranks, fn_R.
# Define dge_all, and put in gene and sample metadata.
# Define the design matrix "design".
#########################################################################

# Minimum log fold change to be considered for differential expression
lfcMin = 1
# Minimum total CPM a gene needs, across all samples, to be expressed "somewhere"
cpmMin = 4


# # Results filenames:  fn_efit, fn_ranks, fn_R
# # Make sure the results file names correctly represents the files going into the DGE object
# cellType = 'PV'; stage = 'Early'
# respath = '/Users/nelsonlab/Documents/Toolboxes/limma_voom/TTX_results/'
# fn_efit = paste(respath, cellType, '_', stage, '_limma_efit.txt', sep='')
# fn_ranks = paste(respath, cellType, '_', stage, '_limma_ranked_genes.csv', sep='')
# fn_R = paste(respath, cellType,'_', stage, '_limma.Rdata', sep='')
day = 30
respath = '/Users/nelsonlab/Documents/Results_temporarily_here/RORb_results/'
fn_efit = paste(respath,'RORb_p',day,'_limma_efit.txt',sep='')
fn_ranks = paste(respath,'RORb_p',day,'_limma_ranked_genes.csv',sep='')
# fn_ranks = paste(respath,'RORb_all_limma_ranked_genes.csv',sep='')
fn_R = paste(respath,'RORb_p',day,'_limma.R',sep='')

# # fileList - files containing sample counts to be included
# setwd('/Volumes/CodingClub1/RNAseq/TTX/counts_m20_q20_no_comment/')
# setwd('/Volumes/DataStorage2/Emma/From_CodingClub1/RNAseq/TTX/counts/counts_m20_q20_no_comment/')
# fileList = list.files(pattern=cellType)
# fileList = fileList[which(regexpr(stage, fileList)>0)]
# if (any(fileList=='EMXTTXEarly_1_fcounts.txt')) {
#     fileList = fileList[-which(fileList=='EMXTTXEarly_1_fcounts.txt')]
# }
setwd('/Users/nelsonlab/Documents/Results_temporarily_here/RORb_results/counts_no_comment/')
fileList = list.files(pattern=paste('p',day,'_',sep='')) # include the underscore to be safe

# # Create the DGEList object (digital gene expression) dge_all
# # "all" is for all genes; we're going to subset later and have "dge"
print('Creating DGEList object')
dge_all = readDGE(fileList,columns=c(1,7))
print('Done.')

# Gene metadata in a dataframe
dge_all$genes = as.data.frame(rownames(dge_all))

# Column names are sample names
samplenames = gsub('_fcounts.txt','',fileList)
colnames(dge_all) = samplenames

# # Get information about samples into the samples dataframe
dge_all$samples$age = as.factor(as.numeric(sub('_.*','',sub('.*p','',fileList))))

# # For the thing you're comparing based on, set the control/reference level by making it the first thing in levels=c('grp1','grp2')
# dge_all$samples$group = factor(substr(samplenames, nchar(cellType)+1, nchar(cellType)+3), levels = c('Ctl','TTX'))
dge_all$samples$gentype = factor(sub('.*BF_RORb','',sub('p.*','',samplenames)), levels = c('HT','KO'))

# # Create design matrix
# group = dge_all$samples$group
# design = model.matrix(~group)
gentype = dge_all$samples$gentype
design = model.matrix(~gentype)

# ################
# # # This is specific to making one large table for the RORb data
app=FALSE
use_colnames=TRUE
# if (day>2) {
#     app=TRUE
#     use_colnames=FALSE
# }
# #################

####################################################################
####################################################################

# Restrict to genes with at least X cpm in at least Y samples total.
dge_cpm = cpm(dge_all)
keep.genes = rowSums(dge_cpm>1)>=cpmMin
dge = dge_all[keep.genes,,keep.lib.sizes=FALSE]
print('Restricted to expressed genes.')
# # Want to have it print this to screen, and the actual numbers.
# dim(dge)[1]/dim(dge_all)[1]

# Calculate normalization factors.  That's now part of the DGEList object.
dge = calcNormFactors(dge, method='TMM')
print('Calculated normalization factors.')
# # Probably want to have a look at these too.
# dge_norm$samples

# Deal with heteroscedascity in the count data
# Get "precision weights"
v = voom(dge, design)
print('Got weights.')
# Fit a linear model to each gene
vfit = lmFit(v, design)
print('Fitted linear model to each gene.')
# Get better gene-wise variability estimates
efit = eBayes(vfit)
print('Estimated gene-wise variability.')

# Summary of numbers of differentially expressed genes
dt=decideTests(efit, lfc=lfcMin)
summary(dt)

# Write results to a file
write.fit(efit, dt, file=fn_efit)
print('Wrote efit to file.')

# Get FDR-adjusted p-values
tt = topTable(efit, n=Inf, sort.by='none')

# Save ranked list of genes
pval_ranks = order(tt$adj.P.Val) # second column for non-reference group, not intercept
pvals_sorted = tt$adj.P.Val[pval_ranks]
lfc_sorted=efit$coefficients[pval_ranks,2] # effect size
# Include gene names
genes = rownames(tt)
genes_sorted = genes[pval_ranks]
# Include a column ("Direction") indicating whether expr went up or down in non-reference group
direction = efit$t[,2]
# (little song and dance to figure out what the non-reference group was called)
comparisonNameLen = nchar(names(attr(design, 'contrasts')))
comparisonNameAndNonRefName = dimnames(design)[[2]][[2]]
nonRefName = substr(comparisonNameAndNonRefName, comparisonNameLen+1, nchar(comparisonNameAndNonRefName))
# Okay, NOW include the freakin column
direction[which(direction>0)] = paste('up', nonRefName, sep='')
direction[which(direction<0)] = paste('dn', nonRefName, sep='')
direction_sorted = direction[pval_ranks]
# # Make sure this is true
# all(direction=='up'|| direction=='dn')
# They're going to be in order of rank and I want that in there
rank = 1:dim(dge)[1]
# Pull out down versus up
tvals_sorted = efit$t[pval_ranks,2]
dn_idx = which(tvals_sorted<0)
up_idx = which(tvals_sorted>0)
new_idx = c(dn_idx, up_idx)
# Put it all in a data frame
genes_by_rank = data.frame(direction_sorted[new_idx],  genes_sorted[new_idx], lfc_sorted[new_idx], rank[new_idx], pvals_sorted[new_idx])
colnames = c('Direction', 'Gene symbol', 'LFC', 'Rank', 'p-value')
# # # for doing one big table of RORb # # #
# genes_by_rank = data.frame(rep(paste('p',day,sep=''),times=length(genes)), direction_sorted[new_idx],  genes_sorted[new_idx], lfc_sorted[new_idx], rank[new_idx], pvals_sorted[new_idx])
# if (use_colnames) {
#     colnames = c('Age', 'Direction', 'Gene symbol', 'LFC', 'Rank', 'p-value')
# } else {
#     colnames = FALSE
# }
# # # end for doing one big table of RORb # # #
write.table(genes_by_rank, quote=FALSE, row.names=FALSE, col.names=colnames, append = app, file=fn_ranks, sep=',')
print('Wrote ranked gene lists to file.')

# Save stuff to an R object
save(file=fn_R, 'dge_all', 'dge', 'design', 'v', 'vfit', 'efit', 'dt', 'genes_by_rank', 'pval_ranks', 'lfcMin')




