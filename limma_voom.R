


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
# day = 2
respath = '/Users/nelsonlab/Documents/Toolboxes/limma_voom/RORb_results/testing/'
fn_efit = paste(respath,'RORb_p',day,'_limma_efit.txt',sep='')
# fn_ranks = paste(respath,'RORb_p',day,'_limma_ranked_genes.csv',sep='')
fn_ranks = paste(respath,'RORb_all_limma_ranked_genes.csv',sep='')
fn_R = paste(respath,'RORb_p',day,'_limma.R',sep='')

# # fileList - files containing sample counts to be included
# setwd('/Volumes/CodingClub1/RNAseq/TTX/counts_m20_q20_no_comment/')
# setwd('/Volumes/DataStorage2/Emma/From_CodingClub1/RNAseq/TTX/counts/counts_m20_q20_no_comment/')
# fileList = list.files(pattern=cellType)
# fileList = fileList[which(regexpr(stage, fileList)>0)]
# if (any(fileList=='EMXTTXEarly_1_fcounts.txt')) {
#     fileList = fileList[-which(fileList=='EMXTTXEarly_1_fcounts.txt')]
# }
setwd('/Users/nelsonlab/Documents/Toolboxes/limma_voom/counts_no_comment/')
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

################
# # This is specific to making one large table for the RORb data
app=FALSE
use_colnames=TRUE
if (day>2) {
    app=TRUE
    use_colnames=FALSE
}
#################

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
# genes_by_rank = data.frame(direction_sorted[new_idx],  genes_sorted[new_idx], lfc_sorted[new_idx], rank[new_idx], pvals_sorted[new_idx])
# colnames = c('Direction', 'Gene symbol', 'LFC', 'Rank', 'p-value')
# # # for doing one big table of RORb # # #
genes_by_rank = data.frame(rep(paste('p',day,sep=''),times=length(genes)), direction_sorted[new_idx],  genes_sorted[new_idx], lfc_sorted[new_idx], rank[new_idx], pvals_sorted[new_idx])
if (use_colnames) {
    colnames = c('Age', 'Direction', 'Gene symbol', 'LFC', 'Rank', 'p-value')
} else {
    colnames = FALSE
}
# # # end for doing one big table of RORb # # #
write.table(genes_by_rank, quote=FALSE, row.names=FALSE, col.names=colnames, append = app, file=fn_ranks, sep=',')
print('Wrote ranked gene lists to file.')

# Save stuff to an R object
save(file=fn_R, 'dge_all', 'dge', 'design', 'v', 'vfit', 'efit', 'dt', 'genes_by_rank', 'pval_ranks', 'lfcMin')



#########################################################################
# # Stuff that will or might be used as I improve the script/pipeline
#########################################################################

# # Partway through ranking up-genes against each other and down-genes against each other
# # Save ranked lists of genes.  When things access the 2nd col it's to get at KO, not intercept.
# # Separate genes that go up in KO samples from those that go down.
# direction = efit$t[,2]
# # Sort each list in ascending order of p-value
# dn_pvals = efit$p.value[which(direction<0),2])
# dn_ranks = order(dn_pvals)
# dn_pvals_sorted = dn_pvals(dn_ranks)
# # Effect size
# dn_lfc = efit$coefficients[dnidx,2]
# dn_lfc_sorted = dn_lfc[dn_ranks]
# # Gene symbols
# dn_genes = rownames(efit$p.value[dnidx])
# dn_genes_sorted = dn_genes[dn_ranks]



# # I haven't looked much at what treat does.  It makes an extreme difference
# # in that there are far fewer genes making the cutoff.  It doesn't include F-values.
# # Let's go with normal LFC thresholding for now.
# # Use treat method to impose minimum LFC
# tfit=treat(vfit, lfc=1);
# print('Used treat method to impose minimum LFC.')


# # Look at the output
# koUp=which(dt==-1)
# koDn=which(dt==1)



# Show me (or save) table of genes with zero counts in all samples.
# Plots, possibly one function each:
# Density of log cpm

########################
# Log cpm density plots
########################
# library(RColorBrewer)
# nsamples=ncol(dge)
# col=brewer.pal(nsamples,'Paired')
# par(mfrow=c(1,2))
# # Plot first sample.  Ylim will actually need to be higher, for the other samples.
# plot( density(dge_lcpm[,1]), col=col[1], lwd=2, las=2, main='', xlab='', ylim=c(0,0.23))
# title(main='Log CPM, all genes')
# abline(v=0,lty=3)
# # Plot rest of samples
# for (i in 2:nsamples) {
# 	den = density(dge_lcpm[,i])
# 	lines(den$x, den$y, col=col[i], lwd=2)
# 	}
# #  Legend is huge, don't bother til you look up how to shrink it
# legend('topright', samplenames, text.col=col, bty='n')

# # In the other panel, plot density of log cpm for expressed genes only
# dge_sub_cpm=cpm(dge_sub)
# dge_sub_lcpm=cpm(dge_sub,log=TRUE)
# plot( density(dge_sub_lcpm[,1]), col=col[1], lwd=2, las=2, main='', xlab='', ylim=c(0,0.25))
# title(main='Log CPM, expressed genes')
# abline(v=0,lty=3)
# for (i in 2:nsamples) {
# 	den = density(dge_sub_lcpm[,i])
# 	lines(den$x, den$y, col=col[i], lwd=2)
# 	}

# ################################
# # Sample distribution box plots
# ################################
# # One with and one without using the scaling factors
# par(mfrow=c(1,2))
# boxplot(dge_sub_lcpm, las=2, col=col, main="")
# title(main="Log cpm, expressed genes, unnormalized",ylab='Log cpm')
# dge_norm_lcpm=cpm(dge_norm,log=TRUE)
# boxplot(dge_norm_lcpm, las=2, col=col, main="")
# title(main="Log cpm, expressed genes, normalized",ylab='Log cpm')
