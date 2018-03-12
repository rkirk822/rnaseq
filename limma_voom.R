



# limma_voom.R
#
# This really needs some cleaning up.
#
# Use limma-voom to rank genes by changed expression across two conditions,
# using one as the reference in the linear models.  Save R objects to an Rdata file, write
# some results from the linear modelling to a text file, and write a ranked gene list to a csv
# file (ranked by BH-adjusted p-value).

# # Creagted dge with just the P2 samples and did: 
# design = model.matrix(~0+group)
# contr.matrix=makeContrasts(HTvsKO = groupHT-groupKO, levels=colnames(design))
# v = voom(dge_norm, design)                                 
# vfit = lmFit(v, design)                                              
# vfit = contrasts.fit(vfit, contrasts=contr.matrix)
# efit = eBayes(vfit)
# # Anton's suggestions: 
# You have defined 0 as the intercept (which is fine when you don't need or don't want to have a control group). 
# Then you defined your contrasts as HT - KO 
#   (i.e. genes that are WT>KO will be positive and genes that are KO>WT will be negative; 
#   this is OK if you keep this in mind, though I find it a bit counter-intuitive).
# If I had WT and KO samples, I would define WT as intercept and then look at the KO coefficients, 
#   but this is just my personal preference and either way is fine. 
#   That way you will be testing if genes are up- or downregulated in KO compared to WT, not vice versa.

# Now this is enough if just look at the genotype effect. 
# If you want to look at the effect of the developmental stage, 
#   I would define them as an additional variable in your model 
#   (you would have to decide if the variable should be discrete or continuous; 
#   I would go for discrete in your case). 
# So you can define the design matrix like this:
#   design = model.matrix(~gentype+stage)
# Where genotype is a factor with levels WT and KO (with WT as reference level), 
#   and stage would be a factor with levels P2, P7, P30 (e.g. with P2=youngest age as reference level).



# # Loading edgeR also loads limma.
library('edgeR')

# respath = '/Users/nelsonlab/Documents/Results_temporarily_here/TTX_results/'
# respath = '/Users/nelsonlab/Documents/Results_temporarily_here/RORb_results/'
# respath = "/Users/emmamyers/Documents/Work_temp/"
respath = "/Users/nelsonlab/Documents/Results_temporarily_here/Aging/"

# fn_ranks = paste(respath, 'TTX_pvals_all_data.csv', sep='')
# fn_ranks = paste(respath, 'Rorb_pvals_test_age.csv', sep='')
# fn_R = paste(respath, "Rorb_pvals_test_age.Rdata", sep="")
fn_ranks = paste(respath, "Aging_limma_ranked_genes.csv", sep="")
fn_R = paste(respath, "Aging_limma_ranked_genes.Rdata", sep="")

########################################################################
# Create the DGEList object and get some metadata in there
#########################################################################

# setwd('/Users/nelsonlab/Documents/Results_temporarily_here/TTX_results/counts_m20_q20_no_comment/')
# setwd("/Users/emmamyers/Documents/Work_temp/counts_no_comment_with_fakes/")
setwd("/Users/nelsonlab/Documents/Results_temporarily_here/Aging/counts_no_comment/")


fileList = list.files(pattern="_fcounts.txt")
# fileList = list.files(pattern="HT")
# fileList = fileList[ - which(regexpr('p2_', fileList) > 0 ) ]
# fileList = fileList[ - which(regexpr('p7_', fileList) > 0 ) ]
# fileList = fileList[ - which(regexpr('p200_', fileList) > 0 ) ]

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

# # For the thing you're comparing based on, set the control/reference level by making it the first thing in levels=c('grp1','grp2')
dge_all$samples$age = factor( sub('_.*','',sub('.*HT','',fileList)), levels = c("p30", "p200") )
# dge_all$samples$group = factor(substr(samplenames, nchar(cellType)+1, nchar(cellType)+3), levels = c('Ctl','TTX'))
# dge_all$samples$gentype = factor(sub('.*BF_RORb','',sub('p.*','',samplenames)), levels = c('HT','KO'))
# dge_all$samples$age = factor( sub('_.*','',substr(fileList, 10, nchar(fileList))), levels = c("p2", "p7", "p30", "p200") )

# # Get any other information about samples into the samples dataframe
# cellType=rep('EMX',times=length(fileList)); cellType[ which( substr(fileList, 1, 2) == "PV" ) ] = "PV"
# dge_all$samples$cellType = cellType
# stage=rep('Early',times=length(fileList)); stage[ which( regexpr('Late', fileList) > 0 ) ] = 'Late'
# dge_all$samples$stage = stage
# dge_all$samples$age = as.factor(as.numeric(sub('_.*','',sub('.*p','',fileList))))


####################################################################
# Couple more things to do with the DGEList object:
# Restrict to expressed genes
# Calculate normalization factors
####################################################################

# Get row indexes of genes you're going to keep
# Get TPM and get gene symbols out of first column and into row names
tpmMin = 10
# countsTable = read.csv('/Users/emmamyers/Documents/Work_temp/Rorb_aging_new_TPM.csv')
countsTable = read.csv("/Users/nelsonlab/Documents/Results_temporarily_here/Aging/Aging_TPM.csv")
rownames(countsTable) = countsTable[,1]
countsTable = countsTable[,-1]
# Get gene means for each group of samples.
group1 = countsTable[, which(regexpr('p30_', colnames(countsTable)) > 0)]
group2 = countsTable[, which(regexpr('p200_', colnames(countsTable)) > 0)]
# group3 = countsTable[, which(regexpr('p2_', colnames(countsTable)) > 0)]
# group4 = countsTable[, which(regexpr('p7_', colnames(countsTable)) > 0)]
# orvec = rowMeans(group1)>tpmMin | rowMeans(group2)>tpmMin | rowMeans(group3)>tpmMin | rowMeans(group4)>tpmMin
keepLogical = cbind( rowMeans(group1) > tpmMin,  rowMeans(group2) > tpmMin )
orvec = keepLogical[,1] | keepLogical[,2]
keepIdx = which(orvec)
# Now you have the row indexes of genes to keep; keep them in dge
dge = dge_all[keepIdx,,keep.lib.sizes=FALSE]
print('Restricted to expressed genes.')


# Calculate normalization factors.  That's now part of the DGEList object.
dge = calcNormFactors(dge, method='TMM')
print('Calculated normalization factors.')
# # Probably want to have a look at these too.
# dge_norm$samples

####################################################################
# Create the design and contrast matrices
####################################################################

# Where genotype is a factor with levels WT and KO (with WT as reference level), 
#   and stage would be a factor with levels P2, P7, P30 (e.g. with P2=youngest age as reference level).

# # Design matrix
comparison = dge$samples$age
# var1 = dge$sample$age
design = model.matrix(~comparison)  # Control as intercept ("keeping" intercept)
# design = model.matrix(~comparison+var1)  # Control as intercept ("keeping" intercept)
# design = model.matrix(~0+comparison)  # 0 as intercept ("dropping" intercept)
colnames(design) = gsub('comparison', '', colnames(design))
colnames(design) = gsub('var1', '', colnames(design))
# gentype = dge_all$samples$gentype
# design = model.matrix(~gentype)

# # Contrast matrix
# # If you left in the intercept, remember your ref condition is called "(Intercept)".
# contrastMat = makeContrasts( HTvsKO = KO - HT, levels=colnames(design) )
# contrastMat = makeContrasts( p30p200 = p200 - p30, levels = colnames(design) )
# contrastMat = makeContrasts(
#     PVEarly = PVTTXEarly - PVCtlEarly,
#     PVLate = PVTTXLate - PVCtlLate,
#     EMXEarly = EMXTTXEarly - EMXCtlEarly,
#     EMXLate = EMXTTXLate - EMXCtlLate,
#     levels = colnames(design))


#####################################################################################
# Get precision weights, linear models, variability estimates,
# coefficients for contrasts in contrast matrix, and adjusted p-values
#####################################################################################

# Get "precision weights" to deal with heteroscedascity in the count data
voomOutput = voom(dge, design)
print('Got weights.')

# Fit a linear model to each gene
vfit = lmFit(voomOutput, design)
print('Fitted linear model to each gene.')

# # Ccoefficients for contrasts
# vfit = contrasts.fit(vfit, contrasts=contrastMat)
# print('Calculated coefficients for desired contrasts.')

# Better gene-wise variability estimates
efit = eBayes(vfit)
print('Estimated gene-wise variability.')

# FDR-adjusted p-values in summary table
summaryTable = topTable(efit, n=Inf, sort.by='none')

####################################################################
# Save stuff to an R object
####################################################################

# save(file=fn_R, "dge_all", "dge", "design", "contrastMat", 
#      "voomOutput", "vfit", "efit", "summaryTable")

# Trying without contrast matrix first
save(file=fn_R, "dge_all", "dge", "design", 
     "voomOutput", "vfit", "efit", "summaryTable")



###################################################################################
# Write to csv the following stuff in the following format, for one comparison: 
# File is names <comparison>_limma_ranked_genes.csv
# Columns names and content, in this order:
# Direction: "dn<non-reference condition>" or "up<non-reference condition>"
# Gene symbol
# LFC - This is the limma-adjusted one, based on its model
# Rank - By p-value
# p-value
# 
# Also, sort so that the genes going down in the non-reference condition are 
# first, followed by those that go up.  Within those group, they're sorted by rank.
####################################################################################

idx = 2 # THIS IS WHAT WON'T WORK WHEN EXCLUDING INTERCEPT
nonRefName = "p200"  # GOTTA FIGURE THIS OUT TOO

# Get ranks; sort the gene symbols, LFCs, and p-values accordingly
pval_rank = order(summaryTable$adj.P.Val)
gene = rownames(summaryTable)[pval_rank]
lfc = efit$coefficients[pval_rank, idx]
pval = summaryTable$adj.P.Val[pval_rank]

# Get values for "Direction" column
dn_idx = which(lfc<0)
up_idx = which(lfc>0)
direction = lfc
direction[up_idx] = paste('up', nonRefName, sep='')
direction[dn_idx] = paste('dn', nonRefName, sep='')

# Going to separate down-regulated from up-regulated genes
new_idx = c(dn_idx, up_idx)

# Put it all in a data frame
# I STILL WANT TO TRIPLE-CHECK IF PVAL_RANK[NEW_IDX] IS CORRECT
genes_by_rank = data.frame(direction[new_idx],  gene[new_idx], lfc[new_idx], pval_rank[new_idx], pval[new_idx])
colnames = c('Direction', 'Gene symbol', 'LFC', 'Rank', 'p-value')

# Write the table
write.table(genes_by_rank, quote=FALSE, row.names=FALSE, col.names=colnames, file=fn_ranks, sep=',')
print('Wrote ranked gene lists to file.')


