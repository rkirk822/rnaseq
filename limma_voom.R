
# BUG:  The "pval_rank" column in the csv files has the wrong values.
# Things are in the right order, though.
#
# Also this really needs some cleaning up.

# limma_voom.R
#
# Use limma-voom to rank genes by changed expression across two conditions,
# using one as the reference in the linear models.  Save R objects to an Rdata file and,
# for each contrast in the contrast matrix, write a ranked gene list to a csv file
# (ranked by BH-adjusted p-value).
#
# If you didn't use a contrast matrix, go down to the end where the csv file gets written
# and adjust accordingly.




# # Loading edgeR also loads limma.
library("edgeR")
source("/Users/nelsonlab/Documents/Toolboxes/rna-seq/file_checks.R")

respath = "/Users/nelsonlab/Documents/Results_temporarily_here/TTX_results/"
# respath = "/Users/nelsonlab/Documents/Results_temporarily_here/RORb_results/"
# respath = "/Users/emmamyers/Documents/Work_temp/"
# respath = "/Users/nelsonlab/Documents/Results_temporarily_here/Aging/"


fn_R = paste(respath, "TTX_limma_ranked_genes.Rdata", sep="")
# fn_R = paste(respath, "Rorb_test_limma_ranked_genes.Rdata", sep="")
# fn_R = paste(respath, "Aging_limma_ranked_genes.Rdata", sep="")

########################################################################
# Create the DGEList object and get some metadata in there
#########################################################################

setwd("/Users/nelsonlab/Documents/Smintheus_stuff_copied_here/TTX_stuff/counts_m20_q20_no_comment/")
# setwd("/Users/emmamyers/Documents/Work_temp/counts_no_comment_with_fakes/")
# setwd("/Users/nelsonlab/Documents/Results_temporarily_here/Aging/counts_no_comment/")
# setwd("/Users/nelsonlab/Documents/Results_temporarily_here/RORb_results/counts_no_comment/")

fileList = list.files(pattern="_fcounts.txt")
# fileList = list.files(pattern="HT")
# fileList = fileList[ - which(regexpr('p2_', fileList) > 0 ) ]
# fileList = fileList[ - which(regexpr('p7_', fileList) > 0 ) ]
# fileList = fileList[ - which(regexpr('p200_', fileList) > 0 ) ]

# # Create the DGEList object (digital gene expression) dge_all
# # "all" is for all genes; we're going to subset later and have "dge"
writeLines("Creating DGEList object...", sep="")
dge_all = readDGE(fileList,columns=c(1,7))
writeLines("Done.")

# Gene metadata in a dataframe
dge_all$genes = as.data.frame(rownames(dge_all))

# Column names are sample names
samplenames = gsub("_fcounts.txt","",fileList)
colnames(dge_all) = samplenames

# # Get sample metadata into the DGEList object.
# # For the thing you're comparing based on, set the control/reference level by making it the first thing in levels=c('grp1','grp2')

# dge_all$samples$age = factor( sub('_.*','',sub('.*HT','',fileList)), levels = c("p30", "p200") )

# dge_all$samples$gentype = factor(sub('.*BF_RORb','',sub('p.*','',samplenames)), levels = c('HT','KO'))
# dge_all$samples$age = factor( sub('_.*','',substr(fileList, 10, nchar(fileList))), levels = c("p2", "p7", "p30", "p200") )
# dge_all$samples$subgroup = factor( paste(dge_all$samples$gentype, dge_all$samples$age, sep=""),
#     levels = c("HTp2", "KOp2", "HTp7", "KOp7", "HTp30", "KOp30") )
# dge_all$samples$age = as.factor(as.numeric(sub('_.*','',sub('.*p','',fileList))))


cellType=rep('EMX',times=length(fileList)); cellType[ which( substr(fileList, 1, 2) == "PV" ) ] = "PV"
dge_all$samples$cellType = cellType
stage=rep('Early',times=length(fileList)); stage[ which( regexpr('Late', fileList) > 0 ) ] = 'Late'
dge_all$samples$stage = stage
dge_all$samples$conditions = factor(substr(samplenames, nchar(cellType)+1, nchar(cellType)+3), levels = c('Ctl','TTX'))
dge_all$samples$subgroup = factor(paste(dge_all$samples$cellType, dge_all$samples$conditions, dge_all$samples$stage, sep=""), levels = c("EMXCtlEarly", "EMXTTXEarly", "EMXCtlLate", "EMXTTXLate", "PVCtlEarly", "PVTTXEarly", "PVCtlLate", "PVTTXLate"))

####################################################################
# Couple more things to do with the DGEList object:
# Restrict to expressed genes
# Calculate normalization factors
####################################################################

# Get row indexes of genes you're going to keep
# Get TPM and get gene symbols out of first column and into row names
tpmMin = 20
# countsTable = read.csv('/Users/emmamyers/Documents/Work_temp/Rorb_aging_new_TPM.csv')
# countsTable = read.csv("/Users/nelsonlab/Documents/Results_temporarily_here/RORb_results/Rorb_TPM.csv")
countsTable = read.csv("/Users/nelsonlab/Documents/Smintheus_stuff_copied_here/TTX_stuff/tpm/TTX_project_TPM.csv")
rownames(countsTable) = countsTable[,1]
countsTable = countsTable[,-1]

# # Get gene means for each group of samples.
# group1 = countsTable[, which(regexpr('p30_', colnames(countsTable)) > 0)]
# # group2 = countsTable[, which(regexpr('p200_', colnames(countsTable)) > 0)]
# group3 = countsTable[, which(regexpr('p2_', colnames(countsTable)) > 0)]
# group4 = countsTable[, which(regexpr('p7_', colnames(countsTable)) > 0)]

group1 = countsTable[, which(regexpr("EMXCtlEarly", colnames(countsTable)) > 0)]
group2 = countsTable[, which(regexpr("EMXTTXEarly", colnames(countsTable)) > 0)]
group3 = countsTable[, which(regexpr("EMXCtlLate", colnames(countsTable)) > 0)]
group4 = countsTable[, which(regexpr("EMXTTXLate", colnames(countsTable)) > 0)]
group5 = countsTable[, which(regexpr("PVCtlEarly", colnames(countsTable)) > 0)]
group6 = countsTable[, which(regexpr("PVTTXEarly", colnames(countsTable)) > 0)]
group7 = countsTable[, which(regexpr("PVCtlLate", colnames(countsTable)) > 0)]
group8 = countsTable[, which(regexpr("PVTTXLate", colnames(countsTable)) > 0)]



# # Get indexes of genes with at least one mean over the minimum.
# orvec = rowMeans(group1)>tpmMin | rowMeans(group3)>tpmMin | rowMeans(group4)>tpmMin
orvec = rowMeans(group1)>tpmMin | rowMeans(group2)>tpmMin | rowMeans(group3)>tpmMin | rowMeans(group4)>tpmMin | rowMeans(group5)>tpmMin | rowMeans(group6)>tpmMin | rowMeans(group7)>tpmMin | rowMeans(group8)>tpmMin
# keepLogical = cbind( rowMeans(group1) > tpmMin,  rowMeans(group2) > tpmMin )
# orvec = keepLogical[,1] | keepLogical[,2]
keepIdx = which(orvec)

# Now you have the row indexes of genes to keep; keep them in dge
dge = dge_all[keepIdx,,keep.lib.sizes=FALSE]
print("Restricted to expressed genes.")


# Calculate normalization factors.  That's now part of the DGEList object.
dge = calcNormFactors(dge, method="TMM")
print("Calculated normalization factors.")
# # Probably want to have a look at these too.
# dge_norm$samples

####################################################################
# Create the design and contrast matrices
####################################################################

# Where genotype is a factor with levels WT and KO (with WT as reference level), 
#   and stage would be a factor with levels P2, P7, P30 (e.g. with P2=youngest age as reference level).

# # Design matrix
comparison = dge$samples$subgroup
# var1 = dge$sample$age
# design = model.matrix(~comparison)  # Control as intercept ("keeping" intercept)
# design = model.matrix(~comparison+var1)  # Control as intercept ("keeping" intercept)
design = model.matrix(~0+comparison)  # 0 as intercept ("dropping" intercept)
colnames(design) = gsub("comparison", "", colnames(design))
colnames(design) = gsub("var1", "", colnames(design))
# gentype = dge_all$samples$gentype
# design = model.matrix(~gentype)

# # Contrast matrix, nonRefName
# # If you left in the intercept, remember your ref condition is called "(Intercept)".
# contrastMat = makeContrasts( HTvsKO = KO - HT, levels=c("HT", "KO") )
# contrastMat = makeContrasts(
#     p2 = KOp2 - HTp2,
#     p7 = KOp7 - HTp7,
#     p30 = KOp30 - HTp30,
#     levels = colnames(design))
# nonRefName = "KO"  # For writing to csv file
# contrastMat = makeContrasts( p30p200 = p200 - p30, levels = colnames(design) )
# nonRefName = "p200"
contrastMat = makeContrasts(
    PVEarly = PVTTXEarly - PVCtlEarly,
    PVLate = PVTTXLate - PVCtlLate,
    EMXEarly = EMXTTXEarly - EMXCtlEarly,
    EMXLate = EMXTTXLate - EMXCtlLate,
    levels = colnames(design))
nonRefName = "TTX"


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
vfit = contrasts.fit(vfit, contrasts=contrastMat)
print('Calculated coefficients for desired contrasts.')

# Better gene-wise variability estimates
efit = eBayes(vfit)
print('Estimated gene-wise variability.')

# FDR-adjusted p-values in summary table
# Note they include all contrasts - if you want to know for individual ones you have to specificy with coef arg, as below
summaryTable = topTable(efit, n=Inf, sort.by='none')

####################################################################
# Save stuff to an R object
####################################################################

# save(file=fn_R, "dge_all", "dge", "design", "contrastMat", 
#      "voomOutput", "vfit", "efit", "summaryTable")

# Trying without contrast matrix first
if (file.exists(fn_R)) {
    writeLines(paste("File already exists; not saving:", fn_R))
} else {
    save(file=fn_R, "dge_all", "dge", "design", "voomOutput", "vfit", "efit", "summaryTable")
    print("Saved various outputs to Rdata object.")
}


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

# # If you only have one contrast and you're leaving the intercept, the stuff you want will be in the 2nd columns of things
# idx = 2

for (idx in 1:dim(contrastMat)[2])

    {
        fn_ranks_this = paste(respath, colnames(contrastMat)[idx], "_limma_ranked_genes.csv", sep="")
        if (file_checks(fn_ranks_this, shouldExist=FALSE, verbose=TRUE))
            {
            summaryTable = topTable(efit, n=Inf, sort.by='none', coef = idx )
            # Get ranks; sort the gene symbols, LFCs, and p-values accordingly
            newOrder = order(summaryTable$adj.P.Val)
            gene = rownames(summaryTable)[newOrder]
            lfc = efit$coefficients[newOrder, idx]
            pval = summaryTable$adj.P.Val[newOrder]

            # Get values for "Direction" column
            dn_idx = which(lfc<0)
            up_idx = which(lfc>0)
            direction = lfc
            direction[up_idx] = paste('up', nonRefName, sep='')
            direction[dn_idx] = paste('dn', nonRefName, sep='')

            # Going to separate down-regulated from up-regulated genes
            new_idx = c(dn_idx, up_idx)
            
            # Get ranks too, even though they're comparing each gene to all the others, not just those changing in the same direction
            ranks = 1:dim(summaryTable)[1]
            
            # Put it all in a data frame
            # I STILL WANT TO TRIPLE-CHECK IF PVAL_RANK[NEW_IDX] IS CORRECT
            genes_by_rank = data.frame(direction[new_idx],  gene[new_idx], lfc[new_idx], ranks[new_idx], pval[new_idx])
            colnames = c('Direction', 'Gene symbol', 'LFC', 'Rank', 'p-value')
            
            writeLines(paste("Writing ranked gene lists to", fn_ranks_this))
            write.table(genes_by_rank, quote=FALSE, row.names=FALSE, col.names=colnames, file=fn_ranks_this, sep=',')
            
            }
        
    }
