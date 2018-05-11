

# sample_check
# 
# Looking for possible bad samples.  Very rough bar plots of mean expression for each sample,
# and each sample's mean correlation with the others.



###################################################################
# Stuff you want to set:  exprMat, sampleNames, corrsPdf, meansPdf
###################################################################
source('/Volumes/CodingClub1/RNAseq/code/read_fcounts.R')
# source('/Volumes/DataStorage2/Emma/From_CodingClub1/RNAseq/code/read_fcounts.R')
fileList=list.files(pattern='fcounts')
fcounts = read_fcounts(fileList)
source('/Volumes/CodingClub1/RNAseq/code/counts_to_tpm.R')
# source('/Volumes/DataStorage2/Emma/From_CodingClub1/RNAseq/code/counts_to_tpm.R')
exprMat = counts_to_tpm(fcounts$Count, fcounts$Length)
sampleNames = sub('_fcounts.txt', '', fileList)
# Check that units you're using are in the filename (counts, TPM, etc)
corrsPdf = '/Users/nelsonlab/Documents/Toolboxes/limma_voom/TTX_results_fixing/Mean_corrs_tpm.pdf'
meansPdf = '/Users/nelsonlab/Documents/Toolboxes/limma_voom/TTX_results_fixing/Mean_expr_tpm.pdf'

####################################################################
# Bar plot of each sample's mean value, across all genes.
####################################################################
sampleMeans = colMeans(exprMat, na.rm = TRUE)
pdf(file = meansPdf)
barplot(sampleMeans, names.arg=sampleNames, cex.names=0.8, las=2)
dev.off()


####################################################################
# Bar plot of each sample's mean correlation with all other samples.
####################################################################
sampleCorrs = cor(exprMat, use='complete.obs')
diag(sampleCorrs) = NA  # don't include self-correlations
meanCorrs = rowMeans(sampleCorrs, na.rm=TRUE)
pdf(file = corrsPdf)
barplot(meanCorrs, names.arg=sampleNames, cex.names=0.8, las=2)
dev.off()
