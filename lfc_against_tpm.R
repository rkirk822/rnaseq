


# lfc_against_tpm.R
#
# Scatterplot for one sample.  One datapoint per gene, showing LFC against TPM.
# Color dots for genes with p<alpha red.

library(plotly)
source("/Users/nelsonlab/Documents/Toolboxes/rna-seq/file_checks.R")

#########################
# Plot stuff
#########################

# Will look at genes significantly DEX with p < each of these values
alpha = 0.01

# Marker size
mSize = 5

#########################
# Data stuff
#########################
projectPath =  "/Users/nelsonlab/Documents/Results_temporarily_here/TTX_results/"
comparison = "EMXEarly"
# Read file with TPM.
tpm = read.csv(paste(projectPath, comparison, "_TPM.csv", sep=""))
rownames(tpm) = tpm[,1]; tpm = tpm[,-1]  # Move the gene symbols into row names
ltpm = log2(tpm) # get log2(tpm)
ltpm[tpm==0] = 0  # make 0s be 0 instead of -Inf
# Read file with pvals and LFC into a dataframe.
# This file is created by limma_voom.R and includes limma-adjusted p-values and LFC
# for all genes.  (The LFC is the limma value too, not directly from TPM.)
dfLimma = read.csv(paste(projectPath, comparison, "_limma_ranked_genes.csv", sep = ""))


# For each sample
for (sIdx in 1:length(colnames(tpm))) {
    
    # Get TPM and LFC (limma stuff doesn't include all genes and is in different order; match them)
    sTPM = ltpm[ match(dfLimma$Gene.symbol,rownames(ltpm)), sIdx ]
    LFC = dfLimma$LFC
    
    # Factor variable indicating which are significant, to color datapoints
    sig = dfLimma$p.value < alpha
    sig[which(sig)] = paste("p < ", alpha, sep="")
    sig[which(sig=="FALSE")] = "Not significant"
    sig = as.factor(sig)
    # Make scatterplot
    p = plot_ly(x = sTPM, y=LFC, type="scatter", mode="markers", hoverinfo = "text", text=paste(dfLimma$Gene.symbol, "</br></br>LFC:", round(dfLimma$LFC,digits=2), "</br>LTPM:",
            round(sTPM,digits=2)),
        color=sig, colors=c("blue", "red")) %>%
        layout(title=colnames(tpm)[sIdx], xaxis=list(title="log2(TPM)"), yaxis=list(title="log2(FC)"))
    
    
    # Save as html file
    resFile = paste(projectPath, "LFC_TPM_", colnames(tpm)[sIdx], ".html", sep="")
    if (file_checks(resFile, shouldExist=FALSE, verbose=TRUE)) {
        writeLines(paste("Creating html file for sample", colnames(tpm)[sIdx]))
        htmlwidgets::saveWidget(as_widget(p), resFile)
    } else { next }
    
}


