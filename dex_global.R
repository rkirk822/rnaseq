
# dex_global.R
#
# Make some plots to get a sense of global changes in expression between conditions.
#
# Line plot showing number of significantly DEX genes as alpha decreases.
# One line for up-regulation, one for down.
# A second line plot that's the same but shows median LFC, in each case.
# Not sure that one's so useful.

library(plotly)
source("/Users/nelsonlab/Documents/Toolboxes/rna-seq/file_checks.R")

#########################
# Plot stuff
#########################

# Will look at genes significantly DEX with p < each of these values, and LFC > lfcMin
alphas = seq(from=0.05, to=0.01, by=-0.01)
lfcMin = 1.5

# Marker size
mSize = 5

#########################
# Data stuff
#########################
projectPath =  "/Users/nelsonlab/Documents/Results_temporarily_here/TTX_results/"
# projectPath =  "/Users/nelsonlab/Documents/Results_temporarily_here/TTX_results/Using_only_current_comparison_for_modelling/"
comparison = "PV_Early"
treatment = "TTX"
# Read file with pvals and LFC into a dataframe.
# This file is created by limma_voom.R and includes limma-adjusted p-values and LFC
# for all genes.  (The LFC is the limma value too, not directly from TPM.)
dfLimma = read.csv(paste(projectPath, comparison, "_limma_ranked_genes.csv", sep = ""))

# Get the values for the two lines to plot
nUp = array(dim=length(alphas))
nDn = array(dim=length(alphas))
lfcUp = vector("list", length(alphas))
lfcDn = vector("list", length(alphas))
for (i in 1:length(alphas)) {
    # Genes with p.value < alpha and Direction is "up"; then again for "dn"
    upIdx = which(dfLimma$p.value<alphas[i] & dfLimma$Direction==paste("up", treatment, sep=""))
    dnIdx = which(dfLimma$p.value<alphas[i] & dfLimma$Direction==paste("dn", treatment, sep=""))
    # Restrict to those that meet LFC threshold
    upIdx = intersect(upIdx, which(dfLimma$LFC > lfcMin))
    dnIdx = intersect(dnIdx, which(abs(dfLimma$LFC) > lfcMin))
    # Number of those genes
    nUp[i] = length(upIdx)
    nDn[i] = length(dnIdx)
    # Their LFCs (all of em)
    lfcUp[[i]] = dfLimma$LFC[upIdx]
    lfcDn[[i]] = dfLimma$LFC[dnIdx]
}


# Write stuff to csv
df = data.frame(nUp, nDn, abs(sapply(lfcUp, median)), abs(sapply(lfcDn, median)))
rownames(df) = alphas
colnames(df) = c('nUp', 'nDn', 'Median abs(LFC) up', 'Median abs(LFC) dn')
fn = paste(projectPath, comparison, "_dex_global.csv", sep="")
if (file_checks (fn, shouldExist=FALSE, verbose=TRUE) ) { write.csv(df, file = fn) }


# Plot
pCounts = plot_ly(x = alphas) %>%
    add_trace(y = nUp, name = "Up-regulated",
        mode="lines+markers", type="scatter", marker=list(size=mSize)) %>%
    add_trace(y = nDn, name = "Down-regulated",
        mode="lines+markers", type="scatter", marker=list(size=mSize)) %>%
    layout(title=paste(comparison, "- Up/down-regulation in condition", treatment),
        xaxis=list(autorange="reversed", title="alpha"),
        yaxis = list(title="DEX gene count"))

htmlwidgets::saveWidget(as_widget(pCounts), paste(projectPath, comparison, "_sig_counts.html", sep=""))

# Plot their median LFCs
pLFCs = plot_ly(x = alphas) %>%
    add_trace(y = sapply(lfcUp, median), name = "Up-regulated",
        mode="lines+markers", type="scatter", marker=list(size=mSize)) %>%
    add_trace(y = abs(sapply(lfcDn, median)), name = "Down-regulated",
        mode="lines+markers", type="scatter", marker=list(size=mSize)) %>%
    layout(title=paste(comparison, "- Up/down-regulation effect size"),
        xaxis=list(autorange="reversed", title="alpha"),
        yaxis = list(title="Median abs(LFC)"))

htmlwidgets::saveWidget(as_widget(pLFCs), paste(projectPath, comparison, "_median_LFCs.html", sep=""))


