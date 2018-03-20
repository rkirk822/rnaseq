

# COUNT_SIG_GENES
#
# Read in files with ranked limma adjusted p-values (i.e. q-values).
# Pick a threshold, say 0.05, and get the number of significantly DEX genes
# for each comparison.  Separately for down- and up-regulation with TTX.
#
# Put in a check that you're not overwriting anything.

setwd('/Users/nelsonlab/Documents/Results_temporarily_here/TTX_results/Old/Using_only_current_comparison_for_modelling/limma_voom_output/')

# Desired FDR
thresh = 0.01
# Results filename
fn = 'TTX_qvals_under_pt01.csv'

# comparisons = c('PVEarly', 'PVLate', 'EMXEarly', 'EMXLate')
comparisons = c('PV_Early', 'PV_Late', 'EMX_Early', 'EMX_Late')
directions = c('dnTTX', 'upTTX')

counts = matrix(data = NA, nrow = length(comparisons), ncol = 2,)
rownames(counts) = comparisons
colnames(counts) = directions

for (c in 1:length(comparisons)) {
    
    # Read in file with q-values
    t = read.csv(paste(comparisons[c], 'limma_ranked_genes.csv', sep = '_'))
    
    # Identify down-regulated sig genes
    dnIdx = which(t$Direction==directions[1])
    dnQ = t$p.value[dnIdx]
    
    # Identify up-regulated sig genes
    upIdx = which(t$Direction==directions[2])
    upQ = t$p.value[upIdx]
    
    # Get counts
    counts[c,1] = length(which(dnQ < thresh))
    counts[c,2] = length(which(upQ < thresh))
}

write.csv(counts, file = fn)
