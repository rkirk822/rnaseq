

# read_counts_reality_check.R
#
# This R script is to make sure that the read counts for RNAseq data look reasonable for selected
# genes (mostly cell-type markers).  The script produces a bar graph for each sample, showing
# selected genes' raw counts, and indicates the 90th percentile in the distribution of non-zero counts.
#
# In the "Stuff to set" section, tell the script how to locate and identify output files from the
# featureCounts tool, and what to name the pdf where the bar plots will go (all plots to one pdf).
#
# Note that percentile ranks are based on an *extremely* right-skewed distribution.  There are very few genes
# with very high counts.
#
# What I'd rather hae is just one bar plot showing mean across samples, with dots for individual samples, but
# I'm not going to figure out how to do that in R right now.

#############################################
# Stuff to set
#############################################
# countsDir, countsSuffix - to locate and identify the files with raw counts
countsDir='/Volumes/CodingClub1/RNAseq/TTX/counts/counts_m20_q20/'
# countsDir='/Users/nelsonlab/Documents/Smintheus_stuff_copied_here/counts/'
countsSuffix='_fcounts.txt'
sampleNames=c('EMXCtlEarly_1', 'EMXCtlEarly_2', 'EMXCtlEarly_3', 'EMXCtlEarly_4', 'EMXCtlLate_1', 'EMXCtlLate_2', 'EMXCtlLate_3', 'EMXCtlLate_4', 'EMXTTXEarly_1', 'EMXTTXEarly_2', 'EMXTTXEarly_3', 'EMXTTXLate_1', 'EMXTTXLate_2', 'EMXTTXLate_3', 'EMXTTXLate_4', 'PVCtlEarly_1', 'PVCtlEarly_2', 'PVCtlEarly_3', 'PVCtlLate_1', 'PVCtlLate_2', 'PVCtlLate_3', 'PVCtlLate_4', 'PVTTXEarly_1', 'PVTTXEarly_2', 'PVTTXEarly_3', 'PVTTXLate_1', 'PVTTXLate_2', 'PVTTXLate_3', 'PVTTXLate_4')
# sampleNames=c('BF_RORbHTp200_1', 'BF_RORbHTp200_2', 'BF_RORbHTp200_3', 'BF_RORbHTp200_4', 'BF_RORbHTp2_1', 'BF_RORbHTp2_2', 'BF_RORbHTp2_3', 'BF_RORbHTp2_4', 'BF_RORbHTp30_1', 'BF_RORbHTp30_2', 'BF_RORbHTp30_3', 'BF_RORbHTp30_4', 'BF_RORbHTp7_1', 'BF_RORbHTp7_2', 'BF_RORbHTp7_3', 'BF_RORbHTp7_4', 'BF_RORbKOp2_1', 'BF_RORbKOp2_2', 'BF_RORbKOp2_3', 'BF_RORbKOp2_4', 'BF_RORbKOp30_1', 'BF_RORbKOp30_2', 'BF_RORbKOp30_3', 'BF_RORbKOp30_4', 'BF_RORbKOp7_1', 'BF_RORbKOp7_2', 'BF_RORbKOp7_3', 'BF_RORbKOp7_4')
# # # For testing
# sampleNames=c('BF_RORbHTp2_1', 'BF_RORbHTp200_2')
# barplotPDF - Name for pdf where all plots will go
barplotPDF='/Volumes/CodingClub1/RNAseq/TTX/figures/TTX_fcounts_reality_check_m20_q20.pdf'
if(file.exists(barplotPDF)) {
    print('Pdf already exists; not gonna overwrite it.')
    stop()
}


#############################################
# Lists of genes expected to be high or low
#############################################
# Genes that should have weak/no expression
# These are from Erin's QC of the data, checking for glial contamination.
# Slcla2, Slcla3, S100a, and Metrn don't appear in the data.
low_genes=c('Aqp4', 'Gfap', 'Fgfr3', 'Gjb6','Vim','Mbp', 'Sox10', 'Mag', 'Mog')  # Slcla2 doens't appear
# Genes that should have strong expression
# From http://docs.abcam.com/pdf/neuroscience/Glutamatergic_cortical_neurogenesis.pdf:
# Bhlhb5 and Bm1 don't appear in the data.
high_genes=c('Cux1','Cux2','Lhx2','Mef2c')
# Check for no RORb in knockout samples
RORb='Rorb'



# Going to feed the gene lists into this function, to get their counts and percentile ranks.
calculate_ranks=function(set_names,all_names,all_counts) {
    # Get counts for gene subset
    set_counts=all_counts[match(set_names,all_names)]
    # Initialize ranks vector
    set_ranks=vector(,length(set_names))
    # Calculate ranks
    for(i in 1:length(set_names)) {
        set_ranks[i]=length(which(all_counts<set_counts[i]))/length(all_counts) * 100
    }
    # # Display counts and percentile ranks
    # # This would be best in a bar graph but I'm not taking the time right now
    # for(i in 1:length(set_names)) {
    # print(paste(set_names[i], ': Read count ', set_counts[i], '; ',
    # 	round(set_ranks[i],digits=2), 'th percentile', sep=''), quote=FALSE)
    # }
    # Return ranks, counts, and names of genes in set
    set_results = list("set_counts"=set_counts, "set_ranks"=set_ranks, "set_names"=set_names)
    return(set_results)
}


# Now go through each sample's featureCounts output, feeding to calculate_ranks function.
# And printing bar plot to this pdf:
pdf(file = barplotPDF)
for (s in 1:length(sampleNames)) {

    # Announce current file
    writeLines(paste('Looking at', sampleNames[s]))

    # Check that it exists
    filename=paste(countsDir,sampleNames[s],countsSuffix, sep='')
    if(!file.exists(filename)) {
        writeLines('This file does not exist; skipping:')
        writeLines(filename)
        next
    }

    # Get counts table for this file.
    counts_table=read.table(filename, header=TRUE)

    # Retrieve read counts for all genes, and their names.
    gnames=counts_table[,1]
    gcounts=counts_table[,7]
    
    # Limit to genes with nonzero counts, since those are all equal and will artificially boost ranks of our test genes
    gcounts_pos=gcounts[which(gcounts>0)]
    gnames_pos=gnames[which(gcounts>0)]
    
    # Calculate deciles so I can understand the distributions a little better
    decile_idx=seq(1,length(gcounts_pos),length.out=10)
    gcounts_pos_sorted=sort(gcounts_pos)
    deciles=gcounts_pos_sorted[round(decile_idx)]
    
    # Display percentage of nonzero reads, and deciles (okay, not displaying those here anymore)
    writeLines(paste('Percentage of genes with zero count: ', round(length(which(gcounts==0))/length(gcounts),digits=2)*100, '%', sep=''))
    # # Okay, it's not so helpful to display for every single sample.
    # writeLines('Deciles of read count distribution, excluding zero counts:')
    # for(d in 1:10) {
    #     writeLines(paste(d,'0th percentile: ',deciles[d],sep=''))
    # }
    

    # Give gcounts_pos and gnames_pos to calculate_ranks function.
    results_high=calculate_ranks(high_genes,gnames_pos,gcounts_pos)
    results_low=calculate_ranks(low_genes,gnames_pos,gcounts_pos)
    results_RORb=calculate_ranks(RORb,gnames_pos,gcounts_pos)
    
    # Make a bar graph showing counts and another showing ranks, each including all selected genes.
    all_set_ranks = c(results_high$set_ranks, results_low$set_ranks, results_RORb$set_ranks)
    all_set_counts = c(results_high$set_counts, results_low$set_counts, results_RORb$set_counts)
    all_set_names = c(results_high$set_names, results_low$set_names, results_RORb$set_names)
    # Colors for genes that should be high, should be low, and RORb
    barcolors=c(rep('red',times=length(high_genes)), rep('blue', times=length(low_genes)), 'green')
    # # Plot percentile ranks - not doing this now cuz not sure how informative it really is
    # barplot(all_set_ranks, names.arg=all_set_names, las=2,ylim=c(0, 100), # col=barcolors,main=paste(sampleNames[s], 'ranks'), ylab=' percentile ranks in distribution of all non-zero read counts for this sample')
    # legend location is based on assumption that low genes will in fact be low
    # legend(10,y=100,c('high','low','rorb'),fill=c('red','blue','green'))
    # Plot raw counts    
    barCounts=barplot(all_set_counts, names.arg=all_set_names, las=2, col=barcolors,main=paste(sampleNames[s],'counts'), ylab=' read counts in distribution of all non-zero read counts for this sample')
    text(length(high_genes)+length(low_genes),y=par('yaxp')[2],'Dashed line = 90th percentile')
    lines(x=barCounts,y=rep(deciles[9],times=length(barCounts)),lty='dashed')
    # text(length(high_genes)+length(low_genes),y=par('yaxp')[2],pos=1,decileStr)
    # legend location would need to be adjusted for this plot; not worth it now.

}

# Close connection to the file where you're putting bar plots, and create the pdf
dev.off()
