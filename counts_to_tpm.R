
# counts_to_tpm.R
# 
# Given a gene (or whatever) by sample read count matrix and a length for each gene, return TPM matrix.
# Follows http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/


counts_to_tpm = function(countMat, geneLengths) {

    rpk = countMat / (geneLengths/1000)          # Reads per kilobase
    scalingFactors = colSums(rpk, na.rm=TRUE) / 10^6   # "Per million" scaling factor
    tpm = t( t(rpk) / scalingFactors)         # Transcripts per million
    
    return(tpm)

    }

# # From https://www.biostars.org/p/171766/
# # Somewhat different values than above; I'd like to work out why.
# # Calculate reads per base
# rpb = countMat / geneLens
# # Divide by that value for all genes in sample, to get "gene fraction"
# gf = t( t(rpb / colSums(rpb)) )
# # Make it per million
# tpm2 = gf * 10^6


