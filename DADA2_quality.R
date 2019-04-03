# DADA2 Big Data pipeline
# Reference: https://benjjneb.github.io/dada2/bigdata_paired.html
# script by Joao Sabino
Sys.time() # for log file

# Starting point (for log file)
sessionInfo()

#########################################################################################################
# Load packages
# 
#########################################################################################################

library(dada2)
sessionInfo()


###############################################################################################
# CHANGE ME!!!
#
###############################################################################################

batchname <- "Batch3B_Lane2"
demuxDir <- "/sc/orga/projects/MECONIUM/Batch3/preprocessing/Batch3B_Lane2"
mainDir <- "/sc/orga/projects/MECONIUM/Batch3/dada2/quality_profile"


###############################################################################################
# DADA2 pipeline
#
###############################################################################################

fnFs <- sort(list.files(demuxDir, pattern=".R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(demuxDir, pattern=".R2.fastq", full.names = TRUE))
if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match.")

# Inspect read quality profiles
# forward reads:
pdf(paste0(batchname, "_quality_profile_forward.pdf"),width =12, height=8,useDingbats=FALSE)      
plotQualityProfile(fnFs[1:4])
dev.off()

pdf(paste0(batchname, "_quality_profile_forward_aggregated.pdf"),width =12, height=8,useDingbats=FALSE)      
plotQualityProfile(fnFs, aggregate=TRUE)
dev.off()

# reverse reads:
pdf(paste0(batchname, "_quality_profile_reverse.pdf"),width =12, height=8,useDingbats=FALSE)      
plotQualityProfile(fnRs[1:4])
dev.off()

pdf(paste0(batchname, "_quality_profile_reverse_aggregated.pdf"),width =12, height=8,useDingbats=FALSE)      
plotQualityProfile(fnRs,aggregate=TRUE)
dev.off()


###############################################################################################
# END 
###############################################################################################

sessionInfo()
Sys.time()
