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

batchname <- "Batch3A_Lane1"
demuxDir <- "/sc/orga/projects/MECONIUM/Batch3/preprocessing/Batch3A_Lane1"
mainDir <- "/sc/orga/projects/MECONIUM/Batch3/dada2/stringent/Batch3A_Lane1"



###############################################################################################
# DADA2 pipeline
#
###############################################################################################

filtDir <- "filtered"
dir.create(file.path(mainDir, filtDir))

fnFs <- sort(list.files(demuxDir, pattern=".R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(demuxDir, pattern=".R2.fastq", full.names = TRUE))
if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match.")

sample.names <- sapply(strsplit(basename(fnFs), ".R1"), `[`, 1)
filtFs <- file.path(mainDir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(mainDir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))


# Filtering
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE, 
                     compress=TRUE, multithread=FALSE) 

head(out)
write.table(out, paste0(batchname, "_filtered.txt"), sep="\t")
saveRDS(out, paste0(batchname, "_filtered.rds"))

# File parsing
filtpathF <- file.path(mainDir, "filtered") 
filtpathR <- file.path(mainDir, "filtered")
filtFs <- list.files(filtpathF, pattern="_F_filt.fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="_R_filt.fastq.gz", full.names = TRUE)
# Extract sample names
sample.namesR <- sapply(strsplit(basename(filtRs), "_R_"), `[`, 1)
sample.names <- sapply(strsplit(basename(filtFs), "_F_"), `[`, 1)
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR)
    mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

saveRDS(mergers, paste0(batchname, "_mergers.rds"))
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, paste0(batchname, "_seqtab.rds"))

# Remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
saveRDS(seqtab_nochim, paste0(batchname, "_seqtab_nochim.rds"))

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(mergers, getN))
colnames(track) <- c("input", "filtered", "mergers")
rownames(track) <- sample.names
head(track)
write.table(track, paste0(batchname, "_track.txt"), sep="\t")


###############################################################################################
# END 
###############################################################################################

sessionInfo()
Sys.time()
