# DADA2 pipeline - Quality 
# script by Joao Sabino
Sys.time() # for log file

# Clear R environment
# This prevents using objects created previously
rm(list=ls())

# Starting point (for log file)
sessionInfo()

#########################################################################################################
# Load packages
# 
#########################################################################################################

library(dada2)

###############################################################################################
# DADA2 pipeline
#
###############################################################################################

seqtabB3AL1 <- readRDS("/sc/orga/projects/MECONIUM/Batch3/dada2/stringent/Batch3A_Lane1/Batch3A_Lane1_seqtab.rds")
namesB3AL1 <- rownames(seqtabB3AL1)
namesB3AL1_short <- gsub("_3A_L1_FR", "", namesB3AL1)
seqtabB3AL1_2 <- seqtabB3AL1
rownames(seqtabB3AL1_2) <- namesB3AL1_short

seqtabB3AL2 <- readRDS("/sc/orga/projects/MECONIUM/Batch3/dada2/stringent/Batch3A_Lane2/Batch3A_Lane2_seqtab.rds")
namesB3AL2 <- rownames(seqtabB3AL2)
namesB3AL2_short <- gsub("_3A_L2_FR", "", namesB3AL2)
seqtabB3AL2_2 <- seqtabB3AL2
rownames(seqtabB3AL2_2) <- namesB3AL2_short

seqtabB3BL1 <- readRDS("/sc/orga/projects/MECONIUM/Batch3/dada2/stringent/Batch3B_Lane1/Batch3B_Lane1_seqtab.rds")
namesB3BL1 <- rownames(seqtabB3BL1)
namesB3BL1_short <- gsub("_3B_L1_FR", "", namesB3BL1)
seqtabB3BL1_2 <- seqtabB3BL1
rownames(seqtabB3BL1_2) <- namesB3BL1_short

seqtabB3BL2 <- readRDS("/sc/orga/projects/MECONIUM/Batch3/dada2/stringent/Batch3B_Lane2/Batch3B_Lane2_seqtab.rds")
namesB3BL2 <- rownames(seqtabB3BL2)
namesB3BL2_short <- gsub("_3B_L2_FR", "", namesB3BL2)
seqtabB3BL2_2 <- seqtabB3BL2
rownames(seqtabB3BL2_2) <- namesB3BL2_short


st.3A <- mergeSequenceTables(seqtabB3AL1_2, seqtabB3AL2_2, repeats="sum")
st.3B <- mergeSequenceTables(seqtabB3BL1_2, seqtabB3BL2_2, repeats="sum")
st.all <- mergeSequenceTables(st.3A, st.3B)

saveRDS(st.all, "/sc/orga/projects/MECONIUM/Batch3/dada2/stringent/seqtab_Batch3_stringent_all.rds")

# Remove chimeras
seqtab_nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
saveRDS(seqtab_nochim, "/sc/orga/projects/MECONIUM/Batch3/dada2/stringent/seqtab_Batch3_all_stringent_nochim.rds")

# Assign taxonomy
taxa <- assignTaxonomy(seqtab_nochim, "/sc/orga/projects/MECONIUM/bin/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(taxa, "/sc/orga/projects/MECONIUM/Batch3/dada2/stringent/tax_Batch3_stringent_final.rds")

taxa <- addSpecies(taxa, "/sc/orga/projects/MECONIUM/bin/silva_species_assignment_v132.fa.gz")
saveRDS(taxa, "/sc/orga/projects/MECONIUM/Batch3/dada2/stringent/tax_Batch3_final_stringent_species.rds")


###############################################################################################
# END 
###############################################################################################

sessionInfo()
Sys.time()