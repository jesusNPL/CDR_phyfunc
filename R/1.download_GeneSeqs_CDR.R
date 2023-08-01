### Libraries
library(phruta)
library(tidyverse)

##### Prepare taxonomy #####

### read taxonomy 
tax <- read.csv("Taxonomy/cdcr-out.csv", na.strings = c("", "NA", " "))

tax <- tax %>% 
  drop_na(spec.name.new)

clades <- tax %>% 
  select(scientificName, spec.name.new) %>% 
  separate(spec.name.new, into = c("genus", "species"), sep = " ")

species <- clades %>% 
  filter(species != "sp.") %>% 
  mutate(Scinames = paste0(genus, "_", species))

genus <- species %>% 
  select(genus) %>% 
  distinct()

scinames <- species %>% 
  select(Scinames) %>% 
  distinct()

spp <- as.character(scinames$Scinames)

##### Genes to be used #####
genes <- c("rbcL", "matK", "ndhF", "atpB", "trnL", "trnK", "ITS")

##### Download genetic data from GENBANK #####

### RBCL - Done!
acc.table.rbcl <- acc.table.retrieve(
  clades  = NULL, # outgroup
  species = "Zygnema_circumcarinatum",#, spp),
  genes = genes[1],
  speciesLevel = TRUE
)

write.csv(acc.table.rbcl, "CDR_phylogeny_v2/CDR_RCBL_accessions.csv", row.names = FALSE) 

spp2 <- gsub("_", " ", spp)

acc_rbcl <- acc.table.rbcl %>% 
  filter(Species %in% spp2)

acc_rbcl %>% 
  distinct(Species)

# Zygnema circumcarinatum is the outgroup
### MATK - Done!
acc.table.matk <- acc.table.retrieve(
  clades  = NULL, # outgroup
  species = "Zygnema_circumcarinatum", #spp),
  genes = genes[2],
  speciesLevel = TRUE
)

write.csv(acc.table.matk, "CDR_phylogeny_v2/CDR_MATK_accessions.csv", row.names = FALSE) 

### NDHF - Done!
acc.table.ndhF <- acc.table.retrieve(
  clades  = NULL, # outgroup
  species = "Zygnema_circumcarinatum", #spp),
  genes = genes[3],
  speciesLevel = TRUE
)

head(acc.table.ndhF)

write.csv(acc.table.ndhF, "CDR_phylogeny_v2/CDR_NDHF_accessions.csv", row.names = FALSE) 

### ATPB 
acc.table.atpb <- acc.table.retrieve(
  clades  = NULL, # outgroup
  species = "Zygnema_circumcarinatum", #spp),
  genes = genes[4],
  speciesLevel = TRUE
)

write.csv(acc.table.atpb, "CDR_phylogeny_v2/CDR_ATPB_accessions.csv", row.names = FALSE) 

### TRNL - Done!
acc.table.trnl <- acc.table.retrieve(
  clades  = NULL, # outgroup
  species = "Zygnema_circumcarinatum", #spp),
  genes = genes[5],
  speciesLevel = TRUE
)

write.csv(acc.table.trnl, "CDR_phylogeny_v2/CDR_TRNL_accessions.csv", row.names = FALSE)

### TRNK -
acc.table.trnk <- acc.table.retrieve(
  clades  = NULL, #"Zygnema circumcarinatum", outgroup
  species = "Zygnema_circumcarinatum", #spp),
  genes = genes[6],
  speciesLevel = TRUE
)

write.csv(acc.table.trnk, "CDR_phylogeny_v2/CDR_TRNK_accessions.csv", row.names = FALSE)

### ITS
acc.table.its <- acc.table.retrieve(
  clades  = NULL,#"Zygnema circumcarinatum", # outgroup
  species = "Zygnema_circumcarinatum", #spp),
  genes = "its",#genes[7],
  speciesLevel = TRUE
)

write.csv(acc.table.trnk, "CDR_phylogeny_v2/CDR_ITS_accessions.csv", row.names = FALSE)

##### Download genetic data #####
setwd("CDR_phylogeny_v3/")

library(phruta)
library(tidyverse)


### Read accession numbers 
cdr_atpb_acc <- read.csv("CDR_ATPB_accessions.csv") 
cdr_matk_acc <- read.csv("CDR_MATK_accessions.csv") 
cdr_ndhf_acc <- read.csv("CDR_NDHF_accessions.csv") 
cdr_rcbl_acc <- read.csv("CDR_RCBL_accessions.csv") 
cdr_trnl_acc <- read.csv("CDR_TRNL_accessions.csv") 
#cdr_trnk_acc <- read.csv("DCR_phylogeny/CDR_TRNK_accessions.csv") 
#cdr_its_acc <- read.csv("DCR_phylogeny/CDR_ITS_accessions.csv") 

# Insert true outgroup
#cdr_rcbl_acc <- bind_rows(acc.table.rbcl, cdr_rcbl_acc)

# Remove these species
taxrem <- c("Zygnema circumcarinatum", "Osmundastrum cinnamomeum", "Stipa comata", 
            "Equisetum arvense", "Equisetum hyemale", "Equisetum laevigatum", 
            "Gymnocarpium dryopteris", "Athyrium filix-femina", "Dryopteris cristata", 
            "Dryopteris carthusiana", "Onoclea sensibilis", "Osmundastrum cinnamomeum", 
            "Osmunda claytoniana", "Osmunda regalis", "Pteridium aquilinum", 
            "Claytosmunda claytoniana")
## ATPB

keep_atpb <- setdiff(cdr_atpb_acc$Species, taxrem)

cdr_atpb_acc <- cdr_atpb_acc %>% 
  filter(Species != taxrem[1] & Species != taxrem[2] & Species != taxrem[3] & 
           Species != taxrem[4] & Species != taxrem[5] & Species != taxrem[6] & 
           Species != taxrem[7] & Species != taxrem[8] &  Species != taxrem[9] & 
           Species != taxrem[10] & Species != taxrem[11] & Species != taxrem[12] & 
           Species != taxrem[13] & Species != taxrem[14] & Species != taxrem[15] & Species != taxrem[16])

cdr_matk_acc <- cdr_matk_acc %>% 
  filter(Species != taxrem[1] & Species != taxrem[2] & Species != taxrem[3] & 
         Species != taxrem[4] & Species != taxrem[5] & Species != taxrem[6] & 
         Species != taxrem[7] & Species != taxrem[8] &  Species != taxrem[9] & 
         Species != taxrem[10] & Species != taxrem[11] & Species != taxrem[12] & 
         Species != taxrem[13] & Species != taxrem[14] & Species != taxrem[15] & Species != taxrem[16])

cdr_ndhf_acc <- cdr_ndhf_acc %>% 
  filter(Species != taxrem[1] & Species != taxrem[2] & Species != taxrem[3] & 
           Species != taxrem[4] & Species != taxrem[5] & Species != taxrem[6] & 
           Species != taxrem[7] & Species != taxrem[8] &  Species != taxrem[9] & 
           Species != taxrem[10] & Species != taxrem[11] & Species != taxrem[12] & 
           Species != taxrem[13] & Species != taxrem[14] & Species != taxrem[15] & Species != taxrem[16])

cdr_rcbl_acc <- cdr_rcbl_acc %>% 
  filter(Species != taxrem[1] & Species != taxrem[2] & Species != taxrem[3] & 
           Species != taxrem[4] & Species != taxrem[5] & Species != taxrem[6] & 
           Species != taxrem[7] & Species != taxrem[8] &  Species != taxrem[9] & 
           Species != taxrem[10] & Species != taxrem[11] & Species != taxrem[12] & 
           Species != taxrem[13] & Species != taxrem[14] & Species != taxrem[15] & Species != taxrem[16])

cdr_trnl_acc <- cdr_trnl_acc %>% 
  filter(Species != taxrem[1] & Species != taxrem[2] & Species != taxrem[3] & 
           Species != taxrem[4] & Species != taxrem[5] & Species != taxrem[6] & 
           Species != taxrem[7] & Species != taxrem[8] &  Species != taxrem[9] & 
           Species != taxrem[10] & Species != taxrem[11] & Species != taxrem[12] & 
           Species != taxrem[13] & Species != taxrem[14] & Species != taxrem[15] & Species != taxrem[16])

### Combine accession numbers for all markers 
cdr_accessions <- bind_rows(cdr_atpb_acc, cdr_matk_acc, cdr_ndhf_acc, cdr_rcbl_acc, 
                            cdr_trnl_acc)#, cdr_trnk_acc, cdr_its_acc)

write_csv(cdr_accessions, 
          file = "CDR_gen_accesions.csv")

##### Download genes #####
sq.retrieve.indirect(acc.table = cdr_accessions, 
                     download.sqs = TRUE)

### Curate sequences
library(msa)

### No removing outliers
sq.curate(filterTaxonomicCriteria = 'Ginkgo|Poa|Pinus|Quercus',
          mergeGeneFiles = NULL, 
          #database = "ncbi", 
          kingdom = NULL, #"Plantae", 
          folder = '0.Sequences',
          removeOutliers = FALSE)

### remove outliers
sq.curate(filterTaxonomicCriteria = 'Ginkgo|Poa|Pinus|Quercus',
          mergeGeneFiles = NULL,
          database = "gbif", 
          kingdom = "Plantae", 
          folder = '0.Sequences',
          removeOutliers = TRUE) # remove outliers in folder "OutliersRemoved"

##### Sequence alignment #####
### Masked
sq.aln(folder = '1.CuratedSequences', 
       FilePatterns = "renamed", 
       mask = TRUE) # remove outliers in folder "OutliersRemoved"

##### PartitionFinder #####

# This will create the partitions necessary to obtain the evolutionary models 
# https://github.com/brettc/partitionfinder/archive/v1.1.1.zip

sq.partitionfinderv1(folderAlignments = "2.Alignments",
                     FilePatterns = "Masked_", 
                     #folderPartitionFinder = "2.1.PartitionFinderv1", 
                     models = "all", 
                     run = TRUE
)

