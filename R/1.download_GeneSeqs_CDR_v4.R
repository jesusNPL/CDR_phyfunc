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

spp <- gsub("_", " ", spp)

### add species not in the original list
tax_add <- readxl::read_excel("Taxonomy/add_spp_to_phylo.xlsx")

spp_add <- as.character(tax_add$SPP)

spp_add <- gsub("_", " ", spp_add)

### Combine species lists

spp_all <- c(spp, spp_add)

spp_all <- sort(unique(spp_all))

setdiff(spp_add, spp)


##### Genes to be used #####
genes <- c("rbcL", "matK", "ndhF", "atpB", "trnL", "trnK", "ITS")

##### Download genetic data from GENBANK #####

### RBCL - Done!
acc.table.rbcl <- acc.table.retrieve(
  clades  = NULL, #"Ginkgo", # outgroup
  species = sort(spp_add),
  genes = genes[1],
  speciesLevel = TRUE
)

write.csv(acc.table.rbcl, "0_Accessions/CDR_RCBL_accessions_spp_add.csv", row.names = FALSE) 

spp2 <- gsub("_", " ", spp)

acc_rbcl <- acc.table.rbcl %>% 
  distinct(Species) 

### MATK - Done!
acc.table.matk <- acc.table.retrieve(
  clades  = NULL, #"Ginkgo", # outgroup
  species = sort(spp_add),
  genes = genes[2],
  speciesLevel = TRUE
)

write.csv(acc.table.matk, "0_Accessions/CDR_MATK_accessions_spp_add.csv", row.names = FALSE) 

### NDHF - Done!
acc.table.ndhF <- acc.table.retrieve(
  clades  = NULL, #"Ginkgo", # outgroup
  species = spp_add,
  genes = genes[3],
  speciesLevel = TRUE
)

head(acc.table.ndhF)

write.csv(acc.table.ndhF, "0_Accessions/CDR_NDHF_accessions_spp_add.csv", row.names = FALSE) 

### ATPB 
acc.table.atpb <- acc.table.retrieve(
  clades  = NULL, #"Ginkgo", # outgroup
  species = sort(spp_add),
  genes = genes[4],
  speciesLevel = TRUE
)

write.csv(acc.table.atpb, "0_Accessions/CDR_ATPB_accessions_spp_add.csv", row.names = FALSE) 

### TRNL - Done!
acc.table.trnl <- acc.table.retrieve(
  clades  = NULL, #"Ginkgo", # outgroup
  species = sort(spp_add),
  genes = genes[5],
  speciesLevel = TRUE
)

write.csv(acc.table.trnl, "0_Accessions/CDR_TRNL_accessions_spp_add.csv", row.names = FALSE)

### TRNK -
acc.table.trnk <- acc.table.retrieve(
  clades  = NULL,#"Ginkgo", # outgroup
  species = sort(spp_add),
  genes = genes[6],
  speciesLevel = TRUE
)

write.csv(acc.table.trnk, "0_Accessions/CDR_TRNK_accessions.csv", row.names = FALSE)

### ITS
acc.table.its <- acc.table.retrieve(
  clades  = NULL, #"Ginkgo", # outgroup
  species = sort(spp_add),
  genes = "its",#genes[7],
  speciesLevel = TRUE
)

write.csv(acc.table.trnk, "0_Accessions/CDR_ITS_accessions_spp_add.csv", row.names = FALSE)

##### Download genetic data #####

library(phruta)
library(tidyverse)

### Read accession numbers 
## ATPB
cdr_atpb_acc <- read.csv("0_Accessions/CDR_ATPB_accessions.csv") 
cdr_atpb_acc_add <- read.csv("0_Accessions/CDR_ATPB_accessions_spp_add.csv") 

cdr_atpb_acc <- bind_rows(cdr_atpb_acc, cdr_atpb_acc_add)

## MATK
cdr_matk_acc <- read.csv("0_Accessions/CDR_MATK_accessions.csv") 
cdr_matk_acc_add <- read.csv("0_Accessions/CDR_MATK_accessions_spp_add.csv") 

cdr_matk_acc <- bind_rows(cdr_matk_acc, cdr_matk_acc_add)

## NDHF
cdr_ndhf_acc <- read.csv("0_Accessions/CDR_NDHF_accessions.csv") 
cdr_ndhf_acc_add <- read.csv("0_Accessions/CDR_NDHF_accessions_spp_add.csv") 

cdr_ndhf_acc <- bind_rows(cdr_ndhf_acc, cdr_ndhf_acc_add)

## RCBL
cdr_rcbl_acc <- read.csv("0_Accessions/CDR_RCBL_accessions.csv") 
cdr_rcbl_acc_add <- read.csv("0_Accessions/CDR_RCBL_accessions_spp_add.csv") 

cdr_rcbl_acc <- bind_rows(cdr_rcbl_acc, cdr_rcbl_acc_add)

## TRNL
cdr_trnl_acc <- read.csv("0_Accessions/CDR_TRNL_accessions.csv") 
cdr_trnl_acc_add <- read.csv("0_Accessions/CDR_TRNL_accessions_spp_add.csv") 

cdr_trnl_acc <- bind_rows(cdr_trnl_acc, cdr_trnl_acc_add)

## TRNK
cdr_trnk_acc <- read.csv("0_Accessions/CDR_TRNK_accessions.csv") 
cdr_trnk_acc_add <- read.csv("0_Accessions/CDR_TRNK_accessions_spp_add.csv") 

cdr_trnk_acc <- bind_rows(cdr_trnk_acc, cdr_trnk_acc_add)
#cdr_its_acc <- read.csv("DCR_phylogeny/CDR_ITS_accessions.csv") 

# Insert true outgroup
#cdr_rcbl_acc <- bind_rows(acc.table.rbcl, cdr_rcbl_acc)

## Remove this species
taxrem <- c("Zygnema circumcarinatum", "Osmundastrum cinnamomeum", "Stipa comata", 
            "Equisetum arvense", "Equisetum hyemale", "Equisetum laevigatum", 
            "Gymnocarpium dryopteris", "Athyrium filix-femina", "Dryopteris cristata", 
            "Dryopteris carthusiana", "Onoclea sensibilis", "Osmundastrum cinnamomeum", 
            "Osmunda claytoniana", "Osmunda regalis", "Pteridium aquilinum", 
            "Claytosmunda claytoniana")
## ATPB
cdr_atpb_acc <- cdr_atpb_acc %>% 
  filter(Species != taxrem[1] & Species != taxrem[2] & Species != taxrem[3] & 
           Species != taxrem[4] & Species != taxrem[5] & Species != taxrem[6] & 
           Species != taxrem[7] & Species != taxrem[8] &  Species != taxrem[9] & 
           Species != taxrem[10] & Species != taxrem[11] & Species != taxrem[12] & 
           Species != taxrem[13] & Species != taxrem[14] & Species != taxrem[15] & Species != taxrem[16])

## MATK
cdr_matk_acc <- cdr_matk_acc %>% 
  filter(Species != taxrem[1] & Species != taxrem[2] & Species != taxrem[3] & 
         Species != taxrem[4] & Species != taxrem[5] & Species != taxrem[6] & 
         Species != taxrem[7] & Species != taxrem[8] &  Species != taxrem[9] & 
         Species != taxrem[10] & Species != taxrem[11] & Species != taxrem[12] & 
         Species != taxrem[13] & Species != taxrem[14] & Species != taxrem[15] & Species != taxrem[16])

## NDHF
cdr_ndhf_acc <- cdr_ndhf_acc %>% 
  filter(Species != taxrem[1] & Species != taxrem[2] & Species != taxrem[3] & 
           Species != taxrem[4] & Species != taxrem[5] & Species != taxrem[6] & 
           Species != taxrem[7] & Species != taxrem[8] &  Species != taxrem[9] & 
           Species != taxrem[10] & Species != taxrem[11] & Species != taxrem[12] & 
           Species != taxrem[13] & Species != taxrem[14] & Species != taxrem[15] & Species != taxrem[16])

## RCBL
cdr_rcbl_acc <- cdr_rcbl_acc %>% 
  filter(Species != taxrem[1] & Species != taxrem[2] & Species != taxrem[3] & 
           Species != taxrem[4] & Species != taxrem[5] & Species != taxrem[6] & 
           Species != taxrem[7] & Species != taxrem[8] &  Species != taxrem[9] & 
           Species != taxrem[10] & Species != taxrem[11] & Species != taxrem[12] & 
           Species != taxrem[13] & Species != taxrem[14] & Species != taxrem[15] & Species != taxrem[16])

## TRNL
cdr_trnl_acc <- cdr_trnl_acc %>% 
  filter(Species != taxrem[1] & Species != taxrem[2] & Species != taxrem[3] & 
           Species != taxrem[4] & Species != taxrem[5] & Species != taxrem[6] & 
           Species != taxrem[7] & Species != taxrem[8] &  Species != taxrem[9] & 
           Species != taxrem[10] & Species != taxrem[11] & Species != taxrem[12] & 
           Species != taxrem[13] & Species != taxrem[14] & Species != taxrem[15] & Species != taxrem[16])

## TRNK
cdr_trnk_acc <- cdr_trnk_acc %>% 
  filter(Species != taxrem[1] & Species != taxrem[2] & Species != taxrem[3] & 
           Species != taxrem[4] & Species != taxrem[5] & Species != taxrem[6] & 
           Species != taxrem[7] & Species != taxrem[8] &  Species != taxrem[9] & 
           Species != taxrem[10] & Species != taxrem[11] & Species != taxrem[12] & 
           Species != taxrem[13] & Species != taxrem[14] & Species != taxrem[15] & Species != taxrem[16])


### Combine accession numbers for all markers 
cdr_accessions <- bind_rows(cdr_atpb_acc, cdr_matk_acc, cdr_ndhf_acc, cdr_rcbl_acc, 
                            cdr_trnl_acc, cdr_trnk_acc)#, cdr_its_acc)

### Save combined accesions
write_csv(cdr_accessions, 
          file = "0_Accessions/CDR_gen_accesions_2024.csv")

cdr_accessions %>% 
  count(Species)

##### Download genes #####
sq.retrieve.indirect(acc.table = cdr_accessions, 
                     download.sqs = TRUE)

### Curate sequences
library(msa)

### No removing outliers
sq.curate(filterTaxonomicCriteria = 'Ginkgo|Poa|Pinus|Quercus',
          mergeGeneFiles = NULL, 
          database = "gbif", 
          kingdom = "Plantae", 
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
sq.aln(folder = '1.CuratedSequences_withOutliers', 
       FilePatterns = "renamed", 
       mask = TRUE) # remove outliers in folder "OutliersRemoved"

### Including the outgroup
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

