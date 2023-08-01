library(tidyverse)
library(Taxonstand)

### read taxonomy CDR
taxCDR <- read.csv("Taxonomy/cdcr-out.csv", na.strings = c("", "NA", " "))

## clean names
taxCDR <- taxCDR %>% 
  drop_na(spec.name.new) %>% 
  filter(specificEpithet != "sp.") %>% 
  distinct(spec.name.new, .keep_all = TRUE)

## search taxonomy using TPL 
spptaxCRD <- TPL(taxCDR$spec.name.new, infra = TRUE, corr = TRUE) # it will return a lot of warnings, please do not pay attention to that. 

### read CDR phylogeny
phyCDR <- ape::read.tree("CDR_phylogeny_v3/4.Timetree10/family-levelCalibration.tre")

## search taxonomy using TPL 
sppPhy <- gsub("_", " ", phyCDR$tip.label)

sppphyCRD <- TPL(sppPhy, infra = TRUE, corr = TRUE) # it will return a lot of warnings, please do not pay attention to that. 

### Join taxCDR and phyCDR
# remove ferns
taxrem <- c("Zygnema circumcarinatum", "Osmundastrum cinnamomeum", "Stipa comata", 
            "Equisetum arvense", "Equisetum hyemale", "Equisetum laevigatum", 
            "Gymnocarpium dryopteris", "Athyrium filix-femina", "Dryopteris cristata", 
            "Dryopteris carthusiana", "Onoclea sensibilis", "Osmundastrum cinnamomeum", 
            "Osmunda claytoniana", "Osmunda regalis", "Pteridium aquilinum", 
            "Claytosmunda claytoniana")

## CDR taxonomy
spptaxCRD <- spptaxCRD %>% 
  rename(scinameCDR = Taxon) %>% 
  mutate(scinameTPL = paste0(New.Genus, "_", New.Species)) %>% 
  select(scinameCDR, scinameTPL)

# Remove ferns from the CDR taxonomy 
spptaxCRD_noferns <- spptaxCRD %>% 
  filter(scinameCDR != taxrem[1] & scinameCDR != taxrem[2] & scinameCDR != taxrem[3] & 
           scinameCDR != taxrem[4] & scinameCDR != taxrem[5] & scinameCDR != taxrem[6] & 
           scinameCDR != taxrem[7] & scinameCDR != taxrem[8] &  scinameCDR != taxrem[9] & 
           scinameCDR != taxrem[10] & scinameCDR != taxrem[11] & scinameCDR != taxrem[12] & 
           scinameCDR != taxrem[13] & scinameCDR != taxrem[14] & scinameCDR != taxrem[15] & 
           scinameCDR != taxrem[16])

## taxonomy phylogeny CDR
sppphyCRD <- sppphyCRD %>% 
  rename(scinamePHY = Taxon) %>% 
  mutate(scinameTPL = paste0(New.Genus, "_", New.Species)) %>% 
  select(scinamePHY, scinameTPL)

## combine both set of species to make a taxon lookup
cdr_lookup <- full_join(spptaxCRD_noferns, sppphyCRD, by = "scinameTPL")

write_csv(cdr_lookup, 
          file = "CDR_phylogeny_v3/5.CkeckTaxonomy/CDR_PHY_species.csv") 

##### Prepare taxonomic data for species imputation #####
library(tidyverse)

tax <- read_csv("Dropbox/Collaborations/CDR/Taxonomy/cdcr-out.csv") 

tax <- tax %>% 
  filter(spec.name.new != taxrem[1] & spec.name.new != taxrem[2] & spec.name.new != taxrem[3] & 
           spec.name.new != taxrem[4] & spec.name.new != taxrem[5] & spec.name.new != taxrem[6] & 
           spec.name.new != taxrem[7] & spec.name.new != taxrem[8] &  spec.name.new != taxrem[9] & 
           spec.name.new != taxrem[10] & spec.name.new != taxrem[11] & spec.name.new != taxrem[12] & 
           spec.name.new != taxrem[13] & spec.name.new != taxrem[14] & spec.name.new != taxrem[15] & 
           spec.name.new != taxrem[16])

tax <- tax %>% 
  drop_na(spec.name.new) %>% 
  filter(family != "Polypodiaceae") %>% 
  filter(specificEpithet != "sp.") %>% 
  distinct(spec.name.new, .keep_all = TRUE)

tax <- tax %>% 
  select(family, genus, spec.name.new)

phy <- read_csv("Dropbox/Collaborations/CDR/CDR_phylogeny_v3/1.CuratedSequences/1.TaxonomyOriginal.csv")

phy <- phy %>% 
  distinct(species, .keep_all = TRUE) %>% 
  select(class, order, family, genus, species)

tax_TACT <- right_join(phy, tax,  
          by = c("family", "genus", "species" = "spec.name.new")) 

tax_TACT <- tax_TACT %>% 
  distinct(species, .keep_all = TRUE)

write_csv(tax_TACT, "Documents/tact_CDR/CDR_taxonomy_v2.csv")

##### Edit tree #####
library(ape)

treeCDR <- read.tree("Documents/tact_CDR/CDR_phylo.tre")

taxon <- taxoCDR %>% 
  select(scinamePHY, scinameFINAL, Notes) %>% 
  filter(scinamePHY != "NA") %>% 
  mutate(tipLabel = gsub(" ", "_", scinameFINAL))


