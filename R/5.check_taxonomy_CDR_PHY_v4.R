##### Check taxonomy in new phylogenetic hypothesis #####
library(data.table)
library(WorldFlora)
library(tidyverse)

### Inferred CDR phylogeny
cdr_phy <- ape::read.nexus("4.TimeTree/CDR_timeTree_secondaryCalibrations.nex")

### Taxonomy used in phylogenetic inference
tax_phy <- read_csv("1.CuratedSequences/1.Taxonomy.csv") %>% 
  rename("sciname_IN_Phylogeny" = species)

## Check for missing species
setdiff(cdr_phy$tip.label, tax_phy$species_names)
# character(0)

##### Combine taxonomy of WFO and phylogeny #####

#### Load WFO ####
WFO.remember(WFO.file = "../WorldFlora/WFO_Backbone/classification.csv")

### Species in phylogeny
spp_IN_phylo <- data.frame(sciname = cdr_phy$tip.label) %>% 
  mutate(speciesPHY = gsub("_", " ", cdr_phy$tip.label))

### check taxonomy
sp_IN_WFO <- data.table(WFO.match(WFO.data = WFO.data, # World Flora data
                                  spec.data = spp_IN_phylo$speciesPHY, # CDR phylogeny
                                  counter = 1, 
                                  verbose = TRUE))

sp_IN_WFO

### this function pares down each input species with multiple matches to the 'accepted' name or what WF determines is the closest match.
### this is a step for more scrutiny if you get odd output on the other end.
sp_IN_WFO <- WFO.one(sp_IN_WFO, 
                     priority = "Accepted")

## Check for missing species
setdiff(sp_IN_WFO$spec.name.ORIG, tax_phy$sciname_IN_Phylogeny)
# character(0)

# no match
sp_IN_WFO[Matched == FALSE, .(spec.name.ORIG)] 

## data.frame for output cleaned up a bit
sp_CDR_WFO <- sp_IN_WFO %>% 
  select(sciname_IN_Phylogeny = spec.name.ORIG, Matched, Fuzzy, scientificName, 
         taxonRank, majorGroup, family, genus, specificEpithet, Old.status) %>% 
  mutate(genus = if_else(taxonRank == 'genus', scientificName, genus), 
         specificEpithet = if_else(taxonRank == 'genus', 'sp.', specificEpithet)) %>% 
  mutate(sciname_IN_WFO = paste(genus, specificEpithet))

# see results
sp_CDR_WFO

write_csv(sp_CDR_WFO, 
          file = "5.Taxonomy/lookup_IN_phylogeny.csv")

### Taxonomy checked across the WFO
tax_wfo <- read_csv("5.Taxonomy/lookup_IN_phylogeny.csv")

### Combine both datasets 

tax_PHY_WFO <- full_join(x = tax_phy, 
          y = tax_wfo, 
          by = c("sciname_IN_Phylogeny")
          )

tax_PHY_WFO <- tax_PHY_WFO %>% 
  select(kingdom, phylum, majorGroup, class, order, 
         family_PHY = family.x, family_WFO = family.y, 
         genus_PHY = genus.x, genus_WFO = genus.y, specificEpithet, 
         sciname_IN_WFO, sciname_IN_Phylogeny
  ) %>% 
  mutate(tipLabel = gsub(" ", "_", sciname_IN_Phylogeny), 
         IN_Phylogeny = 1)

## Check names
setdiff(tax_PHY_WFO$sciname_IN_Phylogeny, tax_PHY_WFO$sciname_IN_WFO)
setdiff(tax_PHY_WFO$sciname_IN_WFO, tax_PHY_WFO$sciname_IN_Phylogeny)

## Verify that no missing species is present in lookup and phylogeny
setdiff(tax_PHY_WFO$tipLabel, cdr_phy$tip.label)

write_csv(tax_PHY_WFO, 
          file = "5.Taxonomy/lookup_IN_Phylogeny_WFO.csv")

##### Load original taxonomy dataset and get taxonomic hierarchy #####

### read taxonomy in CDR
tax_CDR <- read.csv("Taxonomy/cdcr-out.csv", na.strings = c("", "NA", " ")) %>% 
  drop_na(spec.name.new)

### Get missing species in phylogeny
missingSPP <- unique(setdiff(tax_CDR$spec.name.new, tax_PHY_WFO$sciname_IN_WFO))

missingSPP

### Missing species identified by Vini and Maowei
missing_VM <- readxl::read_xlsx("Taxonomy/add_spp_to_phylo.xlsx") %>% 
  mutate(missing = gsub("_", " ", SPP))

missing_VM <- as.character(missing_VM$missing)

### Combine both sets of missing species 
missingSPP <- unique(c(missingSPP, missing_VM))

### Double-check taxonomy
sp_NOphy_IN_WFO <- data.table(WFO.match(WFO.data = WFO.data, # World Flora data
                                        spec.data = missingSPP, # CDR phylogeny
                                        counter = 1, 
                                        verbose = TRUE))

sp_NOphy_IN_WFO

### this function pares down each input species with multiple matches to the 'accepted' name or what WF determines is the closest match.
### this is a step for more scrutiny if you get odd output on the other end.
sp_NOphy_IN_WFO <- WFO.one(sp_NOphy_IN_WFO, 
                           priority = "Accepted")

sp_NOphy_IN_WFO <- sp_NOphy_IN_WFO %>% 
select(sciname_NO_Phylogeny = spec.name.ORIG, Matched, Fuzzy, scientificName, 
       taxonRank, majorGroup, family, genus, specificEpithet, Old.status) %>% 
  mutate(genus = if_else(taxonRank == 'genus', scientificName, genus), 
         specificEpithet = if_else(taxonRank == 'genus', 'sp.', specificEpithet)) %>% 
  mutate(sciname_IN_WFO = paste(genus, specificEpithet), 
         tipLabel = gsub(" ", "_", sciname_IN_WFO), 
         IN_Phylogeny = 0) %>% 
  rename(family_WFO = family, genus_WFO = genus) 

## Check for missing species
setdiff(sp_NOphy_IN_WFO$sciname_NO_Phylogeny, tax_PHY_WFO$sciname_IN_WFO)
setdiff(sp_NOphy_IN_WFO$sciname_NO_Phylogeny, tax_PHY_WFO$sciname_IN_Phylogeny)
setdiff(sp_NOphy_IN_WFO$tipLabel, tax_PHY_WFO$tipLabel)

write_csv(sp_NOphy_IN_WFO, 
          file = "5.Taxonomy/lookup_NO_phylogeny.csv")

##### Get higher ranks taxonomic information #####
source("../R/functions/taxonomy_retrieve.R")

sp_NOphy_IN_WFO <- read_csv("5.Taxonomy/lookup_NO_phylogeny.csv")

NOphy <- taxonomy.retrieve(species_names = as.character(sp_NOphy_IN_WFO$sciname_NO_Phylogeny))

NOphy <- NOphy %>% 
  rename(sciname_NO_Phylogeny = Used_names)

### Combine both datasets 

tax_NO_PHY_WFO <- full_join(x = sp_NOphy_IN_WFO, 
                            y = NOphy, 
                            by = c("sciname_NO_Phylogeny")
)

### Select columns and save results
tax_NO_PHY_WFO <- tax_NO_PHY_WFO %>% 
  select(kingdom, phylum, majorGroup, class, order, 
         family_PHY = family, family_WFO, 
         genus_PHY = genus, genus_WFO, specificEpithet, 
         sciname_IN_WFO, sciname_NO_Phylogeny, tipLabel, IN_Phylogeny)

write_csv(tax_NO_PHY_WFO, 
          file = "5.Taxonomy/lookup_NO_Phylogeny_WFO.csv")

########## ---------- Taxonomic backbone ---------- ##########
library(tidyverse)

##### Create taxonomic backbone to impute missing species #####

### Read taxonomy species in phylogeny 
tax_PHY_WFO <- read_csv("5.Taxonomy/lookup_IN_Phylogeny_WFO.csv")

setdiff(tax_PHY_WFO$tipLabel, cdr_phy$tip.label)

### Read taxonomy missing species or NOT in phylogeny 
tax_NO_PHY_WFO <- read_csv("5.Taxonomy/lookup_NO_Phylogeny_WFO.csv") %>% 
  filter(majorGroup %in% c("A", "G"))

setdiff(tax_NO_PHY_WFO$tipLabel, cdr_phy$tip.label)

##### Select relevant columns and combine data ##### 

### Taxonomy in phylogeny
tax_IN <- tax_PHY_WFO %>% 
  select(majorGroup, class, order, family_WFO, genus_PHY, sciname_IN_Phylogeny, IN_Phylogeny) %>% 
  rename(species = sciname_IN_Phylogeny)

### Taxonomy no in phylogeny
tax_NO <- tax_NO_PHY_WFO %>% 
  select(majorGroup, class, order, family_WFO, genus_PHY, sciname_NO_Phylogeny, IN_Phylogeny) %>% 
  rename(species = sciname_NO_Phylogeny)

CDR_backbone <- bind_rows(tax_IN, tax_NO) %>% 
  distinct(majorGroup, class, order, family_WFO, genus_PHY, species)

write_csv(CDR_backbone, 
          file = "5.Taxonomy/CDR_taxonomic_backbone.csv")

unique(setdiff(CDR_backbone$species, gsub("_", " ", cdr_phy$tip.label)))

CDR_backbone_dups <- bind_rows(tax_IN, tax_NO) %>% 
  select(!"IN_Phylogeny")

write_csv(CDR_backbone_dups, 
          file = "5.Taxonomy/CDR_taxonomic_backbone_dups.csv")

##### Save final taxonomic database #####

cdr_phy_tacted <- ape::read.nexus("4.1.TimeTree_tacted/CDR_timeTree_tacted.nexus.tre")

setdiff(cdr_phy_tacted$tip.label, cdr_phy$tip.label)

### Read taxonomy species in phylogeny 
tax_PHY_WFO <- read_csv("5.Taxonomy/lookup_IN_Phylogeny_WFO.csv") %>% 
  mutate(Imputed = "NO")

setdiff(tax_PHY_WFO$tipLabel, cdr_phy$tip.label)

### Read taxonomy missing species or NOT in phylogeny 
tax_NO_PHY_WFO <- read_csv("5.Taxonomy/lookup_NO_Phylogeny_WFO.csv") %>% 
  filter(majorGroup %in% c("A", "G")) %>% 
  mutate(Imputed = "YES") %>% 
  rename(sciname_IN_Phylogeny = sciname_NO_Phylogeny)

### Combine both taxonomic sets
tax_PHY_WFO_TACTED <- bind_rows(tax_PHY_WFO, tax_NO_PHY_WFO) %>% 
  distinct(tipLabel, .keep_all = TRUE)

setdiff(tax_PHY_WFO_TACTED$tipLabel, cdr_phy_tacted$tip.label)

write_csv(tax_PHY_WFO_TACTED, 
          file = "5.Taxonomy/lookup_PHYLOGENY_WFO_TACTED.csv")
