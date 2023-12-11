# This R script was originally written by Peter Wilfahrt (https://scholar.google.com/citations?user=V1GzZ0MAAAAJ&hl=en)
# Jesus made some updates

##### ---------- CDR-LTR taxonomic resolution ---------- #####
#install.packages('WorldFlora')
library(data.table)
library(WorldFlora)
library(tidyverse)

##### Download taxonomic backbone #####

### The first time you use World Flora, or when you want to update to a new version, 
## you need to run this code. This is ~900 mb file that will be downloaded to the specified folder.
# To download the data you have two options: 
# Option 1) download the data directly using this link "https://www.worldfloraonline.org/downloadData"
# Option 2) using the function "WF.download"

WFO.download(WFO.file = "WorldFlora/WFO_Backbone/classification.txt")

##### Load data #####

### After downloading the first time, run this every session to load into environment
# Note that I downloaded the taxonomic backbone using the option 1. 
# Also, I'm uploading the backbone to the data folder

WFO.remember(WFO.file = "WorldFlora/WFO_Backbone/classification.csv")

### load in species data for checking - CDR plant species
sp.in <- read.table("CDR_phyfunc/taxonomy/Cedar Creek Plant Taxon List.txt", 
                    header = TRUE, 
                    sep = "\t")

sp.in 

##### Check Taxonomy #####

### spec.data should be the species column from 'sp.in' - this should be 'Genus species' formulated
# This part will take some time to finish; be patient!

sp.out.WF <- data.table(WFO.match(WFO.data = WFO.data, # World Flora data
                                  spec.data = sp.in$Species, # CDR plant species
                                  counter = 1, 
                                  verbose = TRUE))

sp.out.WF

### this function pares down each input species with multiple matches to the 'accepted' name or what WF determines is the closest match.
### this is a step for more scrutiny if you get odd output on the other end.
sp.out <- WFO.one(sp.out.WF)

# no match
sp.out[Matched == FALSE, .(spec.name.ORIG)] 

# only two identified taxa here: Polytricum sp. - misspelling, should be Polystrichum sp. (google is better at fuzzy matches when WF fails)
# the other is Botrichium sp. - it should be Botrychium sp.
# the rest of these are unidentified or non-live material

## data.frame for output cleaned up a bit
sp.CDR_WF <- sp.out %>% 
  select(spec.name.ORIG, Matched, Fuzzy, scientificName, taxonRank, family, genus, specificEpithet, Old.status) %>% 
  mutate(genus = if_else(taxonRank == 'genus', scientificName, genus), 
         specificEpithet = if_else(taxonRank == 'genus', 'sp.', specificEpithet)) %>% 
  mutate(spec.name.new = paste(genus, specificEpithet))

sp.CDR_WF

### Save matched taxonomy between CDR-LTR and WFO 
# Note that the columns "spec.name.ORIG" and "spec.name.new" correspond to 
# the original species names as in object "sp.in" and the corrected names, respectively.

write_csv(sp.CDR_WF, "CDR_phyfunc/taxonomy/CDR-LTR_match_WFO.csv")
