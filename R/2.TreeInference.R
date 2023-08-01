# https://docs.ropensci.org/phruta/articles/Phylogenetics_phruta.html

## Open RStudio from terminal using:
# open /Applications/RStudio.app

library(phruta)
library(geiger)

system("raxmlHPC")

exec <- "/raxmlHPC.exe"
##### Phylogenetic inference with phruta and RAxML #####

tree.raxml(folder = '2.Alignments', 
           FilePatterns = 'Masked_', 
           raxml_exec = 'raxmlHPC', 
           Bootstrap = 100, 
           outgroup = "Ginkgo_biloba", 
           partitioned = FALSE)

##### Tree dating in phruta #####
## remove duplicates from the taxonomy table
tax <- read_csv("1.CuratedSequences/1.TaxonomyOriginal.csv")

tax2 <- tax[duplicated(tax$species_names), ]

tax2 <- tax %>% 
  distinct(species_names, .keep_all = TRUE)

write_csv(tax2, file = "1.CuratedSequences/1.Taxonomy.csv")

## Using 100 bootstraps
tree.dating(taxonomyFolder = "1.CuratedSequences", 
            phylogenyFolder = "3.Phylogeny_raxml100", 
            scale = 'treePL')

##### Visualization #####
trCDR <- read.tree("4.Timetree/family-levelCalibration.tre")

write.nexus(trCDR, file = "4.Timetree/CDR_calibratedPhylo.nex")

plot(trCDR, show.tip.label = FALSE)

