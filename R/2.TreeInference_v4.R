## Compile raxml manually
# https://www.metagenomics.wiki/tools/phylogenetic-tree/construction/raxml/install

# https://docs.ropensci.org/phruta/articles/Phylogenetics_phruta.html

## Open RStudio from terminal using:
# open /Applications/RStudio.app

##### Libraries #####
library(phruta)
library(geiger)

system("raxmlHPC")

##### Phylogenetic inference with phruta and RAxML #####

### Run RAxML using alignments with outliers removed
tree.raxml(folder = '2.Alignments', 
           FilePatterns = 'Masked_', 
           raxml_exec = 'raxmlHPC', # '/Users/jpintole/standard-RAxML-master/raxmlHPC-PTHREADS
           Bootstrap = 100, 
           outgroup = "Ginkgo_biloba", 
           partitioned = TRUE)

### Run RAxML using alignments with outliers 
tree.raxml(folder = '2.Alignments_withOutliers', 
           FilePatterns = 'Masked_', 
           raxml_exec = 'raxmlHPC', 
           Bootstrap = 100, 
           outgroup = "Ginkgo_biloba", 
           partitioned = TRUE)

##### ---------- Tree dating in phruta and treePL ---------- #####

##### Using tree inference with outliers removed #####

## remove duplicates from the taxonomy table
tax <- read_csv("1.CuratedSequences/1.TaxonomyOriginal.csv")

tax2 <- tax[duplicated(tax$species_names), ]

tax2 <- tax %>% 
  distinct(species_names, .keep_all = TRUE)

write_csv(tax2, file = "1.CuratedSequences/1.Taxonomy.csv")

### Tree dating using treePL
tree.dating(taxonomyFolder = "1.CuratedSequences", 
            phylogenyFolder = "3.Phylogeny", 
            scale = 'treePL')

##### Using tree inference with outliers removed #####

## remove duplicates from the taxonomy table
tax <- read_csv("1.CuratedSequences_withOutliers/1.TaxonomyOriginal.csv")

tax2 <- tax[duplicated(tax$species_names), ]

tax2 <- tax %>% 
  distinct(species_names, .keep_all = TRUE)

write_csv(tax2, file = "1.CuratedSequences_withOutliers/1.Taxonomy.csv")

### Tree dating using treePL
tree.dating(taxonomyFolder = "1.CuratedSequences_withOutliers", 
            phylogenyFolder = "3.Phylogeny_withOutliers", 
            scale = 'treePL')

##### Visualization #####
trCDR <- read.tree("4.Timetree/family-levelCalibration.tre")

write.nexus(trCDR, file = "4.Timetree/CDR_calibratedPhylo.nex")

plot(trCDR, show.tip.label = FALSE)
axisPhylo()

##### edit tree based on the taxonomy of CDR #####

trCDR <- read.tree("Dropbox/Collaborations/CDR/CDR_phylogeny_v3/4.Timetree/family-levelCalibration.tre")

taxoCDR <- readxl::read_excel("Dropbox/Collaborations/CDR/CDR_phylogeny_v3/5.CkeckTaxonomy/CDR_PHY_species.xlsx")

taxoCDR <- taxoCDR %>% 
  filter(scinamePHY != "NA")

