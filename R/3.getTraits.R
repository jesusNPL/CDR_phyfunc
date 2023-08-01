library(tidyverse)

### Read Trait Data - Díaz et al. 2022
trt_dt <- readxl::read_excel("../../Vasculars/Traits_Diaz_EtAl_2022/Dataset/Species_mean_traits.xlsx", 
                             sheet = "Species_mean_traits")

trt_dt <- trt_dt %>% 
  mutate(scinameTPL = gsub(" ", "_", `Species name standardized against TPL`)) %>% 
  select(`TRY 30 AccSpecies ID`, scinameTPL, everything())

### Read species names 

cdrTax <- read_csv("CDR_phylogeny_v3/4.Timetree/Tacted/OrderLevel/CDR_taxonomy_v3.csv")

## search taxonomy using TPL 
sppTaxCRD_tpl <- Taxonstand::TPL(cdrTax$species, infra = TRUE, corr = TRUE) # it will return a lot of warnings, please do not pay attention to that. 

sppTaxCRD_tpl <- sppTaxCRD_tpl %>% 
  mutate(scinameTPL = paste0(New.Genus, "_", New.Species)) %>% 
  select(Taxon, scinameTPL)

cdrTax_tpl <- full_join(cdrTax, sppTaxCRD_tpl, 
                        by = c("species" = "Taxon"))

### Filter traits in Díaz using CDR taxonomy

trt_CDR <- trt_dt %>% 
  filter(scinameTPL %in% cdrTax_tpl$scinameTPL)

### Join taxonomic data and trait data 
trt_CDR <- full_join(cdrTax_tpl, trt_CDR, 
                     by = "scinameTPL")

write_csv(trt_CDR, 
          file = "CDR_traits/CDR_traits_SDaiz.csv")

##### Prepare data for imputation #####

## Phylogeny 
phy_CDR <- ape::read.tree("CDR_phylogeny_v3/4.Timetree/Tacted/OrderLevel/CDR.tacted_v3.newick.tre")

## Traits from Diaz
trt_CDR_4imputation <- readxl::read_excel("CDR_traits/CDR_traits_SDaiz.xlsx", 
                                          sheet = "CDR_select") 

trt_CDR_4imputation <- trt_CDR_4imputation %>% 
  mutate(Taxa = gsub(" ", "_", species)) %>% 
  select(order, family, genus, species, Taxa, everything())

##### Match trait - phylogeny #####
traits <- trt_CDR_4imputation[trt_CDR_4imputation$Taxa %in% phy_CDR$tip.label, ]
rownames(traits) <- traits$Taxa

phylo <- ape::drop.tip(phy_CDR, setdiff(phy_CDR$tip.label, traits$Taxa))

traits <- traits[match(phylo$tip.label, rownames(traits)), ]

## Get phylogenetic distance matrix derived from tree into a set of orthogonal vectors
pvr <- PVR::PVRdecomp(phylo)

# Extract the PVRs
pvrs <- pvr@Eigen$vectors

# Combine traits and PVRs
traits_pvrs <- cbind(traits, pvrs)

##### Save data for imputation #####
write_csv(traits_pvrs, file = "CDR_traits/CDR_traits_for_imputation.csv")

##### Trait imputation using BHPMF #####
library(BHPMF)

traits_pvrs <- read_csv("CDR_traits/CDR_traits_for_imputation.csv")

traits_pvrs <- traits_pvrs %>% 
  select(order, family, genus, species, # taxonomy
         Leaf_Area, Nmass, LMA, Plant_height, Diaspore_mass, SSD_observed, # traits
         LDMC, SSD_imputed, SSD_combined, # traits
         c1, c2, c3, c4, c5, c6, c7, c8, c9, c10)

### Z transformation
traits.info <- traits_pvrs %>% 
  select(Leaf_Area, Nmass, LMA, Plant_height, Diaspore_mass, SSD_observed, # traits
         LDMC, SSD_imputed, SSD_combined)

traits.info <- apply(X = traits.info, MARGIN = 2, FUN = log10)
traits.info <- apply(X = traits.info, MARGIN = 2, FUN = scale)

### Traits
traits.info <- cbind(traits.info, traits_pvrs[, 14:23])

traits.info <- as.matrix(traits.info)

### Taxonomic hierarchy
hierarchy.info <- traits_pvrs %>% 
  mutate(plant_id = 1:nrow(traits_pvrs)) %>% 
  select(plant_id, species, genus, family, order) %>% 
  as.data.frame()
  
##### RUN BHPMF #####

GapFilling(X = traits.info, 
           hierarchy.info = hierarchy.info, 
           num.samples = 10000, burn = 2000, 
           tuning = TRUE, verbose = TRUE, 
           tmp.dir = "CDR_traits/BHPMF_imputation/", 
           mean.gap.filled.output.path = "CDR_traits/BHPMF_imputation/CDR_mean_gap_filled_BHPMF.csv",
           std.gap.filled.output.path = "CDR_traits/BHPMF_imputation/CDR_std_gap_filled_BHPMF.csv")

#RRMSE for the test data:  1.205673$min.rmse
#[1] 1.370784

#$best.number.latent.features
#[1] 20


##########################################
# Usage 4: Calculate cross validation RMSE
##########################################
# Calculate average RMSE with the default values 
traitNames <- colnames(traits.info)

for(i in 1:length(traitNames)) { 
  
  print(paste0("Running CV for ", traitNames[i]), " please, be patient!")
  
  dir.create(paste0("CDR_traits/BHPMF_imputation/CV/CV_", traitNames[i]))
  
  out <- CalculateCvRmse(as.matrix(traits.info[, i]), 
                          hierarchy.info, 
                          #tuning = TRUE, 
                          num.latent.feats = 20, 
                          tmp.dir = paste0("CDR_traits/BHPMF_imputation/CV/CV_", traitNames[i]), 
                          num.samples = 1000, 
                          burn = 200, 
                          verbose = TRUE)
  
  #avg.rmse <- out1$avg.rmse
  #std.rmse <- out1$std.rmse
  
  save(out, file = paste0("CDR_traits/BHPMF_imputation/CV/", traitNames[i], "_cv.RData"))

}

## Global CV
out3 <- CalculateCvRmse(traits.info, 
                        hierarchy.info, 
                        #tuning = TRUE, 
                        num.latent.feats = 20,
                        tmp.dir = "CDR_traits/BHPMF_imputation/CV/", 
                        num.samples = 10000, 
                        burn = 2000, 
                        verbose = TRUE)

out3

save(out3, file = "CDR_traits/BHPMF_imputation/CV/CV_RMSE_global.RData")
