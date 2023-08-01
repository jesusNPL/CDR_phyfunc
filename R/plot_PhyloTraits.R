
library(tidyverse)
library(scales)

## Read imputed traits
bhpmf_mean <- read.csv("CDR_phyfunc/traits/CDR_mean_gap_filled_BHPMF.csv")

bhpmf_mean <- bhpmf_mean %>% 
  mutate(scinames = gsub(" ", "_", species)) %>% 
  select(plant_id, scinames, everything())

## Read tacted phylogeny
trt <- ape::read.tree("CDR_phyfunc/phylogeny/4.Timetree/Tacted/OrderLevel/CDR.tacted_v3.newick.tre")
trt <- ape::drop.tip(trt, setdiff(trt$tip.label, bhpmf_mean$scinames))

## Order traits according the tips of the phylogeny
traits_order <- bhpmf_mean[match(trt$tip.label, bhpmf_mean$scinames), ]
rownames(traits_order) <- traits_order$scinames

## Visualization
col2 <- alpha(c("blue", "purple", "darkgreen" , "orange", "red", 
                #"yellow", 
                "black", "brown", "darkgoldenrod", "yellow"), 0.7)
names(col2) <- c("Leaf area", "N mass", "LMA", "Plant height", 
                 "Diaspore mass", "SSD", "LDMC", "SSD_imp", "SSD_comb")

pdf(file = "CDR_phyfunc/CDR_phylo_traits.pdf", 
    height = 18, width = 20)
par(mar = c(1, 1, 1, 1))
par(oma = c(3, 3, 3, 3))

plot(ladderize(trt), TRUE, cex = 2, no.margin = TRUE,  
     edge.width = 2, type = "fan", show.tip.label = FALSE, 
     x.lim = c(-500, 500), y.lim = c(-500, 500))

legend("bottomleft", legend = names(col2)[1:7], 
       pch = 21, pt.bg = col2, bty = "n", title = "CDR-traits",
       text.col = "black", cex = 2, pt.cex = c(4.5, 4.5, 4.5, 4.5, 4.5))

tiplabels(pch = 21, bg = col2[1], col = col2[1], cex = abs(na.omit(traits_order$Leaf_Area)), 
          adj = 0.75, offset = 10)
tiplabels(pch = 21, bg = col2[2], col = col2[2], cex = abs(na.omit(traits_order$Nmass)), 
          adj = 0.75, offset = 20)
tiplabels(pch = 21, bg = col2[3], col = col2[3], cex = abs(na.omit(traits_order$LMA)), 
          adj = 0.75, offset = 30)
tiplabels(pch = 21, bg = col2[4], col = col2[4], cex = abs(na.omit(traits_order$Plant_height)), 
          adj = 0.75, offset = 40)
tiplabels(pch = 21, bg = col2[5], col = col2[5], cex = abs(na.omit(traits_order$Diaspore_mass)), 
          adj = 0.75, offset = 50)
tiplabels(pch = 21, bg = col2[6], col = col2[6], cex = abs(na.omit(traits_order$SSD_observed)), 
          adj = 0.75, offset = 60)
tiplabels(pch = 21, bg = col2[7], col = col2[7], cex = abs(na.omit(traits_order$LDMC)), 
          adj = 0.75, offset = 70)
#tiplabels(pch = 21, bg = col2[8], col = col2[8], cex = abs(na.omit(traits_order$SSD_imputed)), 
 #         adj = 0.75, offset = 80)
#tiplabels(pch = 21, bg = col2[9], col = col2[9], cex = abs(na.omit(traits_order$SSD_combined)), 
 #         adj = 0.75, offset = 90)
dev.off()
