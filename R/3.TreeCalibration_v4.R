library(geiger)

##### Tree calibration using reference tree #####

### Read taxonomy
taxonomy <- read.csv("1.CuratedSequences/1.Taxonomy.csv")
row.names(taxonomy) <- taxonomy$species_names

### read target tree - tree estimated using RAxML
# Best-scoring ML tree with support values
TargetTree <- read.tree("3.1.Raxml1000/RAxML_bipartitions.phruta.tre") 

### Read reference tree - Jin and Qian megaphylogeny
megaJQ <- read.tree("../MegaPhylos/Jin_&_Qian_2019_Ecogr/GBOTB.extended.tree")

### Run tree calibration using GEIGER - automagically generating secondary calibrations
resphy <- congruify.phylo(reference = megaJQ, # reference tree
                          target = TargetTree, # target tree, AKA our tree
                          taxonomy = base::as.matrix(taxonomy), 
                          scale = "treePL", # the target tree will be smoothed by treePL
                          ncores = 10 # number of cores 
                          )

##### Save results #####

### Save calibration points

write.csv(resphy$calibrations, 
          file = "4.TimeTree/CDR_calibrationsPoints.csv")

### Save calibrated tree using secondary calibrations

# Nexus format
write.nexus(ladderize(resphy$phy), 
            file = "4.TimeTree/CDR_timeTree_secondaryCalibrations.nex")
# Newick format
write.tree(ladderize(resphy$phy), 
           file = "4.TimeTree/CDR_timeTree_secondaryCalibrations.tre")

##### Calibrate Bootstrap samples #####

### Auxiliary function 
bootstrapCalibration <- function(Taxonomy, # Taxonomy 
                                 bootstrapTrees, # bootstrap trees 
                                 referenceTree, # tree to extract calibrations
                                 nBoots = 25L, # number of trees 
                                 directory = "4.TimeTree/bootstraps/", 
                                 nCores = 20L) { # number of cores

  ### Inits 
  treeList <- list()
  calibrationList <- list()
  
  sampTrees <- sample(bootstrapTrees, nBoots)
  
  for(i in 1:nBoots) { 
    
    print(paste0("Running secondary calibrations for bootstrap N = ", i, " ..."))
    
    tree <- sampTrees[[i]] 
    
    resphy <- congruify.phylo(reference = referenceTree, # reference tree
                              target = tree, # target tree, AKA our tree
                              taxonomy = base::as.matrix(taxonomy), 
                              scale = "treePL", # the target tree will be smoothed by treePL
                              ncores = nCores # number of cores 
    ) 
    
    ### save secondary calibrations 
    cal <- resphy$calibrations 
    cal$bootstrap <- i
    calibrationList[[i]] <- cal
    
    ### save calibrated bootstrap trees
    res <- ladderize(resphy$phy) 
    treeList[[i]] <- res

  }
  
  timeTreeBootstrap <- treeList
  class(timeTreeBootstrap) <- "multiPhylo" 
  
  ### Save calibrated bootstraps 
  # Nexus format
  write.nexus(timeTreeBootstrap, 
              file = paste0(directory, "CDR_calibrated_bootstraps.nex"))
 
  # Newick format
  write.tree(timeTreeBootstrap, 
              file = paste0(directory, "CDR_calibrated_bootstraps.tre"))
  
  ### Save secondary calibrations
  calibrations <- do.call(rbind, calibrationList)
  
  write.csv(calibrations, 
            file = paste0(directory, "CDR_calibrated_bootstraps.csv"), 
            row.names = FALSE)

}

### Read bootstraps
cdr_boots <- read.tree("3.Phylogeny100/RAxML_bootstrap.phruta")

### Reference tree 
cdr_reference <- drop.tip(megaJQ, setdiff(megaJQ$tip.label, cdr_boots[[1]]$tip.label))

### Run calibration

bootstrapCalibration(Taxonomy = taxonomy, 
                     bootstrapTrees = cdr_boots, 
                     referenceTree = cdr_reference, 
                     nBoots = 100L, # number of trees 
                     directory = "4.TimeTree/bootstraps/", 
                     nCores = 20L)

