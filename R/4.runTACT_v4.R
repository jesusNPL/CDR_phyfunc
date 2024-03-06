# open /Applications/RStudio.app

##### Step 1 - Run TACT on BEST RAxML tree #####

### Create a taxonomic tree
system("tact_build_taxonomic_tree CDR_taxonomic_backbone.csv --output CDR_taxonomy.tre")

### Add missing species or 
system("tact_add_taxa --backbone CDR_timeTree_secondaryCalibrations.tre --taxonomy CDR_taxonomy.tre --output CDR_timeTree.tacted")

### Check results 
system("tact_check_results CDR_timeTree_tacted.newick.tre --backbone CDR_timeTree_secondaryCalibrations.tre --taxonomy CDR_taxonomy.tre > checkresults.csv")

##### Load time (secondary) calibrated trees #####

### Sample of 100 trees out of 1000 from RAxML bootstraps
tr_boots <- ape::read.tree("CDR_calibrated_bootstraps.tre")

### Write individual bootstrap tree into working directory
for(i in 1:length(tr_boots)) { 

  # Isolate one bootstrap tree
  tr <- tr_boots[[i]]
  # write one bootstrap tree
  ape::write.tree(tr, file = paste0("CDR_calibrated_BS_", i, ".tre"))
  
}

##### Step 2 - Run TACT for bootstrap trees #####

### Inits
# Number boostraps
nBoots <- 100

### Run TACT for bootstrap trees
for(j in 1:nBoots) { 

  print(paste0("Running tact for bootstrap ", j, " ..."))
  
  system(paste0("tact_add_taxa --backbone CDR_calibrated_BS_", j, 
                ".tre --taxonomy CDR_taxonomy.tre --output CDR_timeTree_BS_", j, 
                ".tacted"))
  
  print(paste0("Check tact results for imputation ", j, " ..."))
  
  system(paste0("tact_check_results CDR_timeTree_BS_", j, 
                ".tacted.newick.tre --backbone CDR_calibrated_BS_", j, 
                ".tre --taxonomy CDR_taxonomy.tre > checkresults_BS_", j, ".csv"))
  
}

##### Step 3 - Read TACTED bootstraps and save them as multiPhylo #####
library(ape)

### Inits
# Number bootstraps
nBoots <- 100 
# Newick Infile
nwFiles <- list.files(path = ".", pattern = ".tacted.newick.tre")
# Nexus Infile
nxFiles <- list.files(path = ".", pattern = ".tacted.nexus.tre")
# Store newick tacted bootstraps
boots_newick_lst <- list()
# Store nexustacted bootstraps
boots_nexus_lst <- list()

### Run loop
for(i in 1:nBoots) { 
  
  print(paste0("Read tacted bootstraps ", i, " ..."))
  
  ### Read tacted trees
  # Newick
  nwBS <- read.tree(nwFiles[i])
  # Nexus
  nxBS <- read.nexus(nxFiles[i])
  
  ### Store trees  
  # Newick
  boots_newick_lst[[i]] <- nwBS
  # Nexus
  boots_nexus_lst[[i]] <- nxBS
  
}

##### Save multiPhylo objects #####

### Newick
class(boots_newick_lst) <- "multiPhylo" 
# Write newick
write.tree(boots_newick_lst, 
           file = "CDR_calibrated_bootstraps_tacted.tre")

### Nexus
class(boots_nexus_lst) <- "multiPhylo" 
# Write nexus
write.nexus(boots_nexus_lst, 
            file = "CDR_calibrated_bootstraps_tacted.nex")

