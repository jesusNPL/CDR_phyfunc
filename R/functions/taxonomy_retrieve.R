#library(taxize)
#library(rgbif)
#library(pbapply)
##### Function from package {phruta} #####
taxonomy.retrieve <-
  function(species_names = NULL,
           database = "gbif",
           kingdom = NULL,
           ranks =
             c("kingdom",
               "phylum",
               "class",
               "order",
               "family",
               "genus",
               "species")) {
    
    
    if (database != "gbif") {
      taxo <- taxize::classification(species_names, db = database, rows = 1)
      invisible(Taxonomy_species <-
                  as.data.frame(do.call(
                    rbind, lapply(seq_along(taxo), function(x) {
                      if (!all(is.na(taxo[[x]]))) {
                        t(data.frame(taxo[[x]][taxo[[x]][, 2] %in% ranks, 1]))
                      } else {
                        sma <- matrix(nrow = 1, ncol = length(ranks) - 1)
                        cbind(sma, species_names[x])
                      }
                    })
                  )))
      row.names(Taxonomy_species) <- NULL
      colnames(Taxonomy_species) <- ranks
      Taxonomy_species$species <-
        sub(" ", "_", Taxonomy_species$species)
      Taxonomy_species
    } else {
      ## If animals and plants
      gbifkey <-
        lapply(species_names, function(x)
          rgbif::name_backbone(name = x, kingdom = kingdom))
      keys <- pbapply::pblapply(seq_along(gbifkey), function(x) {
        if (
          as.character(gbifkey[[x]][which(names(gbifkey[[x]])
                                          == "matchType")]) == "NONE") {
          0
        } else {
          if (length(which(names(gbifkey[[x]]) == "acceptedUsageKey")) == 0) {
            as.character(gbifkey[[x]][which(names(gbifkey[[x]]) == "usageKey")])
          } else {
            as.character(gbifkey[[x]][which(names(gbifkey[[x]])
                                            == "acceptedUsageKey")])
          }
        }
      })
      
      gbif_taxonomy <-
        lapply(unlist(keys), function(x)
          as.data.frame(rgbif::name_usage(key = x)$data))
      Taxonomy_species <-
        lapply(seq_along(gbif_taxonomy), function(y) {
          sub1 <- gbif_taxonomy[[y]]
          cate <-
            t(data.frame(unlist(lapply(seq_along(ranks), function(x) {
              nu <- which(colnames(sub1) == ranks[x])
              if (length(nu) != 1) {
                NA
              } else {
                sub1[, nu]
              }
            }))))
          colnames(cate) <- ranks
          row.names(cate) <- NULL
          cate
        })
      Taxonomy_species <- do.call(rbind.data.frame, Taxonomy_species)
      Taxonomy_species$Used_names <- species_names # argument added
    }
    return(Taxonomy_species) # argument added
  }


#taxonomy.retrieve(species_names = c("Furnarius ruffus", "Tyrannus melancholicus"))
