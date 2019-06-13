
setMethod('availablePairwises',signature = 'FunctionalEnrichment',
          function(x){
            x@explanatory %>%
              .$Pairwise %>%
              unique()
          }
)