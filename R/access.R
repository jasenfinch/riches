
setMethod('availablePairwises',signature = 'FunctionalEnrichment',
          function(x){
            x@explanatory %>%
              .$Pairwise %>%
              unique()
          }
)

#' Get and set enrichment parameters
#' @rdname EnrichmentParameters-accessors
#' @description Retrieve or set enrichment parameters.
#' @param x S4 object of class EnrichmentParameters
#' @param value value to set
#' @export

setMethod('functional',signature = 'EnrichmentParameters',
          function(x){
            x@functional
          })

#' @rdname EnrichmentParameters-accessors
#' @export

setMethod('functional<-',signature = 'EnrichmentParameters',
          function(x,value){
            x@functional <- value
            return(x)
          })