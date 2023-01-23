
setGeneric('organismCompounds',function(FELLA){
  standardGeneric('organismCompounds')
})

setGeneric('plotGraph',function(x,comparison,type = 'diffusion'){
  standardGeneric('plotGraph')
})

setGeneric('availablePairwises',function(x){
  standardGeneric('availablePairwises')
})

#' @rdname EnrichmentParameters-accessors

setGeneric('functional',function(x){
  standardGeneric('functional')
})

#' @rdname EnrichmentParameters-accessors

setGeneric('functional<-',function(x,value){
  standardGeneric('functional<-')
})