
setGeneric('organismCompounds',function(FELLA){
  standardGeneric('organismCompounds')
})

#' @rdname functionalEnrichment
setGeneric('functionalEnrichment',function(analysis,assignment,parameters){
  standardGeneric('functionalEnrichment')
})


setGeneric('plotGraph',function(x,comparison,type = 'diffusion'){
  standardGeneric('plotGraph')
})

setGeneric('availablePairwises',function(x){
  standardGeneric('availablePairwises')
})
