#' enrichmentParameters
#' @description Enrichment parameters.
#' @param organism KEGG organism ID
#' @param model model type to use for extracting explanatory features
#' @param response response variable name to use for extracting explanatory features
#' @param threshold explanatory feature cutoff
#' @importFrom methods new
#' @export

enrichmentParameters <- function(organism, model = 'randomForest', response = 'class', threshold = 0.05){
  new('EnrichmentParameters',
      organism = organism,
      features = list(
        model = model,
        response = response,
        threshold = threshold
      ),
      structural = list()
  )
}