#' enrichmentParameters
#' @description Enrichment parameters.
#' @param organism KEGG organism ID
#' @param model model type to use for extracting explanatory features
#' @param predictor predictor predictor name to use for extracting explanatory features
#' @param threshold explanatory feature cutoff
#' @importFrom methods new
#' @export

enrichmentParameters <- function(organism, model = 'randomForest', predictor = 'class', threshold = 0.05){
  new('EnrichmentParameters',
      organism = organism,
      features = list(
        model = model,
        predictor = predictor,
        threshold = threshold
      ),
      structural = list()
  )
}