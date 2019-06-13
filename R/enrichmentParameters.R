
#' @export

enrichmentParameters <- function(organism){
  new('EnrichmentParameters',
      organism = organism,
      features = list(
        method = 'fs.rf',
        threshold = '0.01'
      ),
      functional = list(),
      structural = list()
  )
}