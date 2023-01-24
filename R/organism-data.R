#' @importFrom FELLA buildGraphFromKEGGREST buildDataFromGraph

buildGraph <- function(organism,
                       methods = availableMethods(),
                       filter_path = NULL,
                       database_directory = NULL,
                       internal_directory = TRUE,
                       damping_factor = 0.85,
                       niter = 100){
  
  organism_graph <- buildGraphFromKEGGREST(
    organism = organism, 
    filter.path = filter_path)
  
  buildDataFromGraph(
    keggdata.graph = organism_graph,
    databaseDir = database_directory,
    internalDir = internal_directory,
    matrices = methods,
    normality = c("diffusion", "pagerank"),
    dampingFactor = damping_factor,
    niter = niter)
}

#' Organism data
#' @description Load or build the organism specific KEGG datafor enrichment analysis.
#' This provides a convenience wrapper around `FELLA::loadKEGGdata`, `FELLA::buildGraphFromKEGGREST` and
#' `FELLA::buildDataFromGraph`. See the documentation of these functions for further information.
#' @param organism the KEGG code for the organism of interest
#' @param methods the enrichment techniques to build. Any returned by `availableMethods`.
#' @param filter_path argument to pass to argument `filter.path` of `FELLA::buildGraphFromKEGGREST`. A vector of regular expressions to match pathways to exclude.
#' @param database_directory argument to pass to argument `databaseDir` of `FELLA::buildDataFromGraph` and `FELLA::loadKEGGdata`. The directory name/path in which to save the KEGG data.
#' @param internal_directory logical. Argument to pass to argument `internalDir` of `FELLA::buildDataFromGraph` and `FELLA::loadKEGGdata`. Should the save directory be internal to the FELLA package directory?
#' @param damping_factor argument to pass to argument `dampingFactor` of `FELLA::buildDataFromGraph`. A value between 0 and 1 for the PageRank damping factor.
#' @param niter argument to pass to argument `niter` of `FELLA::buildDataFromGraph`. A value between 10 and 1000. The number of iterations to estimate the values for CC size.
#' @return An object of S4 class `FELLA.DATA` containing the KEGG data for the specified organism.
#' @examples 
#' ## Load the example organism data available from within the package
#' organismData(
#'   'bdi',
#'   database_directory = system.file(
#'     'bdi',
#'     package = 'riches'),
#'   internal_directory = FALSE
#' )
#' @importFrom FELLA loadKEGGdata
#' @export

organismData <- function(organism,
                          methods = availableMethods(),
                          filter_path = NULL,
                          database_directory = organism,
                          internal_directory = TRUE,
                          damping_factor = 0.85,
                          niter = 100){
  
  organism_graph <- try(
    loadKEGGdata(
      databaseDir = database_directory, 
      internalDir = internal_directory,
      loadMatrix = c('diffusion','pagerank')),
    silent = TRUE
  )
  
  if (inherits(organism_graph,'try-error')){
    message('KEGG data not found. Building...')
    
    buildGraph(organism,
               methods = methods,
               filter_path = filter_path,
               database_directory = database_directory,
               internal_directory = internal_directory,
               damping_factor = damping_factor,
               niter = niter)
    
    organism_graph <- loadKEGGdata(
      databaseDir = database_directory, 
      internalDir = internal_directory,
      loadMatrix = methods[methods != 'hypergeom'])
  }
  
  return(organism_graph)
}

setGeneric('organismCompounds',function(organism_data){
  standardGeneric('organismCompounds')
})

#' @importFrom magrittr %>%
#' @importFrom tidygraph as_tbl_graph
#' @importFrom FELLA getGraph
#' @importFrom dplyr filter
#' @importFrom igraph vertex.attributes
#' @importFrom tibble as_tibble
#' @importFrom purrr map_chr

setMethod('organismCompounds',signature = 'FELLA.DATA',
          function(organism_data){
            organism_data %>%
              getGraph() %>%
              as_tbl_graph() %>%
              activate(nodes) %>%
              filter(com == 5) %>% 
              select(name,NAME) %>% 
              vertex.attributes() %>% 
              as_tibble() %>% 
              mutate(
                NAME = map_chr(NAME,~.x[[1]])
              )
          } 
)

