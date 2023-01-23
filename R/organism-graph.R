#' Organism graph
#' @description Gather organism specific KEGG graph ready for enrichment analysis
#' @param organism KEGG organism code
#' @importFrom FELLA buildGraphFromKEGGREST buildDataFromGraph loadKEGGdata listInternalDatabases
#' @importFrom stringr str_c
#' @export

organismGraph <- function(organism = 'hsa',filter.path = NULL,internalDir = TRUE){
  
  if (!(organism %in% listInternalDatabases())) {
    message('Internal database not found for this organism. Building...')
    
    g <- buildGraphFromKEGGREST(organism = organism, filter.path = filter.path)
    
    buildDataFromGraph(
      keggdata.graph = g,
      databaseDir = organism,
      internalDir = TRUE,
      matrices = c("hypergeom", "diffusion", "pagerank"),
      normality = c("diffusion", "pagerank"),
      dampingFactor = 0.85,
      niter = 1e3)  
  }
  
  FELLA.DATA <- loadKEGGdata(
    organism, 
    internalDir = TRUE,
    loadMatrix = c('diffusion','pagerank'))
  
  return(FELLA.DATA)
}

#' @importFrom magrittr %>%
#' @importFrom tidygraph as_tbl_graph
#' @importFrom FELLA getGraph
#' @importFrom assignments nodes
#' @importFrom dplyr filter

setMethod('organismCompounds',signature = 'FELLA.DATA',
          function(FELLA){
            FELLA %>%
              getGraph() %>%
              as_tbl_graph() %>%
              nodes() %>%
              filter(com == 5)
          } 
)

