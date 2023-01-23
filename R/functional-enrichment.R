#' Available functional enrichment methods
#' @description Methods available for functional enrichment using the FELLA package
#' @return A character vector of available methods.
#' @examples 
#' availableMethods()
#' @export

availableMethods <- function(){
  c('hypergeom','diffusion','pagerank')
}

setClass('FunctionalEnrichment',
         slots = list(
           graph = 'FELLA.DATA',
           hits = 'tbl_df',
           explanatory = 'tbl_df',
           results = 'list'
         ),
         contains = 'RandomForest'
)

#' functional enrichment
#' @rdname functionalEnrichment
#' @description Functional enrichment.
#' @param x object of S4 class RandomForest
#' @param organism
#' @param methods
#' @param adduct_rules_table
#' @param organism_graph
#' @examples 
#' \dontrun{
#' ## Generate enrichment parameters
#' parameters <- enrichmentParameters('bdi')
#' 
#' ## Select only "diffusion" enrichment
#' functional(parameters) <- list(methods = 'diffusion')
#' 
#' ## Run functional enrichment
#' fe <- functionalEnrichment(example_analysis,example_assignment,parameters)
#' }
#' @importFrom FELLA defineCompounds runHypergeom runDiffusion runPagerank
#' @importFrom cheminf metaboliteDB descriptors filterEntries filterMF filterIP entries
#' @importFrom assignments assignments
#' @importFrom dplyr select distinct bind_rows
#' @importFrom purrr map
#' @importFrom metabolyseR analysisResults explanatoryFeatures
#' @importFrom magrittr set_names
#' @export

setGeneric('functionalEnrichment',function(x,
                                           organism,
                                           methods = availableMethods(),
                                           adduct_rules_table = mzAnnotation::adduct_rules(),
                                           organism_graph = organismGraph(organism))
  standardGeneric('functionalEnrichment')
)

#' @rdname functionalEnrichment

setMethod(
  'functionalEnrichment',
  signature = 'RandomForest',
  function(
    x,
    organism,
    methods = availableMethods(),
    adduct_rules_table = mzAnnotation::adduct_rules(),
    organism_graph = organismGraph(organism)
  ){
    
    organism_compounds <- organismCompounds(organism_graph)
    
    organism_database <- metabolites %>%
      filterEntries(organism_compounds$name) 
    
    feature_assignments <- x %>% 
      extractAssignments() %>% 
      filter(!is.na(MF))
    
    mf_hits <- feature_assignments %>%
      dplyr::rowwise() %>% 
      dplyr::group_map(~{
        adduct_rule <- adduct_rules_table %>% 
          filter(Name == .x$adduct) %>% 
          .$Rule %>% 
          rlang::parse_expr()
        
        organism_database %>%
          filterMF(.x$MF) %>%
          filterIP(rule = !!adduct_rule) %>%
          entries() %>%
          mutate(name = .x$name,
                 feature = .x$feature,
                 adduct = .x$adduct)
      }) %>%
      bind_rows()
    
    background_compounds <- mf_hits$ID %>%
      unique()
    
    explanatory_features <- x %>%
      explanatoryFeatures()
    
    enrichment_results <- explanatory_features %>%
      dplyr::group_by(comparison) %>% 
      dplyr::group_map(~{
        message(.x$comparison[1])
        
        explanatory_compounds <- mf_hits %>%
          filter(name %in% .x$feature) %>%
          .$ID %>%
          unique()
        
        if (length(explanatory_compounds) > 0) {
          comparison_enrichment <- defineCompounds(
            compounds = explanatory_compounds,
            compoundsBackground = background_compounds,
            data = organism_graph)
          
          if ('hypergeom' %in% methods) {
            comparison_enrichment <- comparison_enrichment %>%
              runHypergeom(data = organism_graph)
          }
          
          if ('diffusion' %in% methods) {
            comparison_enrichment <- comparison_enrichment %>%
              runDiffusion(data = organism_graph)    
          }
          
          if ('pagerank'%in% methods) {
            comparison_enrichment <- comparison_enrichment %>%
              runPagerank(data = organism_graph)
          }
          
          return(comparison_enrichment)
        } else {
          message('No explanatory features assigned.')
        }
      },.keep = TRUE) %>% 
      set_names(unique(explanatory_features$comparison))
    
    
    new('FunctionalEnrichment',
        graph = organism_graph,
        hits = mf_hits,
        explanatory = explanatory_features,
        results = enrichment_results
    )            
    
  })
