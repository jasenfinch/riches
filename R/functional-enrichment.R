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
           organism_data = 'FELLA.DATA',
           hits = 'tbl_df',
           explanatory = 'tbl_df',
           results = 'list'
         ),
         contains = c('RandomForest',
                      'FELLA.DATA')
)

#' @importFrom methods as show

setMethod('show',signature = 'FunctionalEnrichment',
          function(object){
            show(as(object,'RandomForest'))
            show(as(object,'FELLA.DATA'))
            message()
            cat(length(unique(hits(object)$feature)),'m/z features matched to KEGG compounds.\n')
            cat(length(unique(explanatoryFeatures(object)$feature)),'explanatory m/z features.')
            message()
            
            purrr::iwalk(enrichmentResults(object),
                         ~{
                           message()
                           cat(.y,'\n')
                           show(.x)
                         })
          })

#' FunctionalEnrichment S4 class accessors
#' @rdname functional-accessors
#' @description Accessor methods for the `FunctionalEnrichment` S4 class.
#' @param x object of S4 class `FunctionalEnrichment`
#' @param method the method results to access. One of `availableMethods`.
#' @param nlimit argument to pass to argument `nlimit` of `FELLA::generateResultsTable`. Limits the order of the sub-graph solutions for methods `diffusion` and `pagerank`. 
#' @param ... ignored
#' @return A tibble or a list of objects of `FELLA.USER` S4 class depending on the method used.
#' @examples
#' ## Perform random forest on the example data 
#' random_forest <- assigned_data %>% 
#' metabolyseR::randomForest(
#'   cls = 'class'
#' )
#' 
#' ## Perform functional enrichment analysis
#' enrichment_results <- functionalEnrichment(
#'   random_forest,
#'   'bdi',
#'   methods = 'hypergeom',
#'   organism_data = organismData(
#'     'bdi',
#'     database_directory = system.file(
#'       'bdi',
#'       package = 'riches'),
#'     internal_directory = FALSE
#'   )
#' )
#' 
#' ## Access the m/z feature KEGG compound matches
#' hits(enrichment_results)
#' 
#' ## Access the explanatory features used for functional enrichment
#' explanatoryFeatures(enrichment_results)
#' 
#' ## Access the FELLA.USER functional enrichment object
#' enrichmentResults(enrichment_results)
#'
#' ## Extract a table of enrichment results
#' generateResultsTable(enrichment_results)
#' @export

setGeneric('hits',function(x)
  standardGeneric('hits')
)

#' @rdname functional-accessors

setMethod('hits',signature = 'FunctionalEnrichment',
          function(x){
            x@hits
          })

#' @rdname functional-accessors
#' @export

setMethod('explanatoryFeatures',signature = 'FunctionalEnrichment',
          function(x){
            x@explanatory
          })

#' @rdname functional-accessors
#' @export

setGeneric('enrichmentResults',function(x)
  standardGeneric('enrichmentResults')
)

#' @rdname functional-accessors

setMethod('enrichmentResults',signature = 'FunctionalEnrichment',
          function(x){
            x@results
          })

#' @rdname functional-accessors
#' @export

setGeneric('generateResultsTable',function(
    x,
    method = availableMethods(),
    nlimit = 250,
    ...)
  standardGeneric('generateResultsTable')
)

#' @rdname functional-accessors
#' @importFrom FELLA generateResultsTable

setMethod('generateResultsTable',signature = 'FunctionalEnrichment',
          function(x,
                   method = availableMethods(),
                   nlimit = 250){

            method <- match.arg(
              method,
              choices = availableMethods()
            )

            x %>%
              enrichmentResults() %>%
              map(~{
                if (!is.null(.x)){
                  FELLA::generateResultsTable(
                    method = method,
                    nlimit = nlimit,
                    LabelLengthAtPlot = 50,
                    object = .x,
                    data = x
                  ) %>%
                    tibble::as_tibble() 
                }
              }) %>%
              bind_rows(.id = 'comparison')
          })

#' Functional enrichment
#' @rdname functionalEnrichment
#' @description Perform functional enrichment analyses of explanatory features using the 
#' [{FELLA}](https://bioconductor.org/packages/release/bioc/html/FELLA.html) R package.
#' @param x object of S4 class `RandomForest`
#' @param organism the KEGG code for the organism of interest
#' @param methods the enrichment techniques to build. Any returned by `availableMethods`.
#' @param split split the explanatory features into further groups based on their trends. See details.
#' @param organism_data an object of S4 class `FELLA.DATA`
#' @param adduct_rules_table the adduct ionisation rules for matching m/z features to KEGG compounds. 
#' Format should be as returned from `mzAnnotation::adduct_rules`.
#' @param ... arguments to pass to `metabolyseR::explanatoryFeatures`
#' @details 
#' For argument `split = 'trends'`, the explanatory features can be split into further groups 
#' based on their trends. This is not supported for unsupervised random forest.
#' 
#' For random forest classification, this is for binary comparisons only. Functional enrichment 
#' is performed seperately on the up and down regulated explanatory features for each comparison. The 
#' `up regulated` and `down regulated` groups are based on the trends of log2 ratios between 
#' the comparison classes. `up regulated` explanatory features have a higher median intensity 
#' in the right-hand class compared to the left-hand class of the comparison. The opposite is true
#' for the `down regulated` explanatory features.
#' 
#' For random forest regression, the explanatory features are split based on their Spearman's 
#' correlation coefficient with the response variable prior to functional enrichment analysis
#' giving `positively correlated` and `negatively correlated` subgroups.
#' @return An object of S4 class `FunctionalEnrichment`.
#' @examples 
#' ## Perform random forest on the example data 
#' random_forest <- assigned_data %>% 
#' metabolyseR::randomForest(
#'   cls = 'class'
#' )
#' 
#' ## Perform functional enrichment analysis
#' functionalEnrichment(
#'   random_forest,
#'   'bdi',
#'   methods = 'hypergeom',
#'   organism_data = organismData(
#'     'bdi',
#'     database_directory = system.file(
#'       'bdi',
#'       package = 'riches'),
#'     internal_directory = FALSE
#'   )
#' )
#' 
#' ## An example using split trends
#' ## Perform binary random forest classification on the example data 
#' random_forest <- assigned_data %>% 
#'   metabolyseR::randomForest(
#'     cls = 'class',
#'     binary = TRUE
#'   )
#' 
#' ## Perform functional enrichment analysis
#' functionalEnrichment(
#'   random_forest,
#'   'bdi',
#'   methods = 'hypergeom',
#'   split = 'trends',
#'   organism_data = organismData(
#'     'bdi',
#'     database_directory = system.file(
#'       'bdi',
#'       package = 'riches'),
#'     internal_directory = FALSE
#'   )
#' )
#' @importFrom FELLA defineCompounds runHypergeom runDiffusion runPagerank
#' @importFrom cheminf metaboliteDB descriptors filterEntries filterMF filterIP entries
#' @importFrom dplyr select distinct bind_rows group_keys
#' @importFrom purrr map
#' @importFrom metabolyseR analysisResults explanatoryFeatures type
#' @importFrom magrittr set_names
#' @importFrom mzAnnotation adduct_rules
#' @importFrom methods new
#' @export

setGeneric('functionalEnrichment',function(x,
                                           organism,
                                           methods = availableMethods(),
                                           split = c('none','trends'),
                                           organism_data = organismData(organism),
                                           adduct_rules_table = adduct_rules(),
                                           ...)
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
    split = c('none','trends'),
    organism_data = organismData(organism),
    adduct_rules_table = adduct_rules(),
    ...
  ){
    
    rf_type <- type(x)
    
    methods <- match.arg(methods,
                        choices = availableMethods(),
                        several.ok = TRUE)
    
    split <- match.arg(split,
                       choices = c(
                         'none',
                         'trends'
                       ))
    
    explanatory_features <- explanatoryFeatures(x,...)
    
    if (split == 'trends'){
      explanatory_features <- trends(x,explanatory_features)
    }
    
    mf_hits <- pips(x,
                    organism_data,
                    adduct_rules_table)
    
    background_compounds <- mf_hits$ID %>%
      unique()
    
    if (rf_type == 'classification'){
      explanatory_features <- explanatory_features %>% 
        group_by(response,comparison)
    }
    
    if (rf_type == 'regression'){
      explanatory_features <- explanatory_features %>% 
        group_by(response)
    }
    
    if (split == 'trends'){
      explanatory_features <- explanatory_features %>% 
        group_by(trend,.add = TRUE)
    }
    
    enrichment_results <- explanatory_features %>%
      dplyr::group_map(~{
        message()
        
        if (rf_type != 'unsupervised'){
          message(.x$response[1])
        }
        
        if (rf_type == 'classification'){
          message(.x$comparison[1]) 
        }
        
        if (split == 'trends'){
          message(.x$trend[1]) 
        }
        
        explanatory_compounds <- mf_hits %>%
          filter(name %in% .x$feature) %>%
          .$ID %>%
          unique()
        
        if (length(explanatory_compounds) > 0) {
          comparison_enrichment <- defineCompounds(
            compounds = explanatory_compounds,
            compoundsBackground = background_compounds,
            data = organism_data)
          
          if ('hypergeom' %in% methods) {
            comparison_enrichment <- comparison_enrichment %>%
              runHypergeom(data = organism_data)
          }
          
          if ('diffusion' %in% methods) {
            comparison_enrichment <- comparison_enrichment %>%
              runDiffusion(data = organism_data)    
          }
          
          if ('pagerank'%in% methods) {
            comparison_enrichment <- comparison_enrichment %>%
              runPagerank(data = organism_data)
          }
          
          return(comparison_enrichment)
        } else {
          message('No assigned explanatory m/z features matched to KEGG compounds.')
        }
      },.keep = TRUE)
    
    
    if (rf_type != 'unsupervised') {
      group_names <- group_keys(explanatory_features) %>% 
        tidyr::unite(name) %>% 
        .$name
      
      enrichment_results <- enrichment_results %>% 
        set_names(group_names)
    }
    
    
    results <- new(
      'FunctionalEnrichment',
      x,
      organism_data,
      hits = mf_hits,
      explanatory = explanatory_features,
      results = enrichment_results
    )            
    
    return(results)
  })
