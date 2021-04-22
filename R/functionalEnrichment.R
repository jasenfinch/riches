#' functionalEnrichment
#' @rdname functionalEnrichment
#' @description Functional enrichment.
#' @param analysis S4 object of class Analysis
#' @param assignment S4 object of class Assignment
#' @param parameters S4 object of class EnrichmentParameters
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
#' @importFrom mzAnnotation metaboliteDB descriptors filterACCESSIONS filterMF filterIP getAccessions
#' @importFrom MFassign assignments
#' @importFrom dplyr select distinct bind_rows
#' @importFrom purrr map
#' @importFrom metabolyseR analysisResults explanatoryFeatures
#' @importFrom magrittr set_names
#' @export

setMethod('functionalEnrichment',signature = signature(analysis = 'Analysis',assignment = 'Assignment',parameters = 'EnrichmentParameters'),
          function(analysis,assignment,parameters){
            
            FELLA <- organismNetwork(parameters@organism)
            
            adductRules <- assignment %>%
              .@parameters %>%
              .@adductRules
            
            oc <- organismCompounds(FELLA)
            
            oc <- metabolites %>%
              filterACCESSIONS(oc$name) 
            
            mfs <- assignment %>%
              assignments() %>%
              select(Name,MF,Adduct) %>%
              distinct()
            
            MFhits <- mfs %>%
              split(1:nrow(.)) %>%
              map(~{
                m <- .
                oc %>%
                  filterMF(m$MF) %>%
                  filterIP(adductRules$Rule[adductRules$Name == m$Adduct]) %>%
                  getAccessions() %>%
                  mutate(Name = m$Name,MF = m$MF,Adduct = m$Adduct)
              }) %>%
              bind_rows()
            
            bc <- MFhits$ID %>%
              unique()
            
            explanFeat <- analysis %>%
              analysisResults('modelling') %>%
              {.[[parameters@features$model]][[parameters@features$predictor]]} %>%
              explanatoryFeatures(threshold = parameters@features$threshold)
            
            comparisons <- explanFeat$Comparison %>%
              unique()
            
            enrichRes <- comparisons %>%
              map(~{
                p <- .
                message(str_c('\n',p))
                feat <- explanFeat %>%
                  filter(Comparison == p)
                ec <- MFhits %>%
                  filter(Name %in% feat$Feature) %>%
                  .$ID %>%
                  unique()
                if (length(ec) > 0) {
                  comp <- defineCompounds(
                    compounds = ec,
                    compoundsBackground = bc,
                    data = FELLA)
                  
                  if ('hypergeom' %in% functional(parameters)$methods) {
                    comp <- comp %>%
                      runHypergeom(data = FELLA)
                  }
                  
                  if ('diffusion' %in% functional(parameters)$methods) {
                    comp <- comp %>%
                      runDiffusion(data = FELLA)    
                  }
                  
                  if ('pagerank'%in% functional(parameters)$methods) {
                    comp <- comp %>%
                      runPagerank(data = FELLA)
                  }
                  
                  return(comp)
                } else {
                  message('No explanatory features assigned.')
                }
              }) %>%
              set_names(comparisons)
            
            new('FunctionalEnrichment',
                network = FELLA,
                hits = MFhits,
                explanatory = explanFeat,
                results = enrichRes
            )
          }
)