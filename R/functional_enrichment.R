
extractAssignments2 <- function(features) {
  tibble::tibble(
    name = features
  ) %>% 
    tidyr::separate(
      name,
      c('feature', 'MF', 'isotope', 'adduct'),
      sep = ' ',
      fill = 'right',
      remove = FALSE
    ) %>% 
    dplyr::mutate(
           dplyr::across(
                  .cols = dplyr::everything(),
                  .fns = ~replace(
                                  .x,
                                  .x == "",
                                  NA
                  )
    ))
}

pips2 <- function(features, organism_data, adduct_rules_table) {
  organism_compounds <- organismCompounds(organism_data)
  
  organism_database <- metabolites %>%
    cheminf::filterEntries(organism_compounds$name) 
  
  feature_assignments <- features %>% 
    extractAssignments2() %>% 
    dplyr::filter(!is.na(MF))
  
  feature_assignments %>%
    dplyr::rowwise() %>% 
    dplyr::group_map(~{
      adduct_rule <- adduct_rules_table %>% 
        dplyr::filter(Name == .x$adduct) %>% 
        .$Rule %>% 
        rlang::parse_expr()
      
      organism_database %>%
        cheminf::filterMF(.x$MF) %>%
        cheminf::filterIP(rule = !!adduct_rule) %>%
        cheminf::entries() %>%
        dplyr::mutate(name = .x$name,
               feature = .x$feature,
               adduct = .x$adduct)
    }) %>%
    dplyr::bind_rows() %>% 
    dplyr::left_join(organism_compounds %>% 
                dplyr::rename(ID = name,
                       compound = NAME),
              by = 'ID')
}

functional_enrichment <- function(
          features,
          explanatory_features,
          organism,
          methods = availableMethods(),
          organism_data = organismData(organism),
          adduct_rules_table = adduct_rules(),
          ...
) {

    methods <- match.arg(methods,
                        choices = availableMethods(),
                        several.ok = TRUE)
    
    mf_hits <- pips2(features,
                    organism_data,
                    adduct_rules_table)
    
    background_compounds <- mf_hits$ID %>%
      unique()
    
    explanatory_compounds <- mf_hits %>%
      filter(name %in% explanatory_features) %>%
      .$ID %>%
      unique()
    
    if (length(explanatory_compounds) == 0) {
        message('No assigned explanatory m/z features matched to KEGG compounds.')
        return(invisible())
    }

    enrichment_results <- defineCompounds(
      compounds = explanatory_compounds,
      compoundsBackground = background_compounds,
      data = organism_data)
    
    if ('hypergeom' %in% methods) {
      enrichment_results <- enrichment_results %>%
        runHypergeom(data = organism_data)
    }
    
    if ('diffusion' %in% methods) {
      enrichment_results <- enrichment_results %>%
        runDiffusion(data = organism_data)    
    }
    
    if ('pagerank' %in% methods) {
      enrichment_results <- enrichment_results %>%
        runPagerank(data = organism_data)
    }
    
    return(enrichment_results)
}
