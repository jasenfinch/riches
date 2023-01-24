
#' @importFrom metabolyseR features
#' @importFrom tidyr separate
#' @importFrom dplyr across

extractAssignments <- function(x){
  feature_assignments <- tibble(
    name = features(x)
  ) %>% 
    separate(
      name,
      c('feature','MF','isotope','adduct'),
      sep = ' ',
      fill = 'right',
      remove = FALSE
    ) %>% 
    mutate(across(
      .fns = ~replace(.x,
                      .x == "",
                      NA)
    ))
}

#' @importFrom rlang parse_expr
#' @importFrom dplyr left_join rename

pips <- function(x,organism_data,adduct_rules_table){
  organism_compounds <- organismCompounds(organism_data)
  
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
        parse_expr()
      
      organism_database %>%
        filterMF(.x$MF) %>%
        filterIP(rule = !!adduct_rule) %>%
        entries() %>%
        mutate(name = .x$name,
               feature = .x$feature,
               adduct = .x$adduct)
    }) %>%
    bind_rows() %>% 
    left_join(organism_compounds %>% 
                rename(ID = name,
                       compound = NAME),
              by = 'ID')
}
