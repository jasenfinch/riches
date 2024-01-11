structural_enrichment <- function(structural_classifications, explanatory_features, p_adjust_method = 'bonferroni') {

    structural_classifications <- structural_classifications %>%
        dplyr::mutate(
            Name = replace(
                Name,
                is.na(Name),
                Feature[is.na(Name)]
            )
        )

    if (any(!explanatory_features %in% structural_classifications$Name)) {
        explanatory_features <- explanatory_features[explanatory_features %in% structural_classifications$Name]
        warning(
            "Some explanatory features not found in the structural classifications, these will be drpped.",
            call. = FALSE)
    }

    feature_classifications <- structural_classifications %>% 
      dplyr::select(feature = Name, kingdom:dplyr::last_col(offset = 2))

    background_classification_totals <- feature_classifications %>% 
      tidyr::gather(level,classification,-feature) %>% 
      tidyr::drop_na() %>% 
      dplyr::group_by(level, classification) %>% 
      dplyr::count()

    explanatory_classifications <- feature_classifications %>%
        dplyr::filter(feature %in% explanatory_features) 

    explanatory_classification_totals <- explanatory_classifications %>%
        tidyr::gather(level, classification, -feature) %>% 
        tidyr::drop_na() %>% 
        dplyr::group_by(level, classification) %>% 
        dplyr::count()

    contingency <- explanatory_classification_totals %>% 
        dplyr::rename(`Explanatory & in class` = n) %>% 
        dplyr::left_join(background_classification_totals,
                  by = c('level', 'classification')) %>% 
        dplyr::rename(`Not explanatory & in class` = n) %>% 
        dplyr::mutate(`Not explanatory & in class` = `Not explanatory & in class` - 
                 `Explanatory & in class`) %>%
        dplyr::mutate(n = length(explanatory_features)) %>%
        dplyr::rename(`Explanatory & not in class` = n) %>% 
        dplyr::mutate(`Explanatory & not in class` = `Explanatory & not in class` - 
            `Explanatory & in class`,
            `Not explanatory & not in class` = length(unique(structural_classifications$Name)) - 
            (`Explanatory & in class` + `Not explanatory & in class`) - 
            `Explanatory & not in class`)

    ora_results <- contingency %>% 
        dplyr::rowwise() %>% 
        dplyr::group_map(
            ~dplyr::bind_cols(.x,
                overrepresentation(
                    .x$`Explanatory & in class`,
                    .x$`Not explanatory & in class`,
                    .x$`Explanatory & not in class`,
                    .x$`Not explanatory & not in class`)
            )
        ) %>% 
        dplyr::bind_rows() %>%
        dplyr::mutate(
                adjusted.p.value = p.adjust(p.value, method = p_adjust_method)
              ) %>% 
        dplyr::relocate(adjusted.p.value, .after = p.value) %>%
        dplyr::ungroup() %>% 
        dplyr::arrange(p.value)

    return(ora_results)
}
