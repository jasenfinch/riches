
#' @importFrom broom tidy

overrepresentation <- function(explanatory_in_class,
                               not_explanatory_in_class,
                               explanatory_not_class,
                               not_explanatory_not_class){
  contingency = matrix(
    c(explanatory_in_class,
      not_explanatory_in_class,
      explanatory_not_class,
      not_explanatory_not_class),
    nrow = 2,
    ncol = 2)
  
  fisher.test(
    contingency,
    alternative = 'greater'
  ) %>% 
    tidy()
}

featureInfo <- function(x){
  feature_info <- tibble(
    Name = features(x)
  ) %>% 
    separate(
      Name,
      c('Feature',
        'MF',
        'Isotope',
        'Adduct'),
      sep = ' ',
      fill = 'right',
      remove = FALSE
    )
  
  feature_info[feature_info == ''] <- NA
  
  return(feature_info)
}

setClass('StructuralEnrichment',
         slots = list(
           explanatory = 'tbl_df',
           classifications = 'tbl_df',
           results = 'tbl_df'
         ),
         contains = c('RandomForest'))

setMethod('show',signature = 'StructuralEnrichment',
          function(object){
            show(as(object,'RandomForest'))
            
            cat(length(unique(explanatoryFeatures(object)$feature)),'explanatory m/z features.\n')
            structuralClassifications(object) %>% 
              select(-feature) %>% 
              gather(level,classification) %>%
              drop_na() %>% 
              distinct() %>% 
              nrow() %>% 
              {cat(.,'structural classes total.\n')}
            enrichmentResults(object) %>% 
              filter(adjusted.p.value < 0.05) %>% 
              nrow() %>% 
              cat(.,'significantly enriched structural classes.')
          })

#' StructuralEnrichment S4 class accessors
#' @rdname structural-accessors
#' @description Accessor methods for the `StructuralEnrichment` S4 class.
#' @param x object of S4 class `StructuralEnrichment`
#' @return A tibble containing the accessed information.
#' @examples
#' ## Perform random forest on the example data 
#' random_forest <- assigned_data %>% 
#'   metabolyseR::randomForest(
#'     cls = 'class'
#'   )
#' 
#' ## Perform functional enrichment analysis
#' enrichment_results <- structuralEnrichment(
#'   random_forest,
#'   structural_classifications
#' )
#' 
#' ## The m/z feature structural classifcations
#' structuralClassifications(enrichment_results)
#' 
#' ## Access the explanatory features used for structural enrichment analysis
#' explanatoryFeatures(enrichment_results)
#' 
#' ## Access the structral enrichment analysis results
#' enrichmentResults(enrichment_results)
#' @export

setGeneric('structuralClassifications',function(x)
  standardGeneric('structuralClassifications'))

#' @rdname structural-accessors

setMethod('structuralClassifications',signature = 'StructuralEnrichment',
          function(x){
            x@classifications
          })

#' @rdname structural-accessors
#' @export

setMethod('explanatoryFeatures',signature = 'StructuralEnrichment',
          function(x){
            x@explanatory
          })

#' @rdname structural-accessors

setMethod('enrichmentResults',signature = 'StructuralEnrichment',
          function(x){
            x@results
          })

#' Structural enrichment
#' @rdname structuralEnrichment
#' @description Perform structural enrichment using over-representation analysis of explanatory *m/z* features from random forest.
#' @param x an object of S4 class `RandomForest`
#' @param structural_classifications the structral classifications corresponding to the *m/z* features present in the object specified for argument `x`. This should either be a tibble as returned by `construction::classifications()` or an object of S4 class `Construction`.
#' @param p_adjust_method the p-value adjustment method. One of those returned from `p.adjust.methods`.
#' @param ... arguments to pass to `metabolyseR::explanatoryFeatures()` 
#' @details 
#' Over-representation analysis is performed on the explanatory *m/z* features for each structural class within each experimental class comparison using the Fisher's Exact Test.
#' @return An object of S4 class `StructuralEnrichment`.
#' @examples 
#' ## Perform random forest on the example data 
#' random_forest <- assigned_data %>% 
#'   metabolyseR::randomForest(
#'     cls = 'class'
#'   )
#' 
#' ## Perform structural enrichment analysis using the example structural classifications
#' structuralEnrichment(
#'   random_forest,
#'   structural_classifications
#' )
#' @export

setGeneric('structuralEnrichment',function(x,
                                           structural_classifications,
                                           p_adjust_method = 'bonferroni',
                                           ...)
  standardGeneric('structuralEnrichment'))

#' @rdname structuralEnrichment
#' @importFrom dplyr inner_join group_by count group_map rowwise bind_cols relocate arrange mutate
#' @importFrom tidyr gather drop_na
#' @importFrom tidyselect last_col
#' @importFrom metabolyseR nFeatures

setMethod('structuralEnrichment',signature = c('RandomForest','tbl_df'),
          function(x,
                   structural_classifications,
                   p_adjust_method = 'bonferroni',
                   ...){
            
            rf_type <- type(x)
            
            feature_classifications <- structural_classifications %>% 
              select(-Name) %>% 
              inner_join(featureInfo(x),
                         by = c('Feature', 'MF', 'Isotope', 'Adduct')) %>% 
              select(feature = Name,kingdom:last_col(offset = 2))
            
            background_classification_totals <- feature_classifications %>% 
              gather(level,classification,-feature) %>% 
              drop_na() %>% 
              group_by(level,classification) %>% 
              count()
            
            explanatory_classifications <- explanatoryFeatures(x,...) %>% 
              select(-metric,-value) %>% 
              left_join(feature_classifications,
                        by = 'feature') 
            
            if (rf_type == 'classification'){
              explanatory_classifications <- explanatory_classifications %>% 
                group_by(response,comparison)
            }
            
            if (rf_type == 'regression'){
              explanatory_classifications <- explanatory_classifications %>% 
                group_by(response)
            }
            
            explanatory_totals <- explanatory_classifications %>% 
              count()
            
            if (rf_type == 'classification'){
              explanatory_classification_totals <- explanatory_classifications %>% 
                gather(level,classification,-response,-comparison,-feature) %>% 
                drop_na() %>% 
                group_by(response,comparison,level,classification) %>% 
                count()  
            }
            
            if (rf_type == 'regression'){
              explanatory_classification_totals <- explanatory_classifications %>% 
                gather(level,classification,-response,-feature) %>% 
                drop_na() %>% 
                group_by(response,level,classification) %>% 
                count()
            }
            
            if (rf_type == 'unsupervised'){
              explanatory_classification_totals <- explanatory_classifications %>% 
                gather(level,classification,-feature) %>% 
                drop_na() %>% 
                group_by(level,classification) %>% 
                count()
            }
            
            contingency <- explanatory_classification_totals %>% 
              rename(`Explanatory & in class` = n) %>% 
              left_join(background_classification_totals,
                        by = c('level','classification')) %>% 
              rename(`Not explanatory & in class` = n) %>% 
              mutate(`Not explanatory & in class` = `Not explanatory & in class` - 
                       `Explanatory & in class`)
            
            if (rf_type == 'classification'){
              contingency <- contingency %>%
                left_join(explanatory_totals,by = c('response','comparison'))
            }
            
            if (rf_type == 'regression'){
              contingency <- contingency %>%
                left_join(explanatory_totals,by = c('response'))
            }
            
            if (rf_type == 'unsupervised'){
              contingency <- contingency %>%
                mutate(n = explanatory_totals$n[1])
            }
            
            contingency <- contingency %>% 
              rename(`Explanatory & not in class` = n) %>% 
              mutate(`Explanatory & not in class` = `Explanatory & not in class` - 
                       `Explanatory & in class`,
                     `Not explanatory & not in class` = nFeatures(x) - 
                       (`Explanatory & in class` + `Not explanatory & in class`) - 
                       `Explanatory & not in class`)
            
            ora_results <- contingency %>% 
              rowwise() %>% 
              group_map(
                ~bind_cols(.x,
                           overrepresentation(
                             .x$`Explanatory & in class`,
                             .x$`Not explanatory & in class`,
                             .x$`Explanatory & not in class`,
                             .x$`Not explanatory & not in class`)
                )
              ) %>% 
              bind_rows()
            
            if (rf_type == 'classification'){
              ora_results <- ora_results %>% 
                group_by(response,comparison)
            }
            
            if (rf_type == 'regression'){
              ora_results <- ora_results %>% 
                group_by(response)
            }
            
            ora_results <- ora_results %>% 
              mutate(
                adjusted.p.value = p.adjust(p.value,method = p_adjust_method)
              ) %>% 
              relocate(adjusted.p.value,.after = p.value) %>%
              ungroup() %>% 
              arrange(p.value)
            
            results <- new('StructuralEnrichment',
                           x,
                           explanatory = explanatoryFeatures(x,...),
                           classifications = feature_classifications,
                           results = ora_results)
            
             return(results)
          })

#' @importFrom construction classifications
#' @importFrom dplyr ungroup

setMethod('structuralEnrichment',signature = c('RandomForest','Construction'),
          function(x,
                   structural_classifications,
                   p_adjust_method = 'bonferroni',
                   ...){
            structural_classifications <- classifications(structural_classifications)
            
            structuralEnrichment(x = x,
                                 structural_classifications = structural_classifications,
                                 p_adjust_method = p_adjust_method)
          })