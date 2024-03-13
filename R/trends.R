
trendType <- function(x){
  model_type <- type(x)
  
  switch(
    model_type,
    unsupervised = stop('Trends cannot be used for unsuperised random forest.',call. = FALSE),
    classification = 'ratio',
    regression = 'correlation',
    `t-test` = 'ratio',
    ANOVA = 'ratio',
    `linear regression` = 'correlation'
  )
}

#' @importFrom stringr str_extract_all

checkComparisons <- function(explanatory_features){
  comparisons <- explanatory_features %>% 
    select(comparison,feature) %>% 
    mutate(n_comparisons = comparison %>% 
             str_extract_all('~') %>% 
             nchar())
  
  if (any(comparisons$n_comparisons > 1)){
    non_binary <- comparisons %>% 
      select(comparison,n_comparisons) %>% 
      distinct() %>% 
      filter(n_comparisons > 1) %>% 
      .$comparison
    
    warning(paste0(
      'Trends only possible for binary comparisons. Comparisons ',
      paste(non_binary,collapse = ', '),
      ' will be removed.'),
      call. = FALSE
    )
    
    explanatory_features <- explanatory_features %>% 
      filter(!comparison %in% non_binary)
  }
  
  if (nrow(explanatory_features) == 0){
    stop('No binary comparisons available to calculate trends.',
         call. = FALSE)
  }  
  
  return(explanatory_features)
}

#' @importFrom metabolyseR sinfo keepFeatures dat
#' @importFrom rlang sym
#' @importFrom tidyr spread
#' @importFrom dplyr summarise all_of
#' @importFrom stats median

ratio <- function(x,explanatory_features){
  
  explanatory_features <- checkComparisons(explanatory_features)
  
  comparisons <- explanatory_features %>% 
    select(comparison,feature) %>% 
    separate(
      comparison,
      c('class1','class2'),
      sep = '~',
      remove = FALSE
    )
  
  explanatory_feature_data <- x %>% 
    keepFeatures(
      features = unique(explanatory_features$feature)
    ) %>% 
    {
      d <- .
      
      dat(d) %>% 
        bind_cols(
          sinfo(d) %>% 
            select(all_of(explanatory_features$response[1]))
        )
    } %>% 
    gather(feature,intensity,-!!explanatory_features$response[1]) %>% 
    group_by(!!sym(explanatory_features$response[1]),feature) %>% 
    summarise(median = median(intensity),
              .groups = 'drop') %>% 
    spread(!!explanatory_features$response[1],median)
  
  ratios <- comparisons %>% 
    left_join(explanatory_feature_data,
              by = 'feature') %>% 
    mutate(ratio = !!sym(.$class2[1]) / !!sym(.$class1[1]),
           log2_ratio = log2(ratio)) %>% 
    select(comparison,feature,ratio,log2_ratio) %>%
    mutate(
      trend = ifelse(log2_ratio < 0,'decreased','increased')
    )
  
  explanatory_features <- explanatory_features %>% 
    left_join(ratios,
              by = c('comparison',
                     'feature'))
  
  return(explanatory_features)
  
}

#' @importFrom stats cor

correlation <- function(x,explanatory_features){
  correlations <- x %>% 
    keepFeatures(
      features = unique(explanatory_features$feature)
    ) %>% 
    {
      d <- .
      
      dat(d) %>% 
        bind_cols(
          sinfo(d) %>% 
            select(all_of(explanatory_features$response[1]))
        )
    } %>% 
    gather(feature,intensity,-!!explanatory_features$response[1]) %>% 
    group_by(feature) %>% 
    summarise(correlation = cor(!!sym(explanatory_features$response[1]),
                                intensity,
                                method = 'spearman'),
              .groups = 'drop') %>% 
    mutate(
      trend = ifelse(correlation < 0,'negatively correlated','positively correlated')
    )
  
  explanatory_features %>% 
    left_join(correlations,
              by = 'feature')
}

trends <- function(x,explanatory_features){
  trend_type <- trendType(x)
  
  explanatory_features <- do.call(trend_type,
                                  list(x = x,explanatory_features = explanatory_features))
  
  return(explanatory_features)
}
