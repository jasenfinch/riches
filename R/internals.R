#' @importFrom metabolyseR features
#' @importFrom tidyr separate

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