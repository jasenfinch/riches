
test_that('functional enrichment works for RandomForest class',{
  random_forest <- assigned_data %>%
    metabolyseR::randomForest(
      cls = c('class')
    )
  
  enrichment_results <- functionalEnrichment(
    random_forest,
    'bdi',
    methods = 'hypergeom',
    organism_data = organismData(
      'bdi',
      database_directory = system.file(
        'bdi',
        package = 'riches'),
      internal_directory = FALSE
    )
  )
  
  expect_s4_class(enrichment_results,'FunctionalEnrichment')
  expect_s3_class(hits(enrichment_results),'tbl_df')
  expect_s3_class(explanatoryFeatures(enrichment_results),'tbl_df')
  expect_type(enrichmentResults(enrichment_results),'list')
  expect_s3_class(generateResultsTable(enrichment_results),'tbl_df')
  expect_output(show(enrichment_results),'Random forest')
})
