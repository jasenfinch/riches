
test_that('functional enrichment works for random forest classification',{
  random_forest <- assigned_data %>%
    metabolyseR::randomForest(
      cls = 'class',
      binary = TRUE
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
    ),
    split = 'trends'
  )
  
  expect_s4_class(enrichment_results,'FunctionalEnrichment')
  expect_s3_class(hits(enrichment_results),'tbl_df')
  expect_s3_class(explanatoryFeatures(enrichment_results),'tbl_df')
  expect_type(enrichmentResults(enrichment_results),'list')
  expect_s3_class(generateResultsTable(enrichment_results),'tbl_df')
  expect_output(show(enrichment_results),'Random forest')
})

test_that('functional enrichment works for random forest regression',{
  random_forest <- assigned_data %>%
    metabolyseR::randomForest(
      cls = 'injOrder'
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
    ),
    split = 'trends',
    metric = '%IncMSE'
  )
  
  expect_s4_class(enrichment_results,'FunctionalEnrichment')
  expect_s3_class(hits(enrichment_results),'tbl_df')
  expect_s3_class(explanatoryFeatures(enrichment_results),'tbl_df')
  expect_type(enrichmentResults(enrichment_results),'list')
  expect_s3_class(generateResultsTable(enrichment_results),'tbl_df')
  expect_output(show(enrichment_results),'Random forest')
})

test_that('functional enrichment works for unsupervised random forest',{
  random_forest <- assigned_data %>%
    metabolyseR::randomForest(
      cls = NULL
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
  expect_output(show(enrichment_results),'Unsupervised')
})

test_that('functional enrichment errors if no binary comparisons are found for trends for unsupervised random forest',{
  random_forest <- assigned_data %>%
    metabolyseR::randomForest(
      cls = 'class'
    )
  
  expect_warning(expect_error(functionalEnrichment(
    random_forest,
    'bdi',
    methods = 'hypergeom',
    organism_data = organismData(
      'bdi',
      database_directory = system.file(
        'bdi',
        package = 'riches'),
      internal_directory = FALSE
    ),
    split = 'trends')
  ))
})

test_that('functional enrichment errors with trends for unsupervised random forest',{
  random_forest <- assigned_data %>%
    metabolyseR::randomForest(
      cls = NULL
    )
  
  expect_error(functionalEnrichment(
    random_forest,
    'bdi',
    methods = 'hypergeom',
    organism_data = organismData(
      'bdi',
      database_directory = system.file(
        'bdi',
        package = 'riches'),
      internal_directory = FALSE
    ),
    split = 'trends')
  )
})
