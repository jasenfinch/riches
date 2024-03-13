
test_that("structural enrichment works for random forest classification", {
  random_forest <- assigned_data %>% 
    metabolyseR::randomForest(
      cls = 'class',
      comparisons = list(
        class = 'ABR1~BD21'
      )
    )
  
  structural_enrichment <- structuralEnrichment(
    random_forest,
    structural_classifications,
    split = 'trends'
  )
  
  expect_s4_class(structural_enrichment,'StructuralEnrichment')
  expect_output(show(structural_enrichment),'Random forest')
})

test_that("structural enrichment works for random forest regression", {
  random_forest <- assigned_data %>% 
    metabolyseR::randomForest(
      cls = 'injOrder'
    )
  
  structural_enrichment <- structuralEnrichment(
    random_forest,
    structural_classifications,
    split = 'trends',
    metric = '%IncMSE'
  )
  
  expect_s4_class(structural_enrichment,'StructuralEnrichment')
  expect_output(show(structural_enrichment),'Random forest')
})

test_that("structural enrichment works for unsupervised random forest", {
  random_forest <- assigned_data %>% 
    metabolyseR::randomForest(
      cls = NULL
    )
  
  structural_enrichment <- structuralEnrichment(
    random_forest,
    structural_classifications
  )
  
  expect_s4_class(structural_enrichment,'StructuralEnrichment')
  expect_output(show(structural_enrichment),'Unsupervised')
})

test_that("structural enrichment works for t-test", {
  t_test <- assigned_data %>% 
    metabolyseR::ttest(
      cls = 'class',
      comparisons = list(
        class = 'ABR1~BD21'
      )
    )
  
  structural_enrichment <- structuralEnrichment(
    t_test,
    structural_classifications,
    split = 'trends'
  )
  
  expect_s4_class(structural_enrichment,'StructuralEnrichment')
  expect_output(show(structural_enrichment),'Univariate')
})

test_that("structural enrichment works for anova", {
  anova <- assigned_data %>% 
    metabolyseR::anova(
      cls = 'class',
      comparisons = list(
        class = 'ABR1~BD21'
      )
    )
  
  structural_enrichment <- structuralEnrichment(
    anova,
    structural_classifications,
    split = 'trends'
  )
  
  expect_s4_class(structural_enrichment,'StructuralEnrichment')
  expect_output(show(structural_enrichment),'Univariate')
})

test_that("structural enrichment works for linear regression", {
  lr <- assigned_data %>% 
    metabolyseR::linearRegression(
      cls = 'injOrder',
    )
  
  structural_enrichment <- structuralEnrichment(
    lr,
    structural_classifications,
    split = 'trends',
    value = 'p.value'
  )
  
  expect_s4_class(structural_enrichment,'StructuralEnrichment')
  expect_output(show(structural_enrichment),'Univariate')
})
