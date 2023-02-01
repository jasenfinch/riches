
test_that("structural enrichment works", {
  random_forest <- assigned_data %>% 
    metabolyseR::randomForest(
      cls = 'class'
    )
  
  structural_enrichment <- structuralEnrichment(
    random_forest,
    structural_classifications
  )
  
  expect_s4_class(structural_enrichment,'StructuralEnrichment')
  expect_output(show(structural_enrichment),'Random forest')
})
