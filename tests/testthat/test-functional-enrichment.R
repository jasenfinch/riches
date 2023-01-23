
test_that('functional enrichment works for RandomForest class',{
  parameters <- enrichmentParameters('bdi')
  functional(parameters) <- list(methods = 'diffusion')
  
  fe <- functionalEnrichment(example_analysis,example_assignment,parameters)
  
  expect_s4_class(fe,'FunctionalEnrichment')
})