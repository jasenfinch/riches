test_that("structural_enrichment works", {
    random_forest <- assigned_data %>% 
    metabolyseR::randomForest(
      cls = 'class',
      comparisons = list(
        class = 'ABR1~BD21'
      )
    )

    explanatory_features <- metabolyseR::explanatoryFeatures(random_forest)$feature

    enrichment <- structural_enrichment(
        structural_classifications,
        explanatory_features
    )

    expect_s3_class(enrichment, 'tbl_df')
})
