
test_that("organism data can be loaded", {
  data_path <- system.file('bdi',
                           package = 'riches')
  organism_data <- organismData(
    'bdi',
    database_directory = data_path,
    internal_directory = FALSE)
  
  expect_s4_class(organism_data,'FELLA.DATA')
})

test_that('organism compounds can be extracted',{
  data_path <- system.file('bdi',
                           package = 'riches')
  organism_data <- organismData(
    'bdi',
    database_directory = data_path,
    internal_directory = FALSE)
  
  organism_compounds <- organismCompounds(organism_data)
  
  expect_s3_class(organism_compounds,'tbl_df')
})