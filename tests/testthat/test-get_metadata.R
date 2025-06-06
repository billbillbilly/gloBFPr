testthat::test_that("runs correctly", {
  metadata <- gloBFPr::get_metadata(test=TRUE)
  testthat::expect_type(metadata, "NULL")
})
