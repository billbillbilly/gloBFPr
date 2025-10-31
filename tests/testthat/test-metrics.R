testthat::test_that("runs correctly", {
  res_pop <- gloBFPr::get_pop_density(year = 2025)
  res_morphology <- gloBFPr::get_morphology()
  res_neighbors <- gloBFPr::get_neighbors()
  res_dng <- gloBFPr::get_dng()

  testthat::expect_type(res_pop, "NULL")
  testthat::expect_type(res_morphology, "NULL")
  testthat::expect_type(res_neighbors, "NULL")
  testthat::expect_type(res_dng, "NULL")
})
