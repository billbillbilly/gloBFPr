testthat::test_that("runs correctly", {
  buildings <- search_3dglobdf(bbox=c(-135.459881,57.868314,-133.526287,58.882120),
                               crop = FALSE)

  testthat::expect_type(buildings, "NULL")
})
