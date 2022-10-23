

test_that("Expecting an error when parentDf or to, from, parent is NULL", {
  data("headSCE")

  testthat::expect_error(
    Konditional(
      cells = headSCE,
      r = 50,
      parentDf = NULL,
      to = NULL,
      from = NULL,
      parent = NULL
    )
  )
})


test_that("Konditional values from headSCE are the same as the saved ones", {
  data("headSCE")

  expected_data <- readRDS("CD4_Konditional.rds")

  CD4_Konditional <-
    Konditional(
      cells = headSCE,
      r = 50,
      from = "TC_CD4",
      to = "SC5",
      parent = c("TC_CD4", "TC_CD8"),
      inhom = TRUE,
      edgeCorrect = FALSE,
      window = "convex",
      window.length = NA,
      weightQuantile = .80,
      includeZeroCells = TRUE,
      includeOriginal = TRUE
    )

  testthat::expect_equal(
    CD4_Konditional,
    expected_data
  )
})


