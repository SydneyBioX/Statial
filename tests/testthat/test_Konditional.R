test_that("Expecting an error when parentDf or to, from, parent is NULL", {
  data("kerenSCE")

  testthat::expect_error(
    Kontextual(
      cells = kerenSCE,
      r = 50,
      parentDf = NULL,
      to = NULL,
      from = NULL,
      parent = NULL
    )
  )
})


# test_that("Kontextual values from kerenSCE are the same as the saved ones", {
#   data("kerenSCE")
#
#   expected_data <- readRDS("CD4_Kontextual.rds")
#
#   CD4_Kontextual <-
#     Kontextual(
#       cells = kerenSCE,
#       r = 50,
#       from = "TC_CD4",
#       to = "SC5",
#       parent = c("TC_CD4", "TC_CD8"),
#       inhom = TRUE,
#       edgeCorrect = FALSE,
#       window = "convex",
#       window.length = NA,
#       weightQuantile = .80,
#       includeZeroCells = TRUE,
#       includeOriginal = TRUE
#     )
#
#   testthat::expect_equal(
#     CD4_Kontextual,
#     expected_data
#   )
# })


test_that("Expecting an error when dataset contains wrong column names", {
  data("kerenSCE")

  renamedSCE <- colData(kerenSCE) |>
    data.frame() |>
    rename("imageName" = "imageID")

  testthat::expect_error(
    suppressWarnings(Kontextual(
      cells = renamedSCE,
      r = 50,
      from = "TC_CD4",
      to = "SC5",
      parent = c("TC_CD4", "TC_CD8")
    ))
  )
})


test_that("Fail on invalid input", {
  wrong_input <- c("wrong", "input")

  testthat::expect_error(
    Kontextual(
      cells = wrong_input,
      r = 50,
      from = "TC_CD4",
      to = "SC5",
      parent = c("TC_CD4", "TC_CD8")
    )
  )
})
