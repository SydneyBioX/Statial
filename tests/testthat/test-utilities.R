test_that(".par_lapply() and .par_mapply() give the serial results and leave BLAS threads as they were", {
  skip_on_os("windows")
  X <- setNames(as.list(1:20), letters[1:20])
  blasThreads <- RhpcBLASctl::blas_get_num_procs()
  expect_identical(.par_lapply(X, function(x, y) x * y, y = 3, cores = 2), lapply(X, function(x, y) x * y, y = 3))
  expect_identical(.par_lapply(X, sqrt, cores = 2), lapply(X, sqrt))
  expect_identical(
    .par_mapply(function(a, b, k) a + b + k, a = X, b = rev(X), MoreArgs = list(k = 1), cores = 2),
    mapply(function(a, b, k) a + b + k, a = X, b = rev(X), MoreArgs = list(k = 1), SIMPLIFY = FALSE)
  )
  expect_identical(RhpcBLASctl::blas_get_num_procs(), blasThreads)
  expect_error(.par_lapply(1:4, function(x) if (x == 3) stop("bad task") else x, cores = 2), "bad task")
  expect_identical(RhpcBLASctl::blas_get_num_procs(), blasThreads)
})

test_that("random numbers in .par_lapply() tasks depend on the seed, not on the number of cores", {
  skip_on_os("windows")
  draw <- function(cores) .par_lapply(1:20, function(x) sample(100, 3), cores = .bpparam(cores, RNGseed = 1))
  expect_identical(draw(2), draw(1))
  expect_identical(draw(3), draw(1))
})

test_that("relabelKontextual() is reproducible with set.seed()", {
  skip_on_os("windows")
  data("kerenSCE")
  cells <- kerenSCE[, kerenSCE$imageID == "6"]
  relabelled <- function(cores) {
    set.seed(51773)
    relabelKontextual(cells, nSim = 4, r = 50, from = "Macrophages", to = "Keratin_Tumour",
      parent = c("Macrophages", "CD4_Cell"), cores = cores)
  }
  expect_identical(relabelled(2), relabelled(1))
})
