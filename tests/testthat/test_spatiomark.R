## SpatioMark's C++ neighbour summaries and least-squares state changes against the reference computations they
## replaced (spatstat's closepairs() with tidyr's pivot_wider(), and limma's lmFit()), and kontextualTest().

ref_summary <- function(d, maxDist, distFun) {
  ow <- spatstat.geom::owin(range(d$x), range(d$y))
  X <- spatstat.geom::ppp(d$x, d$y, window = ow, marks = d$cellType)
  cp <- spatstat.geom::closepairs(X, rmax = maxDist, what = "ijd")
  dd <- data.frame(cellID = d$cellID[cp$i], cellType = d$cellType[cp$j], d = cp$d)
  fn <- if (distFun == "min") function(x) min(c(x, maxDist), na.rm = TRUE) else function(x) sum(c(x, 0) > 0, na.rm = TRUE)
  w <- tidyr::pivot_wider(dd, names_from = cellType, values_from = d, values_fn = fn, values_fill = fn(NULL))
  w <- dplyr::left_join(d[, "cellID", drop = FALSE], w, by = "cellID") |> tibble::column_to_rownames("cellID")
  w
}

test_that("neighbour distances and abundances match closepairs() + pivot_wider()", {
  set.seed(3)
  d <- data.frame(x = c(runif(300, 0, 500), 2000), y = c(runif(300, 0, 500), 2000),
                  cellType = c(sample(c("A", "B", "C", "rare"), 300, TRUE, c(.4, .3, .29, .01)), "A"))
  d$x[2] <- d$x[1]; d$y[2] <- d$y[1]                                    # a coincident pair (d = 0)
  d$cellID <- paste0("c", seq_len(nrow(d)))
  for (fun in c("min", "abundance")) {
    ref <- ref_summary(d, 60, fun)
    new <- Statial:::distanceCalculator(d, maxDist = 60, distFun = fun)
    expect_setequal(names(new), names(ref))
    expect_identical(rownames(new), rownames(ref))
    expect_equal(as.matrix(new[, names(ref)]), as.matrix(ref), tolerance = 1e-12)
  }
})

test_that("state changes match lmFit, with and without contamination covariates", {
  set.seed(5)
  n <- 200; markers <- paste0("m", 1:6)
  distances <- data.frame(A = runif(n, 0, 100), B = c(NA, runif(n - 1, 0, 100)), C = rep(c(1, 2), n / 2), D = rep(3, n))
  intensities <- as.data.frame(matrix(rnorm(n * 6), n, dimnames = list(NULL, markers)))
  intensities$m1 <- intensities$m1 + 0.02 * distances$A
  p <- matrix(runif(n * 3), n); contam <- as.data.frame(p / rowSums(p)); names(contam) <- c("t1", "t2", "t3")
  ref_one <- function(x, contaminations) {
    if (length(unique(x)) <= 1) return(NULL)
    design <- if (is.null(contaminations)) data.frame(coef = 1, cellType = x) else {
      ds <- data.frame(coef = 1, cellType = x, contaminations)
      nzv <- vapply(ds, function(z) var(z) > 0, logical(1)); nzv[1] <- TRUE; ds <- ds[, nzv, drop = FALSE]
      q <- qr(ds); keep <- q$pivot[seq_len(q$rank)]; ds[, c(intersect(1:2, keep), setdiff(keep, 1:2)), drop = FALSE] }
    design$cellType[is.na(design$cellType)] <- 0
    fit <- limma::lmFit(t(intensities), design)
    tv <- (fit$coef / fit$stdev.unscaled / fit$sigma)[, "cellType"]
    data.frame(marker = markers, coef = fit$coef[, "cellType"], tval = tv, pval = 2 * pt(-abs(tv), fit$df.residual))
  }
  for (ct in list(NULL, contam)) {
    con <- if (is.null(ct)) data.frame(madeUp = rep(-99, n)) else ct
    new <- Statial:::calculateChangesMarker(distances[, c("A", "C", "D")], intensities, con, 1, "g")
    ref <- dplyr::bind_rows(Filter(Negate(is.null), lapply(distances[, c("A", "C", "D")], ref_one, contaminations = ct)),
                            .id = "otherCellType")
    expect_equal(new$otherCellType, ref$otherCellType)
    expect_equal(new$coef, unname(ref$coef), tolerance = 1e-10)
    expect_equal(new$tval, unname(ref$tval), tolerance = 1e-10)
    expect_equal(new$pval, unname(ref$pval), tolerance = 1e-10)
  }
  # a column with missing values goes through the per-column fit (missing set to 0)
  new <- Statial:::calculateChangesMarker(distances[, "B", drop = FALSE], intensities, data.frame(madeUp = rep(-99, n)), 1, "g")
  expect_equal(new$coef, unname(ref_one(distances$B, NULL)$coef), tolerance = 1e-10)
})

test_that("kontextualTest returns the triples of parentDf", {
  data("kerenSCE")
  kerenSCE <- kerenSCE[, kerenSCE$imageID %in% c("1", "5", "6", "14", "18", "21")]
  kerenSCE$group <- ifelse(kerenSCE$imageID %in% c("1", "5", "6"), "a", "b")
  parentDf <- parentCombinations(all = unique(kerenSCE$cellType), parentList = list(tcells = c("CD4_Cell", "CD8_Cell")))
  res <- kontextualTest(kerenSCE, parentDf, condition = "group", r = 50, from = "Keratin_Tumour")
  tab <- res$cellResults
  expect_setequal(rownames(tab), c("Keratin_Tumour__CD4_Cell__tcells", "Keratin_Tumour__CD8_Cell__tcells"))
  expect_true(all(is.finite(tab$p_value)))
})
