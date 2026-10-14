## Stops with an install hint if suggested package `pkg` is missing.
.need <- function(pkg, why = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is needed", if (!is.null(why)) paste0(" ", why) else "",
         ". Install it with BiocManager::install(\"", pkg, "\").", call. = FALSE)
  }
  invisible(TRUE)
}

## A BiocParallel param for `cores` workers: forked workers (MulticoreParam), or socket workers (SnowParam) on
## Windows. A param object is returned unchanged. `RNGseed` starts the random number streams of the tasks.
.bpparam <- function(cores = 1, RNGseed = NULL) {
  if (is(cores, "BiocParallelParam")) return(cores)
  if (!is.numeric(cores)) stop("`cores` must be a number (or a BiocParallel param object).", call. = FALSE)
  cores <- max(1L, as.integer(cores[1]))
  if (cores == 1L) return(BiocParallel::SerialParam(RNGseed = RNGseed))
  if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(cores, RNGseed = RNGseed)
  } else {
    BiocParallel::MulticoreParam(cores, RNGseed = RNGseed)
  }
}

## The task function and elements of the current .par_lapply() call, found by forked workers in their copy of
## this environment.
.parallel_env <- new.env(parent = emptyenv())
.run_task <- function(i) .parallel_env$task(.parallel_env$X[[i]])

## lapply(X, FUN, ...) over the workers of `cores` (a number or a BiocParallel param).
## - Forked workers start once, after the task exists, so they share the session's data without being sent it,
##   and only task numbers and results pass between processes. Forking a new process for every chunk
##   instead copies much of the session into each one.
## - BLAS runs one thread in each worker during the call. OpenBLAS otherwise starts one thread per core in every
##   worker (128 on a 128-core server), which slows the workers down many times. The setting is changed in this
##   session before the workers start, because changing it inside a worker first restarts every thread.
## - Tasks are handed out in about eight chunks per worker, so that quick and slow tasks balance out.
## - Every element has its own random number stream from the param's RNGseed, so results with a seed are the
##   same for any number of workers, including one.
.par_lapply <- function(X, FUN, ..., cores = 1) {
  BPPARAM <- .bpparam(cores)
  workers <- BiocParallel::bpnworkers(BPPARAM)
  task <- FUN
  if (...length()) {
    # The task's environment holds only FUN and the extra arguments, not X, so that socket workers aren't
    # sent all of X with every chunk.
    taskEnv <- list2env(list(FUN = FUN, args = list(...)), parent = asNamespace("Statial"))
    task <- function(x) do.call(FUN, c(list(x), args))
    environment(task) <- taskEnv
  }
  if (workers <= 1L || length(X) <= 1L) {
    if (!is(BPPARAM, "SerialParam")) {
      BPPARAM <- BiocParallel::SerialParam(RNGseed = BiocParallel::bpRNGseed(BPPARAM))
    }
    return(BiocParallel::bplapply(X, task, BPPARAM = BPPARAM))
  }

  # Only when BLAS has several threads: restoring them after forking restarts every thread, and the new
  # threads spin for a moment, so a single-threaded BLAS is left alone.
  blasThreads <- RhpcBLASctl::blas_get_num_procs()
  if (isTRUE(blasThreads > 1L)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    on.exit(RhpcBLASctl::blas_set_num_threads(blasThreads), add = TRUE)
  }
  # OpenMP likewise: forked workers can hang if this session has used several OpenMP threads.
  ompThreads <- RhpcBLASctl::omp_get_max_threads()
  if (isTRUE(ompThreads > 1L)) {
    RhpcBLASctl::omp_set_num_threads(1L)
    on.exit(RhpcBLASctl::omp_set_num_threads(ompThreads), add = TRUE)
  }
  if (BiocParallel::bpisup(BPPARAM) || !is(BPPARAM, "MulticoreParam")) {
    # Workers that are already running, or are new R sessions, are sent the task with each chunk.
    threadVars <- Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"), unset = NA)
    Sys.setenv(OPENBLAS_NUM_THREADS = 1L, OMP_NUM_THREADS = 1L)
    on.exit(.restore_env(threadVars), add = TRUE)
    return(BiocParallel::bplapply(X, task, BPPARAM = BPPARAM))
  }

  assign("task", task, envir = .parallel_env)
  assign("X", X, envir = .parallel_env)
  on.exit(rm(list = c("task", "X"), envir = .parallel_env), add = TRUE)
  tasks <- BiocParallel::bptasks(BPPARAM)
  if (tasks == 0L) {
    BiocParallel::bptasks(BPPARAM) <- min(length(X), 8L * workers)
    on.exit(BiocParallel::bptasks(BPPARAM) <- tasks, add = TRUE, after = FALSE)
  }
  BiocParallel::bpstart(BPPARAM)
  on.exit(BiocParallel::bpstop(BPPARAM), add = TRUE, after = FALSE)
  out <- BiocParallel::bplapply(seq_along(X), .run_task, BPPARAM = BPPARAM)
  names(out) <- names(X)
  out
}

## mapply(FUN, ..., MoreArgs = MoreArgs, SIMPLIFY = FALSE) over the workers of `cores`.
.par_mapply <- function(FUN, ..., MoreArgs = NULL, cores = 1) {
  args <- list(...)
  n <- max(lengths(args))
  out <- .par_lapply(seq_len(n), function(i) {
    do.call(FUN, c(lapply(args, function(a) a[[(i - 1L) %% length(a) + 1L]]), MoreArgs))
  }, cores = cores)
  names(out) <- if (!is.null(names(args[[1]]))) names(args[[1]]) else if (is.character(args[[1]])) args[[1]]
  out
}

## Put environment variables back as they were (NA: unset).
.restore_env <- function(values) {
  for (name in names(values)) {
    if (is.na(values[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(values[name]))
  }
}
