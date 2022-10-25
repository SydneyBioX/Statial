
#' @noRd
#' 
#' @importFrom BiocParallel SerialParam SnowParam MulticoreParam bpparam
.generateBPParam <- function(cores = 1) {
    seed <- .Random.seed[1]
    
    if (cores == 1) {
        BPparam <- BiocParallel::SerialParam(RNGseed = seed)
    } else { ## Parallel processing is desired.
        ## Also set the BPparam RNGseed if the user ran set.seed(someNumber) themselves.
        if (Sys.info()["sysname"] == "Windows") { # Only SnowParam suits Windows.
            BPparam <- BiocParallel::SnowParam(min(
                cores,
                BiocParallel::snowWorkers("SOCK")
            ),
            RNGseed = seed
            )
        } else if (Sys.info()["sysname"] %in% c("MacOS", "Linux")) {
            BPparam <- BiocParallel::MulticoreParam(min(
                cores,
                BiocParallel::multicoreWorkers()
            ),
            RNGseed = seed
            )
            ## Multicore is faster than SNOW, but it doesn't work on Windows.
        } else { ## Something weird.
            BPparam <- BiocParallel::bpparam() ## BiocParallel will figure it out.
        }
    }
    
    BPparam
}