## Create progress bar (storePB removed — no longer needed without Fortran)
createPB  <- function(...) {
    utils::txtProgressBar(...)
}
