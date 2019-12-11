## Create and store the pb
createPB  <- function(...) {
    pb  <- utils::txtProgressBar(...)
    .Call("storePB", pb, PACKAGE = "glmnet")
    pb
}
