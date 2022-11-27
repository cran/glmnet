#' Make response for coxnet
#'
#' Internal function to make the response y passed to glmnet suitable
#' for coxnet (i.e. glmnet with family = "cox"). Sanity checks are performed
#' here too.
#'
#' If y is a class "Surv" object, this function returns y with no changes. If
#' y is a two-column matrix with columns named 'time' and 'status', it is
#' converted into a "Surv" object.
#'
#' @param y Response variable. Either a class "Surv" object or a two-column
#' matrix with columns named 'time' and 'status'.
#'
#' @return A class "Surv" object.
#'
#' @importFrom survival Surv
response.coxnet <- function(y) {
  if (any(is.na(y))) stop(paste0("NAs encountered in response, not allowed"))

  # if Surv object, check that it is of correct type and perform sanity checks
  # One sanity check is that it have column names. If Surv() is called with a one-column matrix for
  # time, the name is lost
  # if all good, return with no changes
  if (is.Surv(y)) {
    if (attr(y, "type") == "right") {
      if (any(y[, 1] <= 0))
          stop("Non-positive event times encountered; not permitted for Cox family")
      colnames(y) <- c("time","status")
      return(y)
    } else if (attr(y, "type") == "counting") {
      if (any(y[, 1] < 0) || any(y[, 2] <= 0))
        stop(paste("Negative start/non-positive stop times encountered;",
                   "not permitted for Cox family"))
      if (any(y[, 1] >= y[, 2]))
        stop("Some rows have start time >= stop time; not permitted")
      colnames(y) <- c("start","stop","status")
      return(y)
    } else {
      stop("cox.path() only supports 'Surv' objects of type 'right' or 'counting'")
    }
  }

  # if two-column matrix passed, make it into a Surv object
  if (!is.matrix(y) || !all(match(c("time","status"),dimnames(y)[[2]],0)))
    stop(paste0("Cox model requires a matrix with columns 'time' (>0) and ",
                "'status' (binary) as a response; a 'Surv' object suffices"),
         call. = FALSE)
  ty <- as.double(y[,"time"])
  tevent <- as.double(y[,"status"])
  if (any(ty <= 0))
    stop("negative event times encountered; not permitted for Cox family")
  yob <- Surv(ty, tevent)
  colnames(yob) <- c("time","status")
  return(yob)
}
