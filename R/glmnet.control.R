#' Internal glmnet algorithm parameters
#'
#' View and/or change the factory default parameters that control
#' glmnet's numerical algorithms.
#'
#' If called with no arguments, \code{glmnet.control()} returns a list
#' with the current settings of these parameters. Any arguments
#' included in the call set those parameters to the new values, and
#' then silently return the updated settings. Values set via
#' \code{glmnet.control()} persist for the duration of the R session.
#' All parameters listed here can also be supplied per-call via the
#' \code{control = list(...)} argument of \code{\link{glmnet}()}, which
#' does \emph{not} mutate the session state; see that function's
#' documentation for the precedence rules.
#'
#' @section Parameter taxonomy:
#' glmnet has two execution paths, and not every control parameter is
#' consumed by both:
#'
#' \itemize{
#'   \item \strong{Core-engine path} --- \code{family} passed as a
#'     character string (\code{"gaussian"}, \code{"binomial"},
#'     \code{"poisson"}, \code{"multinomial"}, \code{"mgaussian"},
#'     \code{"cox"}) routes to purpose-built C++ glmnetpp kernels.
#'   \item \strong{GLM-family (R-IRLS) path} --- \code{family} passed
#'     as a \code{family()} object (or \code{glmnet.path()} /
#'     \code{glmnet.fit()} called directly) runs an IRLS loop
#'     implemented in R that uses the C++ \code{wls} kernel as its
#'     inner solver.
#' }
#'
#' The table below shows the role and the scope of each parameter.
#' "Scope: engine" means the parameter is only meaningful on the
#' core-engine path; "Scope: R-IRLS" means only on the GLM-family path;
#' "Scope: both" means consulted by both.
#'
#' \tabular{llll}{
#'   \strong{Parameter} \tab \strong{Role} \tab \strong{Scope} \tab \strong{Notes} \cr
#'   \code{fdev}   \tab path termination         \tab both    \tab C++ name \code{sml}; min fractional deviance change \cr
#'   \code{devmax} \tab path termination         \tab both    \tab C++ name \code{rsqmax}; max explained-deviance ratio \cr
#'   \code{mnlam}  \tab path termination         \tab both    \tab min lambda count before early stop \cr
#'   \code{eps}    \tab path termination         \tab both    \tab lambda.min.ratio floor \cr
#'   \code{thresh} \tab CD tolerance             \tab both    \tab inner coordinate-descent convergence \cr
#'   \code{maxit}  \tab CD budget                \tab both    \tab max CD passes across all lambdas \cr
#'   \code{dfmax}  \tab coefficient-count cap    \tab both    \tab default resolves to \code{nvars + 1} at call time \cr
#'   \code{pmax}   \tab coefficient-count cap    \tab both    \tab default resolves to \code{min(dfmax*2+20, nvars)} at call time \cr
#'   \code{big}    \tab numerical guard          \tab both    \tab \code{Inf} substitute for bounds limits \cr
#'   \code{itrace} \tab progress                 \tab both    \tab aliased with \code{trace.it} \cr
#'   \code{trace.it} \tab progress               \tab both    \tab aliased with \code{itrace} \cr
#'   \code{pmin}   \tab probability floor        \tab engine  \tab logistic-family kernels only (\code{lognet}, \code{multnet}) \cr
#'   \code{exmx}   \tab exp() cap                \tab engine  \tab logistic-family kernels only \cr
#'   \code{prec}   \tab bounds-subsolver tolerance \tab engine  \tab active when a coefficient has finite \code{lower} / \code{upper} \cr
#'   \code{mxit}   \tab bounds-subsolver budget  \tab engine  \tab paired with \code{prec} via C++ \code{chg_bnorm()} \cr
#'   \code{epsnr}  \tab Newton-Raphson tolerance \tab R-IRLS  \tab \emph{inert on the core-engine path}; only \code{glmnet.fit()} reads it \cr
#'   \code{mxitnr} \tab Newton-Raphson budget    \tab R-IRLS  \tab \emph{inert on the core-engine path}; only \code{glmnet.fit()} reads it \cr
#' }
#'
#' The naming zoo (three tolerances, three iteration budgets) reflects
#' the fact that glmnet contains three nested loops:
#' \itemize{
#'   \item \strong{Outer path}: iterate over lambda values. Early
#'     termination governed by \code{fdev}, \code{devmax}, \code{mnlam},
#'     \code{eps}.
#'   \item \strong{Middle Newton/IRLS loop} (only on the R-IRLS path):
#'     convergence by \code{epsnr}, budget by \code{mxitnr}.
#'   \item \strong{Inner coordinate descent}: convergence by
#'     \code{thresh}, budget by \code{maxit}.
#'   \item \strong{Bounds subsolver} (only for coefficients with finite
#'     lower/upper): convergence by \code{prec}, budget by \code{mxit}.
#' }
#'
#' @param fdev Path-termination parameter: minimum fractional change
#'   in deviance between consecutive lambdas for the path to continue.
#'   Scope: both paths. Factory default \code{1.0e-5}.
#' @param devmax Path-termination parameter: maximum fraction of null
#'   deviance the path is allowed to explain before stopping. Scope:
#'   both paths. Factory default \code{0.999}.
#' @param eps Path-termination parameter: floor on
#'   \code{lambda.min.ratio} (smallest allowed \eqn{\lambda_{\min} /
#'   \lambda_{\max}}). Scope: both paths. Factory default
#'   \code{1.0e-6}.
#' @param big Numerical guard: finite substitute for \code{Inf} in
#'   \code{lower.limits} / \code{upper.limits}. Scope: both paths.
#'   Factory default \code{9.9e35}.
#' @param mnlam Path-termination parameter: minimum number of lambda
#'   values to compute before early stopping is allowed. Scope: both
#'   paths. Factory default \code{5}.
#' @param pmin Numerical guard: floor on predicted probability in
#'   logistic-family kernels (prevents \code{log(0)}). Implies a max
#'   of \code{1 - pmin}. Scope: engine only (logistic-family kernels).
#'   Factory default \code{1.0e-9}.
#' @param exmx Numerical guard: cap on \code{exp()} arguments in
#'   logistic-family kernels (prevents overflow). Scope: engine only
#'   (logistic-family kernels). Factory default \code{250}.
#' @param prec Bounds-subsolver convergence tolerance. Used by the
#'   engine's \code{bnorm()} subroutine when a coefficient has a
#'   finite \code{lower} or \code{upper} bound. Scope: engine only
#'   (bounds-constrained fits). Factory default \code{1.0e-10}.
#' @param mxit Bounds-subsolver iteration budget. Paired with
#'   \code{prec}. Scope: engine only (bounds-constrained fits).
#'   Factory default \code{100}.
#' @param itrace Progress-bar flag. \code{1} enables the progress bar
#'   in \code{glmnet()} and \code{cv.glmnet()}. Scope: both paths.
#'   Aliased with \code{trace.it} (setting one sets the other).
#'   Factory default \code{0}.
#' @param epsnr Newton-Raphson convergence tolerance for the R-level
#'   IRLS loop in \code{glmnet.fit()}. \strong{Scope: R-IRLS path
#'   only} --- the core-engine kernels do not read this parameter.
#'   Factory default \code{1.0e-6}.
#' @param mxitnr Newton-Raphson iteration budget for the R-level IRLS
#'   loop in \code{glmnet.fit()}. \strong{Scope: R-IRLS path only}.
#'   Factory default \code{25}.
#' @param thresh Inner coordinate-descent convergence threshold. Each
#'   inner CD loop continues until the maximum change in the objective
#'   after any coefficient update is less than \code{thresh} times the
#'   null deviance. Scope: both paths. Factory default \code{1e-7}.
#'   Can also be set per-call via \code{control = list(thresh = ...)}
#'   in \code{glmnet()}.
#' @param maxit Inner coordinate-descent iteration budget: maximum
#'   total passes over the data across all lambda values. Scope: both
#'   paths. Factory default \code{100000}. Can also be set per-call
#'   via \code{control = list(maxit = ...)}.
#' @param dfmax Coefficient-count cap: limit the maximum number of
#'   variables in the model at any lambda. Scope: both paths. Factory
#'   default \code{NULL}, which resolves at call time to \code{nvars +
#'   1} (i.e., unconstrained). Can also be set per-call via
#'   \code{control = list(dfmax = ...)}.
#' @param pmax Coefficient-count cap: limit the maximum number of
#'   variables ever to be nonzero anywhere on the path. Scope: both
#'   paths. Factory default \code{NULL}, which resolves at call time
#'   to \code{min(dfmax * 2 + 20, nvars)}. Can also be set per-call
#'   via \code{control = list(pmax = ...)}.
#' @param trace.it Progress-bar flag (alias of \code{itrace}). Scope:
#'   both paths. Factory default \code{0}.
#' @param factory If \code{TRUE}, reset all parameters to their
#'   factory defaults. Default is \code{FALSE}.
#' @return A list with named elements as in the argument list.
#' @author Jerome Friedman, Kenneth Tay, Trevor Hastie\cr
#'   Maintainer: Trevor Hastie \email{hastie@@stanford.edu}
#' @seealso \code{\link{glmnet}}
#' @keywords models regression
#' @examples
#'
#' glmnet.control(fdev = 0)  # continue along path even though not much changes
#' glmnet.control(thresh = 1e-8)  # tighten CD convergence for session
#' glmnet.control()  # view current settings
#' glmnet.control(factory = TRUE)  # reset all the parameters to their default
#'
#' @export glmnet.control
## Factory defaults for all 17 algorithm-control parameters. Used by
## glmnet.control(factory = TRUE).
.GLMNET_CONTROL_FACTORY <- list(
    fdev     = 1e-5,
    eps      = 1e-6,
    big      = 9.9e35,
    mnlam    = 5L,
    devmax   = 0.999,
    pmin     = 1e-9,
    exmx     = 250,
    itrace   = 0L,
    prec     = 1e-10,
    mxit     = 100L,
    epsnr    = 1e-6,
    mxitnr   = 25L,
    thresh   = 1e-7,
    maxit    = 100000L,
    dfmax    = NULL,
    pmax     = NULL,
    trace.it = 0L
)

glmnet.control <-
  function (fdev = 1e-05, devmax = 0.999, eps = 1e-06, big = 9.9e+35,
            mnlam = 5, pmin = 1e-09, exmx = 250, prec = 1e-10, mxit = 100,
            itrace = 0, epsnr = 1e-06, mxitnr = 25,
            thresh = 1e-7, maxit = 100000, dfmax = NULL,
            pmax = NULL, trace.it = 0,
            factory = FALSE)
{
    inquiry <- !nargs()

    if (factory) {
        glmnet_control_set(.GLMNET_CONTROL_FACTORY)
    } else {
        ## Build update list from supplied args only. Anything left at
        ## its formal default is left alone. Note this also fixes a
        ## pre-existing 5.0 bug: setting `prec` alone would silently
        ## reset `mxit` to 100 because `chg_bnorm` paired them.
        updates <- list()
        if (!missing(fdev))     updates$fdev     <- as.double(fdev)
        if (!missing(devmax))   updates$devmax   <- as.double(devmax)
        if (!missing(eps))      updates$eps      <- as.double(eps)
        if (!missing(big))      updates$big      <- as.double(big)
        if (!missing(mnlam))    updates$mnlam    <- as.integer(mnlam)
        if (!missing(pmin))     updates$pmin     <- as.double(pmin)
        if (!missing(exmx))     updates$exmx     <- as.double(exmx)
        if (!missing(prec))     updates$prec     <- as.double(prec)
        if (!missing(mxit))     updates$mxit     <- as.integer(mxit)
        if (!missing(itrace))   updates$itrace   <- as.integer(itrace)
        if (!missing(epsnr))    updates$epsnr    <- as.double(epsnr)
        if (!missing(mxitnr))   updates$mxitnr   <- as.integer(mxitnr)
        if (!missing(thresh))   updates$thresh   <- as.double(thresh)
        if (!missing(maxit))    updates$maxit    <- as.integer(maxit)
        if (!missing(dfmax))    updates$dfmax    <- if (is.null(dfmax)) NULL else as.integer(dfmax)
        if (!missing(pmax))     updates$pmax     <- if (is.null(pmax))  NULL else as.integer(pmax)
        if (!missing(trace.it)) updates$trace.it <- as.integer(trace.it)
        if (length(updates)) glmnet_control_set(updates)
    }

    ## Single read-back of the C++ state. dfmax/pmax come back as
    ## NA_integer_ when unset; map to NULL to preserve the historical
    ## R-side return shape. Use `["k"] <- list(NULL)` rather than
    ## `$k <- NULL` because the latter deletes the slot.
    value <- glmnet_control_get()
    if (is.na(value$dfmax)) value["dfmax"] <- list(NULL)
    if (is.na(value$pmax))  value["pmax"]  <- list(NULL)

    if (inquiry) value else invisible(value)
}

## Set of deprecation messages already emitted this session, so
## .resolve_control() emits each at most once per arg name. (Without
## this, downstream packages that wrap glmnet() in tryCatch(warning =
## ...) -- notably stm -- treat each fit as a failure and grow nlambda
## without bound.)
.deprecation_emitted <- new.env(parent = emptyenv())

## Canonical list of valid control keys (the 17 names returned by
## glmnet.control()). Used by the test helper .control_with() in
## tests/testthat/helper-functions.R. The package itself no longer
## validates control= keys -- glmnet_control_set() silently drops
## unknown names by design.
.VALID_CONTROL_KEYS <- c(
    "fdev", "eps", "big", "mnlam", "devmax", "pmin", "exmx", "itrace",
    "prec", "mxit", "epsnr", "mxitnr",
    "thresh", "maxit", "dfmax", "pmax", "trace.it"
)

## Resolve the layered algorithm-control configuration for one glmnet()
## invocation. Precedence (low to high): session defaults -> control=
## list -> deprecated top-level args.
##
## Arguments:
##   control    -- the user's control= list (may be empty/NULL).
##   nvars      -- number of predictors; needed to resolve dfmax/pmax
##                 data-dependent defaults to integers.
##   deprecated -- named list of any deprecated top-level args the
##                 caller supplied (NULL entries mean "not supplied").
##                 Each non-NULL entry triggers a once-per-session
##                 .Deprecated() warning and is folded into control.
##
## Returns a list:
##   control -- 17-key resolved list (dfmax/pmax as integers).
##   restore -- niladic function to revert C++ session state, or NULL
##              if no overrides were applied. Caller registers via
##              on.exit().
.resolve_control <- function(control, nvars, deprecated = list()) {
    ## Fold deprecated top-level args into control= with one warning
    ## per arg per session. Deprecated args win over control= entries.
    deprecated <- deprecated[!vapply(deprecated, is.null, logical(1))]
    for (nm in names(deprecated)) {
        if (!isTRUE(.deprecation_emitted[[nm]])) {
            .deprecation_emitted[[nm]] <- TRUE
            .Deprecated(msg = sprintf(
                "Passing '%s' to glmnet() is deprecated. Use control = list(%s = ...) instead.",
                nm, nm))
        }
    }
    control <- modifyList(if (is.null(control)) list() else control, deprecated)

    ## Snapshot session state. If the caller supplied any overrides,
    ## push them to the C++ globals in one shot and register a single
    ## restore thunk. The C++ setter silently ignores unknown keys, so
    ## no R-side validation is needed.
    gc <- glmnet.control()
    if (length(control)) {
        glmnet_control_set(control)
        restore <- function() glmnet_control_set(gc)
        ## Build the post-override view by overlaying the known control
        ## keys onto gc. Avoids a second glmnet.control() round-trip.
        valid    <- intersect(names(control), .VALID_CONTROL_KEYS)
        resolved <- modifyList(gc, control[valid])
    } else {
        restore  <- NULL
        resolved <- gc
    }

    ## Resolve dfmax/pmax data-dependent defaults to integers.
    resolved$dfmax <- if (is.null(resolved$dfmax)) as.integer(nvars + 1L)
                      else                          as.integer(resolved$dfmax)
    resolved$pmax  <- if (is.null(resolved$pmax))  as.integer(min(resolved$dfmax * 2L + 20L, nvars))
                      else                          as.integer(resolved$pmax)

    list(control = resolved, restore = restore)
}
