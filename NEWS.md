# glmnet 5.0 Series

## 5.0 (2026-04-27)

First release since v4.1-10 (July 2025). Major rearchitecture with
Cox speedups, a new per-call algorithm-control mechanism, and two new
vignettes.

### Architecture
* All Fortran removed; glmnetpp and coxdev are now git submodules.
* Cox is a first-class GLM type in glmnetpp; Cox deviance handled by
  coxdev (stratified and unstratified).
* Cox speedups vs CRAN 4.1-10 on current hardware: ~2× dense, 5–11×
  sparse, 10–70× stratified. Non-Cox families unchanged and
  bit-identical.

### Cox interface
* New `cox.ties = c("breslow", "efron")`;
  default is `"breslow"` for now. Calls with `family = "cox"` that
  don't set `cox.ties` emit a warning — the default switches to
  `"efron"` in v5.1 to match `survival::coxph`.

### Algorithm control
* New `control = list(...)` argument to `glmnet()` (and to
  `cv.glmnet`, `glmnet.path`, `glmnet.fit`, `elnet.fit`) overrides
  any of 17 algorithm-control parameters for a single call without
  mutating session state. Unknown keys error.
* Top-level `thresh`, `maxit`, `dfmax`, `pmax`, `trace.it` in
  `glmnet()` are deprecated; use `control = list(...)` or
  `glmnet.control()`.
* `dfmax` and `pmax` signature defaults are now `NULL`, resolving at
  call time to `nvars + 1` and `min(2*dfmax + 20, nvars)`.

### Vignettes
* New `glmnet-history.Rmd` — package history, milestones,
  attributions.
* New `coxdev.Rmd` — the coxdev deviance/gradient library.

### Other behavior changes
* `predict.glmnet(type = "nonzero", s = ...)` now always returns a
  list (one slot per `s` value, each slot an integer vector of
  active-coefficient indices). Earlier versions returned a
  `data.frame` when the per-`s` index counts happened to be uniform,
  and a list otherwise. This makes the return shape consistent with
  the no-`s` form and with the function's documentation, but is a
  silent breaking change for callers that branched on
  `is.data.frame()` of the result. To get a stable union of indices
  across `s` values regardless of glmnet version, use
  `sort(unique(unlist(predict(fit, type = "nonzero", s = lambda))))`.

# glmnet 4.1 Series

## 4.1-10 (2025-07-17)

Adjusted plotting again. Default is now using -log(lambda) as default
for plotting, which is option xvar='lambda'. Also fixed labeling
issues with predict and coef

## 4.1-9 (2025-06-02)

Changed default plotting xvar = "lambda". Corrected some minor issues
in citation file.

## 4.1-7 (2023-03-23)

Added DOI for JSS 2023 paper and corrected some typos in documentation
(`nfold` -> `nfolds`) and vignette.

## 4.1-6 (2022-11-27)

Removed unneeded legacy fortran code, leaving only coxnet. Fixed up
Matrix as() sequences

## 4.1-5

Relatively minor changes to bugs in survival functions and bigGlm,
and some improved failure messages.

## 4.1-4 (2022-04-15)

Most of the Fortran code has been replaced by C++ by James Yang,
leading to speedups in all cases. The exception is the Cox routine for
right censored data, which is still under development.

## 4.1-3 (2021-11-02)

Some of the Fortran in glmnet has been replaced by C++, written by the
newest member of our team, James Yang.
* the `wls` routines (dense and sparse), that are the engines under the `glmnet.path`
  function when we use programmable families, are now written in C++, and lead
  to speedups of around 8x.
* the family of elnet routines (sparse/dense, covariance/naive) for
  `glmnet(...,family="gaussian")` are all in C++, and lead to speedups
  around 4x.

## 4.1-2 (2021-06-24)

A new feature added, as well as some minor fixes to documentation.
* The exclude argument has come to life. Users can now pass a function
  that can take arguments x, y and weights, or a subset of these, for
  filtering variables. Details in documentation and vignette.
* Prediction with single `newx` observation failed before. This is
  fixed.
* Labeling of predictions from `cv.glmnet` improved.
* Fixed a bug in mortran/fortran that caused program to loop ad infinitum

## 4.1-1 (2021-02-21)

Fixed some bugs in the coxpath function to do with sparse X.
* when some penalty factors are zero, and X is sparse, we should not
  call GLM to get the start
* apply does not work as intended with sparse X, so we now use matrix
  multiplies instead in computing lambda_max
* added documentation for `cv.glmnet` to explain implications of
  supplying `lambda`

## 4.1 (2021-01-11)

Expanded scope for the Cox model.
* We now allow (start, stop) data in
  addition to the original right-censored all start at zero option.
* Allow for strata as in `survival::coxph`
* Allow for sparse X matrix with Cox models (was not available before)
* Provide method for `survival::survfit`

Vignettes are revised and reorganized.
Additional index information stored on `cv.glmnet` objects, and
included when printed.

# glmnet 4.0 Series

## 4.0-2 (2020-06-16)

* Biggest change. Cindex and auc calculations now use the `concordance`
  function from package `survival`
* Minor changes. Allow coefficient warm starts for glmnet.fit. The print
  method for glmnet now really prints %Dev rather than the fraction.

## 4.0 (2020-05-14)

Major revision with added functionality. Any GLM family can be used
now with `glmnet`, not just the built-in families. By passing a
"family" object as the family argument (rather than a character
string), one gets access to all families supported by `glm`. This
development was programmed by our newest member of the `glmnet` team,
Kenneth Tay.

# glmnet 3.0 Series

## 3.0-3

Bug fixes

* `Intercept=FALSE` with "Gaussian" is fixed. The `dev.ratio` comes out
  correctly now. The mortran code was changed directly in 4
  places. look for "standard". Thanks to Kenneth Tay.

## 3.0-2 (2019-12-11)

Bug fixes

* `confusion.glmnet` was sometimes not returning a list because of
  apply collapsing structure
* `cv.mrelnet` and `cv.multnet` dropping dimensions inappropriately
* Fix to `storePB` to avoid segfault. Thanks [Tomas  Kalibera](https://github.com/kalibera)!
* Changed the help for `assess.glmnet` and cousins to be more helpful!
* Changed some logic in `lambda.interp` to avoid edge cases  (thanks
  David Keplinger)

## 3.0-1 (2019-11-15)

Minor fix to correct Depends in the DESCRIPTION to R (>= 3.6.0)

## 3.0 (2019-11-09)

 This is a major revision with much added functionality, listed
 roughly in order of importance. An additional vignette called `relax`
 is supplied to describe the usage.

* `relax` argument added to `glmnet`. This causes the
   models in the path to be refit without regularization. The
   resulting object inherits from class `glmnet`, and has an
   additional component, itself a glmnet object, which is the relaxed
   fit.
* `relax` argument to `cv.glmnet`. This allows selection from a
   mixture of the relaxed fit and the regular fit. The mixture is
   governed by an argument `gamma` with a default of 5 values between
   0 and 1.
* `predict`, `coef` and `plot` methods for `relaxed` and `cv.relaxed`
   objects.
* `print` method for `relaxed` object, and new `print` methods for
   `cv.glmnet` and `cv.relaxed` objects.
* A progress bar is provided via an additional argument
   `trace.it=TRUE` to `glmnet` and `cv.glmnet`. This can also be set
   for the session via `glmnet.control`.
* Three new functions `assess.glmnet`, `roc.glmnet` and
   `confusion.glmnet` for displaying the performance of models.
* `makeX` for building the `x` matrix for input to `glmnet`. Main
   functionality is *one-hot-encoding* of factor variables, treatment
   of `NA` and creating sparse inputs.
* `bigGlm` for fitting the GLMs of `glmnet` unpenalized.

In addition to these new features, some of the code in `glmnet` has
been tidied up, especially related to CV.

# glmnet 2.0 Series

## 2.0-20

* Fixed a bug in internal function `coxnet.deviance` to do with input
  `pred`, as well as saturated `loglike` (missing) and weights
* added a `coxgrad` function for computing the gradient

## 2.0-19

* Fixed a bug in coxnet to do with ties between death set and risk set

## 2.0-18 (2019-05-20)

* Added an option alignment to `cv.glmnet`, for cases when wierd
  things happen

## 2.0-17

* Further fixes to mortran to get clean fortran; current mortran src is in
  `inst/mortran`

## 2.0-16 (2018-04-02)

* Additional fixes to mortran; current mortran src is in
  `inst/mortran`
* Mortran uses double precision, and variables are initialized to
  avoid `-Wall` warnings
* cleaned up repeat code in CV by creating a utility function

## 2.0-15

* Fixed up the mortran so that generic fortran compiler can
  run without any configure

## 2.0-13 (2017-09-22)

* Cleaned up some bugs to do with exact prediction
* `newoffset` created problems all over - fixed these

## 2.0-11

* Added protection with `exact=TRUE` calls to `coef` and
  `predict`. See help file for more details

## 2.0-10 (2017-05-06)

* Two iterations to fix to fix native fortran registration.

## 2.0-8 (2017-04-30)

* included native registration of fortran

## 2.0-7

* constant `y` blows up `elnet`; error trap included
* fixed `lambda.interp` which was returning `NaN` under degenerate
  circumstances.

## 2.0-6

* added some code to extract time and status gracefully from a `Surv` object

## 2.0-3 (2016-02-23)

* changed the usage of `predict` and `coef` with `exact=TRUE`. The
  user is strongly encouraged to supply the original `x` and `y`
  values, as well as any other data such as weights that were used in
  the original fit.

## 2.0-1 (2015-04-08)

* Major upgrade to CV; let each model use its own lambdas, then predict at original set.
* fixed some minor bugs

# glmnet 1.0 Series

## 1.9-9

* fixed subsetting bug in `lognet` when some weights are zero and `x` is sparse

## 1.9-8 (2014-05-24)

* fixed bug in multivariate response model (uninitialized variable), leading to valgrind issues
* fixed issue with multinomial response matrix and zeros
* Added a link to a glmnet vignette

## 1.9-6

* fixed bug in `predict.glmnet`, `predict.multnet` and `predict.coxnet`,
  when `s=` argument is used with a vector of values. It was not doing
  the matrix multiply correctly
* changed documentation of glmnet to explain logistic response matrix

## 1.9-5 (2013-08-04)

* added parallel capabilities, and fixed some minor bugs

## 1.9-3 (2013-03-02)

* added `intercept` option

## 1.9-1 (2013-02-10)

* added upper and lower bounds for coefficients
* added `glmnet.control` for setting systems parameters
* fixed serious bug in `coxnet`

## 1.8-5 (2013-01-04)

* added `exact=TRUE` option for prediction and coef functions

## 1.8 (2012-07-03)

* Major new release
* added `mgaussian` family for multivariate response
* added `grouped` option for multinomial family

## 1.7-4 (2012-04-27)

* nasty bug fixed in fortran - removed reference to dble
* check class of `newx` and make `dgCmatrix` if sparse

## 1.7-1 (2011-09-23)

* `lognet` added a classnames component to the object
* `predict.lognet(type="class")` now returns a character vector/matrix

## 1.6 (2011-04-24)

* `predict.glmnet` : fixed bug with `type="nonzero"`
* `glmnet`: Now x can inherit from `sparseMatrix` rather than the very
   specific `dgCMatrix`, and this will trigger sparse mode for glmnet

## 1.5 (2010-11-04)

* `glmnet.Rd` (`lambda.min`) : changed value to 0.01 if `nobs < nvars`,
  (`lambda`) added warnings to avoid single value,
  (`lambda.min`): renamed it `lambda.min.ratio`
* `glmnet` (`lambda.min`) : changed	value to 0.01 if `nobs < nvars`
    (`HessianExact`) : changed the sense (it was wrong),
   (`lambda.min`): renamed it `lambda.min.ratio`. This allows it to be
    called `lambda.min` in a call though
* `predict.cv.glmnet` (new function) : makes predictions directly from
  the saved `glmnet` object on the cv object
* `coef.cv.glmnet` (new function) : as above
* `predict.cv.glmnet.Rd` : help functions for the above
* `cv.glmnet` : insert `drop(y)` to avoid 1 column matrices; now include a `glmnet.fit` object for later predictions
* `nonzeroCoef` : added a special case for a single variable in `x`; it was dying on this
* `deviance.glmnet` : included
* `deviance.glmnet.Rd` : included

## 1.4 (2010-06-16)

* Note that this starts from version `glmnet_1.4`.
