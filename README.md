
<!-- README.md is generated from the source: README.Rmd -->

# Lasso and Elastic-Net Regularized Generalized Linear Models <img src="man/figures/logo.png" width="100" align="right" />

<!-- badges: start -->
<!-- NOTE on badges below: Manually comment out R_CMD_check badge as -->
<!-- the repo is private and will frail CRAN README.md checks -->
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/glmnet)](https://cran.r-project.org/package=glmnet)[![](https://cranlogs.r-pkg.org/badges/glmnet)](https://CRAN.R-project.org/package=glmnet)
<!-- badges: end -->

We provide extremely efficient procedures for fitting the entire lasso
or elastic-net regularization path for linear regression (gaussian),
multi-task gaussian, logistic and multinomial regression models (grouped
or not), Poisson regression and the Cox model. The algorithm uses
cyclical coordinate descent in a path-wise fashion. Details may be found
in Friedman, Hastie, and Tibshirani ([2010](#ref-glmnet)), Simon et al.
([2011](#ref-coxnet)), Tibshirani et al. ([2012](#ref-strongrules)),
Simon, Friedman, and Hastie ([2013](#ref-block)), Tay, Narasimhan and
Hastie ([2023](#ref-glm)).

There are several vignettes under `Articles` which are helpful for
getting started with glmnet, and getting into some of the details. Since
the glmnet project is more than 18 years old, there is a document that
outlines the history.

The companion book [*glmnet is (almost) all you
need*](https://glmnet.stanford.edu/book/), by Tay, Craig, Narasimhan,
Tibshirani, and Hastie, brings the package’s many capabilities together
in one place, with worked examples. Accessible directly using the `Book`
link in the menu bar. This book is somewhat dynamic, since it is in
advanced draft mode, so will change slightly from time to time.

Version 5.0 is a major release with several new features, including:

  - All kernels are written in C++, so Fortran is now gone from glmnet;
  - All aspects of the Cox model are implemented efficiently, and the
    “Efron” method for tie-breaking is implemented;
  - `glmnet.control()` has been streamlined and reorganized, and absorbs
    some of the many arguments in earlier versions of glmnet.

These build on earlier versions, that included features such as:

  - Relaxed fitting to allow models in the path to be refit without
    regularization;
  - Assessment functions for displaying performance of models on test
    data. These include all the measures available via `cv.glmnet`, as
    well as confusion matrices and ROC plots for classification models;
  - Functions for building the `x` input matrix for `glmnet` that allow
    for *one-hot-encoding* of factor variables, appropriate treatment of
    missing values, and an option to create a sparse matrix if
    appropriate.  
  - A function for fitting unpenalized a single version of any of the
    GLMs of `glmnet`.
  - The family argument allows for any GLM family used by `glm()`, in
    addition to the built-in glmnet families.
  - The scope for survival modeling has been enlarged, allowing for
    (start, stop) data and strata. We provide a much-requested method
    for `survival:survfit`.
  - Any model can use an `x` in sparse matrix format.

## References

<div id="refs" class="references">

<div id="ref-glmnet">

Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2010.
“Regularization Paths for Generalized Linear Models via Coordinate
Descent.” *Journal of Statistical Software, Articles* 33 (1): 1–22.
<https://doi.org/10.18637/jss.v033.i01>.

</div>

<div id="ref-block">

Simon, Noah, Jerome Friedman, and Trevor Hastie. 2013. “A Blockwise
Descent Algorithm for Group-Penalized Multiresponse and Multinomial
Regression.”

</div>

<div id="ref-coxnet">

Simon, Noah, Jerome Friedman, Trevor Hastie, and Robert Tibshirani.
2011. “Regularization Paths for Cox’s Proportional Hazards Model via
Coordinate Descent.” *Journal of Statistical Software, Articles* 39 (5):
1–13. <https://doi.org/10.18637/jss.v039.i05>.

</div>

<div id="ref-strongrules">

Tibshirani, Robert, Jacob Bien, Jerome Friedman, Trevor Hastie, Noah
Simon, Jonathan Taylor, and Ryan Tibshirani. 2012. “Strong Rules for
Discarding Predictors in Lasso-Type Problems.” *Journal of the Royal
Statistical Society: Series B (Statistical Methodology)* 74 (2): 245–66.
<https://doi.org/10.1111/j.1467-9868.2011.01004.x>.

</div>

<div id="ref-glm">

Kenneth Tay, J, Narasimhan, Balasubramanian, Hastie, Trevor. 2023.
“Elastic Net Regularization Paths for All Generalized Linear Models.”
*Journal of Statistical Software, Articles* 106 (1): 1–31.
<https://doi.org/10.18637/jss.v106.i01>.

</div>

</div>
