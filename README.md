
<!-- README.md is generated from the source: README.Rmd -->

# Lasso and Elastic-Net Regularized Generalized Linear Models <img src="man/figures/logo.png" width="100" align="right" />

<!-- [![Travis-CI Build -->

<!-- Status](https://travis-ci.org/trevorhastie/glmnet.svg?branch=master)](https://travis-ci.org/trevorhastie/glmnet) -->

<!-- [![Coverage -->

<!-- Status](https://img.shields.io/codecov/c/github/trevorhastie/glmnet/master.svg)](https://codecov.io/github/trevorhastie/glmnet?branch=master) -->

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/glmnet)](https://cran.r-project.org/package=glmnet)[![](https://cranlogs.r-pkg.org/badges/glmnet)](https://CRAN.R-project.org/package=glmnet)

We provide extremely efficient procedures for fitting the entire lasso
or elastic-net regularization path for linear regression (gaussian),
multi-task gaussian, logistic and multinomial regression models (grouped
or not), Poisson regression and the Cox model. The algorithm uses
cyclical coordinate descent in a path-wise fashion. Details may be found
in Friedman, Hastie, and Tibshirani ([2010](#ref-glmnet)), Simon et al.
([2011](#ref-coxnet)), Tibshirani et al. ([2012](#ref-strongrules)),
Simon, Friedman, and Hastie ([2013](#ref-block)).

Version 3.0 is a major release with several new features, including:

  - Relaxed fitting to allow models in the path to be refit without
    regularization. CV will select from these, or from specified
    mixtures of the relaxed fit and the regular fit;
  - Progress bar to monitor computation;
  - Assessment functions for displaying performance of models on test
    data. These include all the measures available via `cv.glmnet`, as
    well as confusion matrices and ROC plots for classification models;
  - print methods for CV output;
  - Functions for building the `x` input matrix for `glmnet` that allow
    for *one-hot-encoding* of factor variables, appropriate treatment of
    missing values, and an option to create a sparse matrix if
    appropriate.
  - A function for fitting unpenalized a single version of any of the
    GLMs of `glmnet`.

Version 4.0 is a major release that allows for any GLM family, besides
the built-in families.

Version 4.3 is a major release that expands the scope for survival
modeling, allowing for (start, stop) data, strata, and sparse X inputs.
It also provides a much-requested method for `survival:survfit`.

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

</div>
