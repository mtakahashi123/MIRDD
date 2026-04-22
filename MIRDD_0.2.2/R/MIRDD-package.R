#' MIRDD: Diagnostic Tool by Multiple Imputation for Regression Discontinuity Designs
#'
#' @description
#' This package implements the method proposed in Takahashi (2023), which provides
#' a novel framework for the regression discontinuity design (RDD) by reinterpreting
#' the estimation of treatment effects as a missing data problem. While standard RDD
#' relies on observations near the cutoff, MIRDD utilizes multiple imputation to
#' account for missing potential outcomes, offering a diagnostic tool to assess
#' the validity and robustness of RDD estimates.
#'
#' @details
#' The main function in this package is \code{\link{MIdiagRDD}}.
#'
#' @references
#' Takahashi, M. 2023. Multiple imputation regression discontinuity designs:
#' Alternative to regression discontinuity designs to estimate the local average
#' treatment effect at the cutoff. Communications in Statistics - Simulation and
#' Computation 52 (9): 4293-4312. \doi{10.1080/03610918.2021.1960374}
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \doi{10.1080/03610918.2021.1960374}
#' }
#'
"_PACKAGE"
