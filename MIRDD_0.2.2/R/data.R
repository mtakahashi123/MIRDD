#' Lee (2008) Dataset
#'
#' @description A dataset used in Takahashi (2023).
#'
#' @usage data(lee2008)
#'
#' @format A data frame with 6558 observations on the following 3 variables.
#' \describe{
#'   \item{y1}{a numeric vector (dependent variable). Democrat vote share election at t + 1.}
#'   \item{x1}{a numeric vector (running variable). The cutoff point is where x1 = 0. Democratic vote share at t.}
#'   \item{x2}{a numeric vector (additional covariate). Democratic vote share at t - 1}
#' }
#'
#' @source Lee, D. S. 2008. Randomized experiments from non-random selection in U.S. House elections. Journal of Econometrics 142 (2): 675-697.
#'
#' @references Takahashi, M. 2023. Multiple imputation regression discontinuity designs: Alternative to regression discontinuity designs to estimate the local average treatment effect at the cutoff. Communications in Statistics - Simulation and Computation 52 (9): 4293-4312.
"lee2008"

#' Ludwig and Miller (2007) Modified Dataset
#'
#' @description A modified dataset used in Takahashi (2023).
#'
#' @usage data(LudwigMiller2007Modified)
#'
#' @format A data frame with 1037 observations on the following 11 variables.
#' \describe{
#'   \item{y1}{a numeric vector (dependent variable)}
#'   \item{x1}{a numeric vector (running variable). The cutoff point is where x1 = log(59.1984).}
#'   \item{z1}{a numeric vector (additional covariate)}
#'   \item{z2}{a numeric vector(additional covariate)}
#'   \item{z3}{a numeric vector(additional covariate)}
#'   \item{z4}{a numeric vector(additional covariate)}
#'   \item{z5}{a numeric vector(additional covariate)}
#'   \item{z6}{a numeric vector(additional covariate)}
#'   \item{z7}{a numeric vector(additional covariate)}
#'   \item{z8}{a numeric vector(additional covariate)}
#'   \item{z9}{a numeric vector(additional covariate)}
#' }
#'
#' @source Ludwig, J., and D. L. Miller. 2007. Does Head Start improve children's life chances? Evidence from a regression discontinuity design. Quarterly Journal of Economics 122 (1): 159-208.
#'
#' @references Takahashi, M. 2023. Multiple imputation regression discontinuity designs: Alternative to regression discontinuity designs to estimate the local average treatment effect at the cutoff. Communications in Statistics - Simulation and Computation 52 (9): 4293-4312.
"LudwigMiller2007Modified"