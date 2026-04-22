#' Diagnostic Tool by Multiple Imputation for Regression Discontinuity Designs
#'
#' @description
#' Estimates average treatment effects at the cutoff based on sharp regression discontinuity
#' designs (RDD) and multiple imputation regression discontinuity designs (MIRDD).
#' It provides diagnostic tools for RDD by comparing results with those from MIRDD.
#'
#' @param y A numeric vector of the outcome variable.
#' @param x A numeric vector of the running variable (forcing variable).
#' @param cut A numeric value indicating the cutoff point in x. The user must supply a specific number.
#' @param seed A seed number for reproducibility. Default is NULL.
#' @param M1 Number of imputations for MIRDD. Default is 100.
#' @param M2 Number of imputations for visualization (plots 3, 4, 9, and 10). Default is 5. These datasets are the subsets of M1 imputed datasets. Thus, M2 cannot be larger than M1.
#' @param M3 Number of imputed datasets for plots 5 to 10. Default is 1. These datasets are the subsets of M1 imputed datasets. Thus, M3 cannot be larger than M1.
#' @param p2s1 Integer for Amelia's p2s argument (0 or 1), where 0 for no screen printing and 1 for screen printing of multiple imputation process. Default is 1.
#' @param emp Amelia's empirical (ridge) prior. Default is 0. A reasonable upper bound is 0.1.
#' @param bw Bandwidth selection method for rdrobust. Options are "mserd" (default), "msesum", "cerrd", and "cersum". "mserd" is one common MSE-optimal bandwidth selector. "msesum" is one common MSE-optimal bandwidth selector for the sum of regression estimates. "cerrd" is one common CER-optimal bandwidth selector. "cersum" is one common CER-optimal bandwidth selector for the sum of regression estimates. MSE is Mean Squared Error. CER is Coverage Error Rate.
#' @param ker Kernel function for rdrobust. Options are "triangular" (default option), "epanechnikov", and "uniform".
#' @param h Number for bandwidth. Default is NULL (data-driven).
#' @param type Inference type: "Conventional" (default), "Bias-Corrected", or "Robust".
#' @param p1 Polynomial order (1 or 2) for rdrobust and MIRDD. Default is 1 (local linear regression). Can take either 1 (local linear regression) or 2 (local quadratic regression). When specified larger than 2, it will be considered 2.
#' @param conf Confidence level (0-100). Default is 95.
#' @param upper If 1 (default), treatment is x >= cut. If 0, treatment is x < cut.
#' @param covs1 Optional covariates. If two additional covariates z1 and z2 need to be used, then covs1 = data.frame(z1, z2).
#' @param up Optional upper bound for imputed values.
#' @param lo Optional lower bound for imputed values.
#'
#' @return
#' \item{Estimate}{Estimated quantities of the average treatment effects (ATE) at the cutoff.}
#' \item{Std.Error}{Standard error of the estimate.}
#' \item{CI.LL}{Lower limit of the 95\% confidence interval.}
#' \item{CI.UL}{Upper limit of the 95\% confidence interval.}
#' \item{size}{Sub-sample size to estimate the ATE at the cutoff.}
#' \item{bandwidth}{Length of the bandwidth used for RDD analysis.}
#'
#' In addition to the data frame, a series of diagnostic plots are generated:
#' \describe{
#'   \item{1. MIRDD, RDD, Naive}{A diagnostic plot to visualize the relationship among the three estimators. Red vertical line is RDD, black solid line is naive, and histogram is MI.}
#'   \item{2. MIRDD and RDD}{A diagnostic plot to visualize the relationship between the two estimators. Red vertical line is RDD and histogram is MI.}
#'   \item{3. Densities (Control)}{A diagnostic plot to visualize the densities of observed and imputed data. Gray solid curve is the density of observed data in the control group. Blue solid curve is the density of observed data in the treatment group. Red dashed lines are the densities of imputed data in the control group.}
#'   \item{4. Densities (Treatment)}{A diagnostic plot to visualize the densities of observed and imputed data. Gray solid curve is the density of observed data in the control group. Blue solid curve is the density of observed data in the treatment group. Red dashed lines are the densities of imputed data in the treatment group.}
#'   \item{5. Observed Values}{A diagnostic plot to visualize the scatterplot of observed data. Gray circles are observed data in the control group. Blue triangles are observed data in the treatment group.}
#'   \item{6. Observed & Imputed Values}{A diagnostic plot to visualize the scatterplot of observed and imputed data. Red circles are imputed data in the control group. Red triangles are imputed data in the treatment group. These imputed data are overlaid on the observed data in Figure 5.}
#'   \item{7. Observed & Imputed (Control)}{A diagnostic plot to clearly visualize the scatterplot of observed and imputed data in the control group only.}
#'   \item{8. Observed & Imputed (Treatment)}{A diagnostic plot to clearly visualize the scatterplot of observed and imputed data in the treatment group only.}
#'   \item{9. Around Cutoff (Control)}{A diagnostic plot to clearly visualize the scatterplot, around the cutoff point, of observed and imputed data in the control group only. Five solid lines are the estimated linear regression lines based on multiply imputed data.}
#'   \item{10. Around Cutoff (Treatment)}{A diagnostic plot to clearly visualize the scatterplot, around the cutoff point, of observed and imputed data in the treatment group only. Five solid lines are the estimated linear regression lines based on multiply imputed data.}
#'   \item{11. Local Slope (Control)}{A diagnostic plot to visualize the distribution of the coefficients of the estimated linear regression models around the cutoff point in the control group.}
#'   \item{12. Local Slope (Treatment)}{A diagnostic plot to visualize the distribution of the coefficients of the estimated linear regression models around the cutoff point in the treatment group.}
#' }
#'
#' @references
#' Takahashi, M. 2023. Multiple imputation regression discontinuity designs: Alternative to regression discontinuity designs to estimate the local average treatment effect at the cutoff. Communications in Statistics - Simulation and Computation 53(9): 4293-4312. \doi{10.1080/03610918.2021.1960374}
#'
#' Calonico, S., Cattaneo, M.D., and Titiunik, R. 2015. rdrobust: An R Package for robust nonparametric inference in regression-discontinuity designs. R Journal 7(1): 38-51. \doi{10.32614/RJ-2015-004}
#'
#' Honaker, J., King, G., and Blackwell, M. 2011. Amelia II: A program for missing data. Journal of Statistical Software 45(7): 1-47. \doi{10.18637/jss.v045.i07}
#'
#' @import Amelia
#' @import rdrobust
#' @importFrom stats na.omit sd var qt lm density coef
#' @importFrom graphics layout hist abline plot par lines points
#' @export
#'
#' @examples
#' # Example usage with dummy data
#' x <- runif(100, -1, 1)
#' y <- 0.5 * x + (x >= 0) + rnorm(100, 0, 0.1)
#' MIdiagRDD(y = y, x = x, cut = 0)
MIdiagRDD <- function(y, x, cut, seed = NULL, M1 = 100, M2 = 5, M3 = 1, p2s1 = 1, emp = 0,
                      bw = "mserd", ker = "triangular", h = NULL, type = "Conventional",
                      p1 = 1, conf = 95, upper = 1, covs1 = NULL, up = NULL, lo = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Saving and Restoring Drawing Settings
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  # --- 1. Data Formatting ---
  df1 <- stats::na.omit(data.frame(x = x, y = y))
  p1 <- min(p1, 2)
  x1_all <- df1$x
  y3_all <- df1$y
  d_all  <- if (upper == 1) as.numeric(x1_all >= cut) else as.numeric(x1_all < cut)

  y0_all <- ifelse(d_all == 0, y3_all, NA)
  y1_all <- ifelse(d_all == 1, y3_all, NA)

  data3_all <- data.frame(y3 = y3_all, x1 = x1_all, d = d_all)
  if (!is.null(covs1)) data3_all <- cbind(data3_all, covs1)

  # --- 2. Bandwidth ---
  rbd1 <- rdrobust::rdrobust(y3_all, x1_all, c = cut, kernel = ker, bwselect = bw, p = p1, covs = covs1)

  if (!is.null(h)) {
    if (length(h) == 1) {
      h1 <- h
      h2 <- h
    } else {
      h1 <- h[1]
      h2 <- h[2]
    }
  } else {
    if (type == "Conventional") {
      h1 <- rbd1$bws[1, 1]
      h2 <- rbd1$bws[1, 2]
    } else {
      # Bias-Corrected or Robust
      h1 <- rbd1$bws[2, 1]
      h2 <- rbd1$bws[2, 2]
    }
  }

  # Extracting data within bandwidth
  in_bw  <- (cut - h1 <= x1_all) & (x1_all <= cut + h2)
  data2a <- data.frame(y0 = y0_all[in_bw], y1 = y1_all[in_bw], x1 = x1_all[in_bw])
  if (!is.null(covs1)) {
    covs_sub <- if(is.vector(covs1)) covs1[in_bw] else covs1[in_bw, , drop=FALSE]
    data2a   <- cbind(data2a, covs_sub)
  }

  # If Bias-Corrected or Robust, then use x1^2
  if (type != "Conventional") {
    data2a$x2 <- (data2a$x1)^2
  }

  MIn <- nrow(data2a)
  da  <- if (upper == 1) as.numeric(data2a$x1 >= cut) else as.numeric(data2a$x1 < cut)

  # --- 3. Naive Estimator ---
  y1n_v <- stats::na.omit(data2a$y1); y0n_v <- stats::na.omit(data2a$y0)
  n1n <- length(y1n_v); n0n <- length(y0n_v)
  se1n <- stats::sd(y1n_v) / sqrt(n1n); se0n <- stats::sd(y0n_v) / sqrt(n0n)
  dfna <- (se1n^2 + se0n^2)^2 / (se1n^4/(n1n-1) + se0n^4/(n0n-1))
  naive1 <- mean(y1n_v) - mean(y0n_v)
  naive2 <- se1n + se0n

  # --- 4. MI: Amelia ---
  a.out <- Amelia::amelia(data2a, p2s = p2s1, m = M1, empri = emp * MIn)
  if (is.null(up)) up <- max(y3_all, na.rm = TRUE)
  if (is.null(lo)) lo <- min(y3_all, na.rm = TRUE)

  tauhata <- matrix(NA, MIn, M1); tauhat0 <- matrix(NA, MIn, M1); tauhat1 <- matrix(NA, MIn, M1)
  nMIa1 <- numeric(M1); nMIa2 <- numeric(M1)

  for (i in 1:M1) {
    y1i <- a.out$imputations[[i]]$y1; y0i <- a.out$imputations[[i]]$y0
    y1i[y1i > up] <- up; y1i[y1i < lo] <- lo
    y0i[y0i > up] <- up; y0i[y0i < lo] <- lo
    tauhata[, i] <- y1i - y0i; tauhat0[, i] <- y0i; tauhat1[, i] <- y1i
    nMIa1[i] <- mean(tauhata[, i]); nMIa2[i] <- stats::var(tauhata[, i]) / MIn
  }
  mia1 <- mean(nMIa1); mia2 <- sqrt(mean(nMIa2) + (1 + 1/M1) * stats::var(nMIa1))

  # --- 5. RDD ---
  if (type == "Bias-Corrected") {
    modelRDD1 <- rdrobust::rdrobust(data3_all$y3, data3_all$x1, c = cut, all = TRUE, h = c(h1, h2), kernel = ker, p = p1, covs = covs1)
    rddb1 <- modelRDD1$Estimate[2]
    rddb2 <- modelRDD1$se[1]
    rdd_n_eff <- sum(modelRDD1$N_b)
  } else if (type == "Robust") {
    modelRDD1 <- rdrobust::rdrobust(data3_all$y3, data3_all$x1, c = cut, all = TRUE, h = c(h1, h2), kernel = ker, p = p1, covs = covs1)
    rddb1 <- modelRDD1$Estimate[2]
    rddb2 <- modelRDD1$se[2]
    rdd_n_eff <- sum(modelRDD1$N_b)
  } else {
    # Conventional (Default)
    modelRDD1 <- rdrobust::rdrobust(data3_all$y3, data3_all$x1, c = cut, all = TRUE, h = c(h1, h2), kernel = ker, p = p1, covs = covs1)
    rddb1 <- modelRDD1$Estimate[1]
    rddb2 <- modelRDD1$se[1]
    rdd_n_eff <- sum(modelRDD1$N_h)
  }

  # --- 6. Graph Plotting ---
  graphics::layout(matrix(1:12, 3, 4, byrow = TRUE))

  graphics::hist(nMIa1, xlim = range(c(nMIa1, rddb1, naive1)), main = "1.MIRDD, RDD, Naive", xlab = "ATE at cutoff")
  graphics::abline(v = c(rddb1, naive1), col = c(2, 1), lwd = 2)
  graphics::hist(nMIa1, xlim = range(c(nMIa1, rddb1)), main = "2.MIRDD and RDD", xlab = "ATE at cutoff")
  graphics::abline(v = rddb1, col = 2, lwd = 2)

  y0n_o <- data2a$y0[da == 0]; y1n_o <- data2a$y1[da == 1]
  rng_x_den <- range(c(tauhat0, tauhat1, y0n_o, y1n_o))
  graphics::plot(stats::density(y0n_o), xlim = rng_x_den, lwd = 2, col = 8, main = "3.Densities (Ctrl)")
  for (i in 1:M2) graphics::lines(stats::density(tauhat0[, i][a.out$missMatrix[, 1]]), col = 2, lty = 2)
  graphics::plot(stats::density(y1n_o), xlim = rng_x_den, lwd = 2, col = 4, main = "4.Densities (Trt)")
  for (i in 1:M2) graphics::lines(stats::density(tauhat1[, i][a.out$missMatrix[, 2]]), col = 2, lty = 2)

  miss0 <- a.out$missMatrix[, 1]; miss1 <- a.out$missMatrix[, 2]
  rng_ys <- range(y3_all); rng_xs <- range(x1_all)

  graphics::plot(x1_all[d_all==0], y0_all[d_all==0], xlim = rng_xs, ylim = rng_ys, col = 8, pch = 1, main = "5.Observed Values", xlab="x", ylab="y")
  graphics::points(x1_all[d_all==1], y1_all[d_all==1], col = 4, pch = 2)
  graphics::abline(v = cut, lwd = 1)

  graphics::plot(x1_all[d_all==0], y0_all[d_all==0], xlim = rng_xs, ylim = rng_ys, col = 8, pch = 1, main = "6.Observed & Imputed", xlab="x", ylab="y")
  graphics::points(x1_all[d_all==1], y1_all[d_all==1], col = 4, pch = 2)
  for (i in 1:M3) {
    graphics::points(data2a$x1[miss1], tauhat1[miss1, i], col = 2, pch = 2, cex = 1.0)
    graphics::points(data2a$x1[miss0], tauhat0[miss0, i], col = 2, pch = 1, cex = 1.0)
  }
  graphics::abline(v = cut, lwd = 1)

  graphics::plot(x1_all[d_all==0], y0_all[d_all==0], xlim = rng_xs, ylim = rng_ys, col = 8, pch = 1, main = "7.Obs & Imp (Ctrl)", xlab="x", ylab="y")
  for (i in 1:M3) graphics::points(data2a$x1[miss0], tauhat0[miss0, i], col = 2, pch = 1, cex = 1.0)
  graphics::abline(v = cut, lwd = 1)

  graphics::plot(x1_all[d_all==1], y1_all[d_all==1], xlim = rng_xs, ylim = rng_ys, col = 4, pch = 2, main = "8.Obs & Imp (Trt)", xlab="x", ylab="y")
  for (i in 1:M3) graphics::points(data2a$x1[miss1], tauhat1[miss1, i], col = 2, pch = 2, cex = 1.0)
  graphics::abline(v = cut, lwd = 1)

  graphics::plot(data2a$x1[da==0], data2a$y0[da==0], xlim = c(cut-h1, cut+h2), ylim = rng_ys, col = 8, pch = 1, main = "9.Local Reg (Ctrl)", xlab="x", ylab="y")
  for (i in 1:M2) {
    graphics::points(data2a$x1[miss0], tauhat0[miss0, i], col = 2, pch = 1, cex = 1.0)
    graphics::abline(stats::lm(tauhat0[, i] ~ data2a$x1), lwd = 1.0)
  }

  graphics::plot(data2a$x1[da==1], data2a$y1[da==1], xlim = c(cut-h1, cut+h2), ylim = rng_ys, col = 4, pch = 2, main = "10.Local Reg (Trt)", xlab="x", ylab="y")
  for (i in 1:M2) {
    graphics::points(data2a$x1[miss1], tauhat1[miss1, i], col = 2, pch = 2, cex = 1.0)
    graphics::abline(stats::lm(tauhat1[, i] ~ data2a$x1), lwd = 1.0)
  }

  b0 <- numeric(M1); b1 <- numeric(M1)
  for (i in 1:M1) {
    b0[i] <- stats::coef(stats::lm(tauhat0[, i] ~ data2a$x1))[2]
    b1[i] <- stats::coef(stats::lm(tauhat1[, i] ~ data2a$x1))[2]
  }
  graphics::hist(b0, main = "11.Slope (Ctrl)"); graphics::hist(b1, main = "12.Slope (Trt)")

  # --- 7. Output ---
  cv_func <- function(df_val) stats::qt((1 - (1 - conf/100)/2), df_val)
  est <- c(mia1, rddb1, naive1)
  se  <- c(mia2, rddb2, naive2)
  dff <- c(MIn - 1, rdd_n_eff - 1, dfna)

  output <- data.frame(
    Method = c("MIRDD", "RDD", "Naive"),
    Estimate = round(est, 4),
    Std.Error = round(se, 4),
    CI.LL = round(est - cv_func(dff) * se, 4),
    CI.UL = round(est + cv_func(dff) * se, 4),
    Size = c(MIn, rdd_n_eff, n1n + n0n),
    Bandwidth = round(h1, 4)
  )
  return(output)
}
