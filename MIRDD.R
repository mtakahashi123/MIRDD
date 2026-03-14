#' MIRDD: Multiple Imputation Regression Discontinuity Design
#'
#' @import Amelia
#' @import rdrobust
#' @importFrom stats na.omit sd var qt lm density coef
#' @importFrom graphics layout hist abline plot par lines points
#' @export
MIdiagRDD <- function(y, x, cut, seed = 1, M1 = 100, M2 = 5, M3 = 1, 
                      p2s1 = 1, emp = 0, bw = "mserd", ker = "triangular", 
                      bwidth = 1, p1 = 1, conf = 95, upper = 1, 
                      covs1 = NULL, up = NULL, lo = NULL) {
  
  set.seed(seed)
  
  # Saving and Restoring Drawing Settings
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  # --- 1. Data Formatting ---
  df1 <- na.omit(data.frame(x = x, y = y))
  p1 <- min(p1, 2) 
  x1_all <- df1$x
  y3_all <- df1$y
  d_all  <- if (upper == 1) as.numeric(x1_all >= cut) else as.numeric(x1_all < cut)
  
  y0_all <- ifelse(d_all == 0, y3_all, NA)
  y1_all <- ifelse(d_all == 1, y3_all, NA)
  
  data3_all <- data.frame(y3 = y3_all, x1 = x1_all, d = d_all)
  if (!is.null(covs1)) data3_all <- cbind(data3_all, covs1)
  
  # --- 2. Selecting a Bandwidth ---
  rbd1 <- rdrobust::rdbwselect(data3_all$y3, data3_all$x1, c = cut, bwselect = bw, 
                               kernel = ker, p = p1, covs = covs1)
  h1 <- rbd1$bws[1] * bwidth
  h2 <- rbd1$bws[2] * bwidth
  
  # Extracting data within bandwidth
  in_bw  <- (cut - h1 <= x1_all) & (x1_all <= cut + h2)
  data2a <- data.frame(y0 = y0_all[in_bw], y1 = y1_all[in_bw], x1 = x1_all[in_bw])
  if (!is.null(covs1)) {
    covs_sub <- if(is.vector(covs1)) covs1[in_bw] else covs1[in_bw, , drop=FALSE]
    data2a   <- cbind(data2a, covs_sub)
  }
  if (p1 != 1) data2a$x2 <- data2a$x1^2
  
  MIn <- nrow(data2a)
  da  <- if (upper == 1) as.numeric(data2a$x1 >= cut) else as.numeric(data2a$x1 < cut)
  
  # --- 3. Naive Estimator ---
  y1n_v <- na.omit(data2a$y1); y0n_v <- na.omit(data2a$y0)
  n1n <- length(y1n_v); n0n <- length(y0n_v)
  se1n <- sd(y1n_v) / sqrt(n1n); se0n <- sd(y0n_v) / sqrt(n0n)
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
    nMIa1[i] <- mean(tauhata[, i]); nMIa2[i] <- var(tauhata[, i]) / MIn
  }
  mia1 <- mean(nMIa1); mia2 <- sqrt(mean(nMIa2) + (1 + 1/M1) * var(nMIa1))

  # --- 5. RDD ---
  modelRDD1 <- rdrobust::rdrobust(data3_all$y3, data3_all$x1, c = cut, all = TRUE, h = c(h1, h2), kernel = ker, p = p1, covs = covs1)
  rddb1 <- modelRDD1$Estimate[1] # Conventional Estimate
  rddb2 <- modelRDD1$se[1]       # Conventional SE
  rdd_n_eff <- sum(modelRDD1$N_h) # Effective sample size within the bandwidth

  # --- 6. Graph Plotting ---
  layout(matrix(1:12, 3, 4, byrow = TRUE))
  
  hist(nMIa1, xlim = range(c(nMIa1, rddb1, naive1)), main = "1.MI, RDD, Naive", xlab = "LATE")
  abline(v = c(rddb1, naive1), col = c(2, 1), lwd = 2)
  hist(nMIa1, xlim = range(c(nMIa1, rddb1)), main = "2.MI and RDD", xlab = "LATE")
  abline(v = rddb1, col = 2, lwd = 2)

  y0n_o <- data2a$y0[da == 0]; y1n_o <- data2a$y1[da == 1]
  rng_x_den <- range(c(tauhat0, tauhat1, y0n_o, y1n_o))
  plot(density(y0n_o), xlim = rng_x_den, lwd = 2, col = 8, main = "3.Densities (Ctrl)")
  for (i in 1:M2) lines(density(tauhat0[, i][a.out$missMatrix[, 1]]), col = 2, lty = 2)
  plot(density(y1n_o), xlim = rng_x_den, lwd = 2, col = 4, main = "4.Densities (Trt)")
  for (i in 1:M2) lines(density(tauhat1[, i][a.out$missMatrix[, 2]]), col = 2, lty = 2)

  miss0 <- a.out$missMatrix[, 1]; miss1 <- a.out$missMatrix[, 2]
  rng_ys <- range(y3_all); rng_xs <- range(x1_all)

  plot(x1_all[d_all==0], y0_all[d_all==0], xlim = rng_xs, ylim = rng_ys, col = 8, pch = 1, main = "5.Observed Values", xlab="x", ylab="y")
  points(x1_all[d_all==1], y1_all[d_all==1], col = 4, pch = 2)
  abline(v = cut, lwd = 1)

  plot(x1_all[d_all==0], y0_all[d_all==0], xlim = rng_xs, ylim = rng_ys, col = 8, pch = 1, main = "6.Observed & Imputed", xlab="x", ylab="y")
  points(x1_all[d_all==1], y1_all[d_all==1], col = 4, pch = 2)
  for (i in 1:M3) {
    points(data2a$x1[miss1], tauhat1[miss1, i], col = 2, pch = 2, cex = 0.5)
    points(data2a$x1[miss0], tauhat0[miss0, i], col = 2, pch = 1, cex = 0.5)
  }
  abline(v = cut, lwd = 1)

  plot(x1_all[d_all==0], y0_all[d_all==0], xlim = rng_xs, ylim = rng_ys, col = 8, pch = 1, main = "7.Obs & Imp (Ctrl)", xlab="x", ylab="y")
  for (i in 1:M3) points(data2a$x1[miss0], tauhat0[miss0, i], col = 2, pch = 1, cex = 0.5)
  abline(v = cut, lwd = 1)

  plot(x1_all[d_all==1], y1_all[d_all==1], xlim = rng_xs, ylim = rng_ys, col = 4, pch = 2, main = "8.Obs & Imp (Trt)", xlab="x", ylab="y")
  for (i in 1:M3) points(data2a$x1[miss1], tauhat1[miss1, i], col = 2, pch = 2, cex = 0.5)
  abline(v = cut, lwd = 1)

  plot(data2a$x1[da==0], data2a$y0[da==0], xlim = c(cut-h1, cut+h2), ylim = rng_ys, col = 8, pch = 1, main = "9.Local Reg (Ctrl)", xlab="x", ylab="y")
  for (i in 1:M2) {
    points(data2a$x1[miss0], tauhat0[miss0, i], col = 2, pch = 1, cex = 0.5)
    abline(lm(tauhat0[, i] ~ data2a$x1), lwd = 0.5)
  }

  plot(data2a$x1[da==1], data2a$y1[da==1], xlim = c(cut-h1, cut+h2), ylim = rng_ys, col = 4, pch = 2, main = "10.Local Reg (Trt)", xlab="x", ylab="y")
  for (i in 1:M2) {
    points(data2a$x1[miss1], tauhat1[miss1, i], col = 2, pch = 2, cex = 0.5)
    abline(lm(tauhat1[, i] ~ data2a$x1), lwd = 0.5)
  }

  b0 <- numeric(M1); b1 <- numeric(M1)
  for (i in 1:M1) {
    b0[i] <- coef(lm(tauhat0[, i] ~ data2a$x1))[2]
    b1[i] <- coef(lm(tauhat1[, i] ~ data2a$x1))[2]
  }
  hist(b0, main = "11.Slope (Ctrl)"); hist(b1, main = "12.Slope (Trt)")

  # --- 7. Output ---
  cv_func <- function(df_val) qt((1 - (1 - conf/100)/2), df_val)
  est <- c(mia1, rddb1, naive1)
  se  <- c(mia2, rddb2, naive2)
  dff <- c(MIn - 1, rdd_n_eff - 1, dfna)
  
  output <- data.frame(
    Method = c("MI", "RDD", "Naive"),
    Estimate = round(est, 4),
    Std.Error = round(se, 4),
    CI.LL = round(est - cv_func(dff) * se, 4),
    CI.UL = round(est + cv_func(dff) * se, 4),
    Size = c(MIn, rdd_n_eff, n1n + n0n),
    Bandwidth = round(h1, 4)
  )
  return(output)
}
