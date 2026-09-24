setwd("C:/Users/DAIRX/OneDrive - AbbVie Inc (O365)/Documents/tehTMLE/Yixin_Revision/Rcodes")

## ==================================================================
## Simulation driver: ORIGINAL (non-moderate-overlap) DGP, tilted h(x)
## only, using tmle_h_methods.R for all estimation (MLE/IPW/AIPW/TMLE,
## EIF-Wald + target-step-bootstrap inference). Propensity-score
## truncation ON ([0.025, 0.975]).
##
## DGP matches the simulation-section writeup: n1 = 200 treated,
## n0 = 300 external controls;
##   X1 ~ Bernoulli(0.68) / Bernoulli(0.32)      [treated / control]
##   X2 ~ N(0.75, 1.00^2) / N(-0.55, 0.85^2)
##   X3 ~ N(0.45, 0.80^2) / N(-0.35, 1.15^2)
##   X4 ~ Uniform(-1, 1)  [same in both arms -> drops out of the true g]
## Outcome model (matches FixedX_Rcodes_8.R, the script that produced
## the original table3_tilted.csv/table4_tilted.csv/ps_overlap.pdf):
##   logit Q0(X,0) = -1.20 + 0.60 X1 + 2.10 sin(1.40 X2) + 1.10 tanh(X3)
##                   + 0.85 (X3^2 - 1) + 0.60 X2 X4
##   logit Q0(X,1) = logit Q0(X,0) + 0.60 + 0.75 X1 + 1.75 sin(1.70 X2)
##                   - 1.30 (X3^2 - 1) + 1.00 X2 X4
## True propensity score (X4 drops out of the likelihood ratio, same
## covariate marginals as above):
##   g(X) = expit(-1.14 + 1.51 X1 + 1.51 X2 + 0.19 X2^2 + 0.97 X3 - 0.40 X3^2)
## Tilted h(x) = exp(-0.18 X2^2 - 0.10 X3^2), unchanged from prior scripts.
## ==================================================================

rm(list = ls())

source("tmle_h_methods.R")


## ------------------------------------------------------------------
## SECTION 1: WORKING MODELS (Sec 5.2.3) -- functional forms only,
## independent of the DGP's specific means/variances.
## ------------------------------------------------------------------

Qform_wrong   <- Y ~ A + X1 + X2 + X3 + X4
Qform_correct <- Y ~ A + X1 + sin(1.4*X2) + tanh(X3) + I(X3^2) + X2:X4 +
                     A:X1 + A:sin(1.7*X2) + A:I(X3^2) + A:X2:X4
gform_wrong   <- A ~ X1 + X2 + X3 + X4
gform_correct <- A ~ X1 + X2 + X3 + I(X2^2) + I(X3^2)


## ------------------------------------------------------------------
## SECTION 2: DATA GENERATION -- original-overlap DGP
## ------------------------------------------------------------------

simulate_dgp_original <- function(n1 = 200, n0 = 300) {
  treated <- data.frame(A = 1L, X1 = rbinom(n1, 1, 0.68), X2 = rnorm(n1, 0.75, 1.00),
                         X3 = rnorm(n1, 0.45, 0.80), X4 = runif(n1, -1, 1))
  external <- data.frame(A = 0L, X1 = rbinom(n0, 1, 0.32), X2 = rnorm(n0, -0.55, 0.85),
                          X3 = rnorm(n0, -0.35, 1.15), X4 = runif(n0, -1, 1))
  dat <- rbind(treated, external)

  eta0 <- with(dat, -1.20 + 0.60*X1 + 2.10*sin(1.40*X2) + 1.10*tanh(X3) +
                 0.85*(X3^2 - 1) + 0.60*X2*X4)
  contrast <- with(dat, 0.60 + 0.75*X1 + 1.75*sin(1.70*X2) -
                      1.30*(X3^2 - 1) + 1.00*X2*X4)

  Q0_true <- plogis(eta0)
  Q1_true <- plogis(eta0 + contrast)
  Y <- rbinom(nrow(dat), 1, ifelse(dat$A == 1L, Q1_true, Q0_true))
  h <- exp(-0.18 * dat$X2^2 - 0.10 * dat$X3^2)

  list(Y = Y, A = dat$A, W = dat[, c("X1", "X2", "X3", "X4")], h = h,
       Q0_true = Q0_true, Q1_true = Q1_true)
}

## True g(x) via Bayes' rule (X4 is shared by both arms and drops out
## of the likelihood ratio). Reference/plot use only -- NOT used to fit
## anything (g is always estimated from data via gform_wrong/correct).
true_g_original <- function(X1, X2, X3, n1 = 200, n0 = 300) {
  f1 <- dbinom(X1, 1, 0.68) * dnorm(X2, 0.75, 1.00) * dnorm(X3, 0.45, 0.80)
  f0 <- dbinom(X1, 1, 0.32) * dnorm(X2, -0.55, 0.85) * dnorm(X3, -0.35, 1.15)
  (n1 * f1) / (n1 * f1 + n0 * f0)
}


## ------------------------------------------------------------------
## SECTION 3: MONTE CARLO TRUTH FOR THE TILTED h
## Computed once (h does not depend on g), via chunked large-sample
## draws from the same design mixture (n1 : n0 = 2:3).
## ------------------------------------------------------------------

approximate_true_tau_tilted <- function(n1_large = 2e6, n0_large = 3e6,
                                         chunk_size = 250000, seed = 20260625) {
  set.seed(seed)
  numerator <- 0; denominator <- 0

  draw_chunk <- function(n, arm) {
    if (arm == 1L) {
      X1 <- rbinom(n, 1, 0.68); X2 <- rnorm(n, 0.75, 1.00); X3 <- rnorm(n, 0.45, 0.80)
    } else {
      X1 <- rbinom(n, 1, 0.32); X2 <- rnorm(n, -0.55, 0.85); X3 <- rnorm(n, -0.35, 1.15)
    }
    X4 <- runif(n, -1, 1)
    eta0 <- -1.20 + 0.60*X1 + 2.10*sin(1.40*X2) + 1.10*tanh(X3) +
      0.85*(X3^2 - 1) + 0.60*X2*X4
    contrast <- 0.60 + 0.75*X1 + 1.75*sin(1.70*X2) -
      1.30*(X3^2 - 1) + 1.00*X2*X4
    Q0 <- plogis(eta0); Q1 <- plogis(eta0 + contrast)
    h <- exp(-0.18*X2^2 - 0.10*X3^2)
    c(num = sum(h * (Q1 - Q0)), den = sum(h))
  }

  for (arm in c(1L, 0L)) {
    remaining <- if (arm == 1L) n1_large else n0_large
    while (remaining > 0) {
      n_now <- min(chunk_size, remaining)
      z <- draw_chunk(n_now, arm)
      numerator <- numerator + z["num"]; denominator <- denominator + z["den"]
      remaining <- remaining - n_now
    }
  }
  as.numeric(numerator / denominator)
}


## ------------------------------------------------------------------
## SECTION 4: ONE REPLICATION (tilted h only, 4 nuisance scenarios)
## All estimation is delegated to tmle_h_methods.R's fit_Q()/fit_g()/
## estimate_tau_h_given_fits() -- no estimator logic is redefined here.
## ------------------------------------------------------------------

one_simulation_tilted <- function(true_tau, n1 = 200, n0 = 300, n_boot = 1000,
                                   gbounds = c(0.025, 0.975)) {
  d <- simulate_dgp_original(n1, n0)
  Y <- d$Y; A <- d$A; W <- d$W; h <- d$h

  Q_wrong   <- fit_Q(Y, A, W, method = "glm", Qform = Qform_wrong)
  Q_correct <- fit_Q(Y, A, W, method = "glm", Qform = Qform_correct)
  g_wrong   <- fit_g(A, W, method = "glm", gform = gform_wrong,   gbounds = gbounds)
  g_correct <- fit_g(A, W, method = "glm", gform = gform_correct, gbounds = gbounds)

  scenarios <- list(
    "Q misspecified; g correct" = list(Q = Q_wrong,   g = g_correct),
    "Q correct; g misspecified" = list(Q = Q_correct, g = g_wrong),
    "Both misspecified"         = list(Q = Q_wrong,   g = g_wrong),
    "Both correctly specified"  = list(Q = Q_correct, g = g_correct)
  )

  rows <- lapply(names(scenarios), function(sc_name) {
    Qfit <- scenarios[[sc_name]]$Q; gfit <- scenarios[[sc_name]]$g

    out <- estimate_tau_h_given_fits(Y, A, h, Qfit = Qfit, gfit = gfit, n_boot = n_boot)

    data.frame(
      Scenario = sc_name, H = "tilted", Truth = true_tau,
      MLE = out$MLE, IPW = out$IPW, AIPW = out$AIPW, TMLE = out$TMLE,
      TMLE_EIF_SE = out$TMLE_se,
      TMLE_EIF_Cover = as.integer(out$TMLE_CI[1] <= true_tau & true_tau <= out$TMLE_CI[2]),
      TMLE_Bootstrap_SE = out$Bootstrap_SE,
      TMLE_Bootstrap_Normal_Cover = as.integer(out$Bootstrap_CI_normal[1] <= true_tau & true_tau <= out$Bootstrap_CI_normal[2]),
      TMLE_Bootstrap_Percentile_Cover = as.integer(out$Bootstrap_CI_percentile[1] <= true_tau & true_tau <= out$Bootstrap_CI_percentile[2])
    )
  })

  do.call(rbind, rows)
}


## ------------------------------------------------------------------
## SECTION 5: DRIVER + SUMMARY TABLES
## ------------------------------------------------------------------

methods_used <- c("MLE", "IPW", "AIPW", "TMLE")

run_tilted_study <- function(B = 1000, n_boot = 1000, n1 = 200, n0 = 300,
                              gbounds = c(0.025, 0.975), seed_offset = 840000) {
  true_tau <- approximate_true_tau_tilted()

  seeds <- seed_offset + seq_len(B)
  do.call(rbind, lapply(seeds, function(s) {
    set.seed(s)
    one_simulation_tilted(true_tau, n1 = n1, n0 = n0, n_boot = n_boot, gbounds = gbounds)
  }))
}

summarize_table3 <- function(results) {
  do.call(rbind, lapply(unique(results$Scenario), function(sc) {
    d <- results[results$Scenario == sc, ]
    do.call(rbind, lapply(methods_used, function(m) {
      est <- d[[m]]
      dev <- est - d$Truth
      data.frame(Scenario = sc, Method = m,
                 Mean_Truth = mean(d$Truth), Mean = mean(est),
                 Abs_Bias = abs(mean(dev)), EmpSD = sd(est), RMSE = sqrt(mean(dev^2)))
    }))
  }))
}

summarize_table4 <- function(results) {
  do.call(rbind, lapply(unique(results$Scenario), function(sc) {
    d <- results[results$Scenario == sc, ]
    rbind(
      data.frame(Scenario = sc, Interval = "EIF-Wald",
                 Mean_SE = mean(d$TMLE_EIF_SE), Coverage = mean(d$TMLE_EIF_Cover)),
      data.frame(Scenario = sc, Interval = "Target-step bootstrap, normal",
                 Mean_SE = mean(d$TMLE_Bootstrap_SE), Coverage = mean(d$TMLE_Bootstrap_Normal_Cover)),
      data.frame(Scenario = sc, Interval = "Target-step bootstrap, percentile",
                 Mean_SE = NA_real_, Coverage = mean(d$TMLE_Bootstrap_Percentile_Cover))
    )
  }))
}


## ------------------------------------------------------------------
## SECTION 6: RUN + WRITE table3_tilted_boot1000.csv / table4_tilted_boot1000.csv
## (renamed from table3_tilted.csv/table4_tilted.csv so the originals,
## produced with n_boot = 300, are not overwritten)
## ------------------------------------------------------------------

gbounds_run <- c(0.025, 0.975)   # truncation ON

results_tilted <- run_tilted_study(B = 1000, n_boot = 300, gbounds = gbounds_run)
write.csv(summarize_table3(results_tilted), "table3_tilted.csv")
write.csv(summarize_table4(results_tilted), "table4_tilted.csv")


## ------------------------------------------------------------------
## SECTION 7: ps_overlap.pdf -- estimated propensity overlap (g correct)
## ------------------------------------------------------------------

plot_ps_overlap <- function(A, g1W, breaks = seq(0, 1, by = 0.05),
                             gbounds = c(0.025, 0.975),
                             main = "Estimated propensity score overlap (g correct)") {
  col_treated <- adjustcolor("#2a78d6", alpha.f = 0.55)  # categorical slot 1 (blue)
  col_control <- adjustcolor("#eb6834", alpha.f = 0.55)  # categorical slot 2 (orange)

  h1 <- hist(g1W[A == 1], breaks = breaks, plot = FALSE)
  h0 <- hist(g1W[A == 0], breaks = breaks, plot = FALSE)
  ymax <- max(h1$density, h0$density) * 1.05

  plot(h1, freq = FALSE, col = col_treated, border = "white",
       xlim = range(breaks), ylim = c(0, ymax),
       main = main, xlab = "Estimated P(A = 1 | X)", ylab = "Density")
  plot(h0, freq = FALSE, col = col_control, border = "white", add = TRUE)

  if (!is.null(gbounds)) {
    abline(v = gbounds, lty = 2, col = "#898781")   # propensity truncation bounds
  }

  legend("topright", legend = c("Treated (A = 1)", "Control (A = 0)"),
         fill = c(col_treated, col_control), border = NA, bty = "n")
}

d_plot <- simulate_dgp_original(n1 = 200, n0 = 300)
g_correct_plot <- fit_g(d_plot$A, d_plot$W, method = "glm", gform = gform_correct)

pdf("ps_overlap.pdf", width = 7, height = 5)
plot_ps_overlap(d_plot$A, g_correct_plot$g1W)
dev.off()
