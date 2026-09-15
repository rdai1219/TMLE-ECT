
rm(list = ls())

## ==================================================================
## Simulation driver: moderate-overlap DGP (Sec 5.2.1), all 5 h-choices
## (tilted/ATE/ATT/ATC/ATO), propensity-score truncation ON ([0.025,0.975]).
## ==================================================================

results_dir <- file.path("..", "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

source("tmle_h_methods.R")


## ------------------------------------------------------------------
## SECTION 1: WORKING MODELS (Sec 5.2.3)
## ------------------------------------------------------------------

Qform_wrong   <- Y ~ A + X1 + X2 + X3 + X4
Qform_correct <- Y ~ A + X1 + sin(1.4*X2) + tanh(X3) + I(X3^2) + X2:X4 +
                     A:X1 + A:sin(1.7*X2) + A:I(X3^2) + A:X2:X4
gform_wrong   <- A ~ X1 + X2 + X3 + X4
gform_correct <- A ~ X1 + X2 + X3 + I(X2^2) + I(X3^2)


## ------------------------------------------------------------------
## SECTION 2: NAMED-ESTIMAND h(x) CHOICES 
## ------------------------------------------------------------------

true_g_5_2_1 <- function(X1, X2, X3, p1 = 0.4, p0 = 0.6) {
  f1 <- dbinom(X1, 1, 0.68) * dnorm(X2, 0.75, 1.00) * dnorm(X3, 0.45, 0.80)
  f0 <- dbinom(X1, 1, 0.45) * dnorm(X2, 0.00, 0.85) * dnorm(X3, -0.05, 1.15)
  (p1 * f1) / (p1 * f1 + p0 * f0)
}

h_fun_ATE <- function(X1, X2, X3, X4, g) rep(1, length(g))
h_fun_ATT <- function(X1, X2, X3, X4, g) g
h_fun_ATC <- function(X1, X2, X3, X4, g) 1 - g
h_fun_ATO <- function(X1, X2, X3, X4, g) g * (1 - g)


## ------------------------------------------------------------------
## SECTION 3: DATA GENERATION 
## ------------------------------------------------------------------

simulate_dgp_5_2_1 <- function(n1 = 200, n0 = 300) {
  treated <- data.frame(A = 1L, X1 = rbinom(n1, 1, 0.68), X2 = rnorm(n1, 0.75, 1.00),
                         X3 = rnorm(n1, 0.45, 0.80), X4 = runif(n1, -1, 1))
  external <- data.frame(A = 0L, X1 = rbinom(n0, 1, 0.45), X2 = rnorm(n0, 0.00, 0.85),
                          X3 = rnorm(n0, -0.05, 1.15), X4 = runif(n0, -1, 1))
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

approximate_true_tau_5_2_1 <- function(n1_large = 2e6, n0_large = 3e6,
                                        chunk_size = 250000, seed = 20260625) {
  set.seed(seed)
  numerator <- 0; denominator <- 0

  draw_chunk <- function(n, arm) {
    if (arm == 1L) {
      X1 <- rbinom(n, 1, 0.68); X2 <- rnorm(n, 0.75, 1.00); X3 <- rnorm(n, 0.45, 0.80)
    } else {
      X1 <- rbinom(n, 1, 0.45); X2 <- rnorm(n, 0.00, 0.85); X3 <- rnorm(n, -0.05, 1.15)
    }
    X4 <- runif(n, -1, 1)
    eta0 <- -1.20 + 0.60*X1 + 2.10*sin(1.40*X2) + 1.10*tanh(X3) +
      0.85*(X3^2 - 1) + 0.60*X2*X4
    contrast <- 0.60 + 0.75*X1 + 1.75*sin(1.70*X2) - 1.30*(X3^2 - 1) + 1.00*X2*X4
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

approximate_true_tau_general <- function(h_fun, n1_large = 2e6, n0_large = 3e6,
                                          chunk_size = 250000, seed = 20260625) {
  set.seed(seed)
  numerator <- 0; denominator <- 0

  draw_chunk <- function(n, arm) {
    if (arm == 1L) {
      X1 <- rbinom(n, 1, 0.68); X2 <- rnorm(n, 0.75, 1.00); X3 <- rnorm(n, 0.45, 0.80)
    } else {
      X1 <- rbinom(n, 1, 0.45); X2 <- rnorm(n, 0.00, 0.85); X3 <- rnorm(n, -0.05, 1.15)
    }
    X4 <- runif(n, -1, 1)
    eta0 <- -1.20 + 0.60*X1 + 2.10*sin(1.40*X2) + 1.10*tanh(X3) +
      0.85*(X3^2 - 1) + 0.60*X2*X4
    contrast <- 0.60 + 0.75*X1 + 1.75*sin(1.70*X2) - 1.30*(X3^2 - 1) + 1.00*X2*X4
    Q0 <- plogis(eta0); Q1 <- plogis(eta0 + contrast)
    g_true <- true_g_5_2_1(X1, X2, X3)
    h <- h_fun(X1, X2, X3, X4, g_true)
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

## Per-replication truth for h = eta(ghat) (ATT/ATC/ATO)
mc_truth_given_ghat <- function(g_fit, h_fun, mc_n = 10000,
                                 gbounds = c(0.025, 0.975),
                                 p1 = 0.4, p0 = 0.6) {
  n1_mc <- round(mc_n * p1); n0_mc <- round(mc_n * p0)

  draw_arm <- function(n, arm) {
    if (arm == 1L) {
      X1 <- rbinom(n, 1, 0.68); X2 <- rnorm(n, 0.75, 1.00); X3 <- rnorm(n, 0.45, 0.80)
    } else {
      X1 <- rbinom(n, 1, 0.45); X2 <- rnorm(n, 0.00, 0.85); X3 <- rnorm(n, -0.05, 1.15)
    }
    data.frame(X1 = X1, X2 = X2, X3 = X3, X4 = runif(n, -1, 1))
  }

  Wmc <- rbind(draw_arm(n1_mc, 1L), draw_arm(n0_mc, 0L))

  eta0 <- with(Wmc, -1.20 + 0.60*X1 + 2.10*sin(1.40*X2) + 1.10*tanh(X3) +
                 0.85*(X3^2 - 1) + 0.60*X2*X4)
  contrast <- with(Wmc, 0.60 + 0.75*X1 + 1.75*sin(1.70*X2) -
                      1.30*(X3^2 - 1) + 1.00*X2*X4)
  Q0 <- plogis(eta0); Q1 <- plogis(eta0 + contrast)

  ghat <- predict_ghat(g_fit, Wmc)          # from tmle_h_methods.R
  ghat <- pmin(pmax(ghat, gbounds[1]), gbounds[2])

  h <- h_fun(Wmc$X1, Wmc$X2, Wmc$X3, Wmc$X4, ghat)
  sum(h * (Q1 - Q0)) / sum(h)
}


## ------------------------------------------------------------------
## SECTION 4: ONE REPLICATION
## ------------------------------------------------------------------

one_simulation <- function(h_choice, true_tau = NULL, n1 = 200, n0 = 300, n_boot = 300,
                            gbounds = c(0.025, 0.975), mc_n = 10000) {
  h_choice <- match.arg(h_choice, c("tilted", "ATE", "ATT", "ATC", "ATO"))

  d <- simulate_dgp_5_2_1(n1, n0)
  Y <- d$Y; A <- d$A; W <- d$W

  Q_wrong   <- fit_Q(Y, A, W, method = "glm", Qform = Qform_wrong)
  Q_correct <- fit_Q(Y, A, W, method = "glm", Qform = Qform_correct)
  g_wrong   <- fit_g(A, W, method = "glm", gform = gform_wrong,   gbounds = gbounds)
  g_correct <- fit_g(A, W, method = "glm", gform = gform_correct, gbounds = gbounds)

  h_fun <- switch(h_choice, ATT = h_fun_ATT, ATC = h_fun_ATC, ATO = h_fun_ATO, NULL)

  h_spec <- switch(h_choice,
                    tilted = d$h,
                    ATE    = rep(1, length(Y)),
                    ATT    = function(g1W, g0W) g1W,
                    ATC    = function(g1W, g0W) g0W,
                    ATO    = function(g1W, g0W) g1W * g0W)

  ## Truth per distinct g fit (shared by the two scenarios using it).
  truth_correct <- if (!is.null(h_fun)) {
    mc_truth_given_ghat(g_correct$fit, h_fun, mc_n = mc_n, gbounds = gbounds)
  } else true_tau
  truth_wrong <- if (!is.null(h_fun)) {
    mc_truth_given_ghat(g_wrong$fit, h_fun, mc_n = mc_n, gbounds = gbounds)
  } else true_tau

  scenarios <- list(
    "Q misspecified; g correct" = list(Q = Q_wrong,   g = g_correct, truth = truth_correct),
    "Q correct; g misspecified" = list(Q = Q_correct, g = g_wrong,   truth = truth_wrong),
    "Both misspecified"         = list(Q = Q_wrong,   g = g_wrong,   truth = truth_wrong),
    "Both correctly specified"  = list(Q = Q_correct, g = g_correct, truth = truth_correct)
  )

  rows <- lapply(names(scenarios), function(sc_name) {
    Q <- scenarios[[sc_name]]$Q; g <- scenarios[[sc_name]]$g; tt <- scenarios[[sc_name]]$truth

    out <- estimate_tau_h_given_fits(Y, A, h_spec, Qfit = Q, gfit = g, n_boot = n_boot)

    data.frame(
      Scenario = sc_name, H = h_choice, Truth = tt,
      MLE = out$MLE, IPW = out$IPW, AIPW = out$AIPW, TMLE = out$TMLE,
      TMLE_EIF_SE = out$TMLE_se,
      TMLE_EIF_Cover = as.integer(out$TMLE_CI[1] <= tt & tt <= out$TMLE_CI[2]),
      TMLE_Bootstrap_SE = out$Bootstrap_SE,
      TMLE_Bootstrap_Normal_Cover = as.integer(out$Bootstrap_CI_normal[1] <= tt & tt <= out$Bootstrap_CI_normal[2]),
      TMLE_Bootstrap_Percentile_Cover = as.integer(out$Bootstrap_CI_percentile[1] <= tt & tt <= out$Bootstrap_CI_percentile[2])
    )
  })

  do.call(rbind, rows)
}


## ------------------------------------------------------------------
## SECTION 5: DRIVER (one h-choice at a time) + SUMMARY TABLES
## ------------------------------------------------------------------

methods_used <- c("MLE", "IPW", "AIPW", "TMLE")

run_h_study <- function(h_choice, B = 1000, n_boot = 300, n1 = 200, n0 = 300,
                         gbounds = c(0.025, 0.975), mc_n = 10000, seed_offset = 840000) {
  h_choice <- match.arg(h_choice, c("tilted", "ATE", "ATT", "ATC", "ATO"))

  true_tau <- switch(h_choice,
                      tilted = approximate_true_tau_5_2_1(),
                      ATE    = approximate_true_tau_general(h_fun_ATE),
                      NULL)   # ATT/ATC/ATO: recomputed every replication inside one_simulation()

  seeds <- seed_offset + seq_len(B)
  do.call(rbind, lapply(seeds, function(s) {
    set.seed(s)
    one_simulation(h_choice, true_tau = true_tau, n1 = n1, n0 = n0, n_boot = n_boot,
                    gbounds = gbounds, mc_n = mc_n)
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
## SECTION 6: RUN ALL 5 h-choices, truncation ON
## ([0.025, 0.975], matching the tmle package default)
## ------------------------------------------------------------------

gbounds_run <- c(0.025, 0.975)   # truncation ON for this run

out_csv <- function(x, filename) write.csv(x, file.path(results_dir, filename))

results_tilted <- run_h_study("tilted", B = 1000, n_boot = 300, gbounds = gbounds_run)
out_csv(summarize_table3(results_tilted), "table3_tilted_trunc_moderateoverlap.csv")
out_csv(summarize_table4(results_tilted), "table4_tilted_trunc_moderateoverlap.csv")

results_ATE <- run_h_study("ATE", B = 1000, n_boot = 300, gbounds = gbounds_run)
out_csv(summarize_table3(results_ATE), "table3_ATE_trunc_moderateoverlap.csv")
out_csv(summarize_table4(results_ATE), "table4_ATE_trunc_moderateoverlap.csv")

results_ATT <- run_h_study("ATT", B = 1000, n_boot = 300, gbounds = gbounds_run)
out_csv(summarize_table3(results_ATT), "table3_ATT_trunc_moderateoverlap.csv")
out_csv(summarize_table4(results_ATT), "table4_ATT_trunc_moderateoverlap.csv")

results_ATC <- run_h_study("ATC", B = 1000, n_boot = 300, gbounds = gbounds_run)
out_csv(summarize_table3(results_ATC), "table3_ATC_trunc_moderateoverlap.csv")
out_csv(summarize_table4(results_ATC), "table4_ATC_trunc_moderateoverlap.csv")

results_ATO <- run_h_study("ATO", B = 1000, n_boot = 300, gbounds = gbounds_run)
out_csv(summarize_table3(results_ATO), "table3_ATO_trunc_moderateoverlap.csv")
out_csv(summarize_table4(results_ATO), "table4_ATO_trunc_moderateoverlap.csv")
