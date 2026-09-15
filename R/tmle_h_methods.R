## ==================================================================
## tmle_h_methods.R
##
## The "methods" module for the TMLE-ECT framework 
##



## ------------------------------------------------------------------
## LAYER 1: PURE ESTIMATORS
## ------------------------------------------------------------------

## MLE 
mle_h <- function(h, Q1W, Q0W) {
  sum(h * (Q1W - Q0W)) / sum(h)
}

## IPW
ipw_h <- function(Y, A, h, g1W, g0W = 1 - g1W) {
  sum(h * (A * Y / g1W - (1 - A) * Y / g0W)) / sum(h)
}

## AIPW
## estimate plus an EIF-based SE (same EIF form used by TMLE below).
aipw_h <- function(Y, A, h, Q1W, Q0W, g1W, g0W = 1 - g1W) {
  N <- length(Y)
  h_star <- h / mean(h)
  terms <- Q1W - Q0W + A / g1W * (Y - Q1W) - (1 - A) / g0W * (Y - Q0W)
  estimate <- mean(h_star * terms)
  EIF <- h_star * (terms - estimate)
  list(estimate = estimate, se = sqrt(mean(EIF^2) / N))
}

## Shared TMLE fluctuation-fitting step: the target-step logistic
## regression logit[Q(X,A;eps)] = logit[Qhat(X,A)] + eps*H(X,A)
fit_eps_h <- function(Y, H, logitQAW, indices = seq_along(Y)) {
  dat <- data.frame(Y = Y[indices], H = H[indices], off = logitQAW[indices])
  fit <- suppressWarnings(glm(Y ~ -1 + H + offset(off), data = dat, family = binomial()))
  eps <- unname(coef(fit)["H"])
  if (!is.finite(eps)) eps <- 0     # guards against degenerate/separated (bootstrap) fits
  eps
}

## TMLE: initial step + target step

tmle_h <- function(Y, A, h, QAW, Q1W, Q0W, g1W, g0W = 1 - g1W) {
  N <- length(Y)
  h_star <- h / mean(h)
  HAW <- A / g1W - (1 - A) / g0W        # base clever covariate
  H <- h_star * HAW                      # H(X,A) of Sec 3.3.1
  logitQAW <- qlogis(QAW)

  eps <- fit_eps_h(Y, H, logitQAW)

  QAW_star <- plogis(logitQAW + eps * H)
  Q1W_star <- plogis(qlogis(Q1W) + eps * (h_star / g1W))
  Q0W_star <- plogis(qlogis(Q0W) - eps * (h_star / g0W))

  estimate <- sum(h * (Q1W_star - Q0W_star)) / sum(h)
  EIF <- h_star * (Q1W_star - Q0W_star - estimate) + h_star * HAW * (Y - QAW_star)
  se <- sqrt(mean(EIF^2) / N)
  ci <- estimate + c(-1, 1) * qnorm(0.975) * se

  list(estimate = estimate, se = se, CI = ci, eps = eps,
       Q1W_star = Q1W_star, Q0W_star = Q0W_star, EIF = EIF,
       H = H, logitQAW = logitQAW)
}

## Target-step-only stratified bootstrap
tmle_bootstrap_h <- function(Y, A, h, Q1W, Q0W, g1W, g0W = 1 - g1W,
                              tmle_fit, n_boot = 300,
                              strata = A, subset = seq_along(Y)) {
  h_star <- h / mean(h)
  H <- tmle_fit$H; logitQAW <- tmle_fit$logitQAW

  strata_id <- split(subset, strata[subset])   # one index vector per stratum

  boot_est <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    boot_id <- unlist(lapply(strata_id, function(idx) sample(idx, length(idx), replace = TRUE)))
    eps_b <- fit_eps_h(Y, H, logitQAW, indices = boot_id)
    Q1_b <- plogis(qlogis(Q1W) + eps_b * (h_star / g1W))
    Q0_b <- plogis(qlogis(Q0W) - eps_b * (h_star / g0W))
    boot_est[b] <- sum(h * (Q1_b - Q0_b)) / sum(h)
  }

  se <- sd(boot_est)
  list(estimates = boot_est, se = se,
       CI_normal     = tmle_fit$estimate + c(-1, 1) * qnorm(0.975) * se,
       CI_percentile = unname(quantile(boot_est, c(0.025, 0.975))))
}


## ------------------------------------------------------------------
## LAYER 2: NUISANCE FITTING (Q and g), can choose glm or super learner
## ------------------------------------------------------------------

fit_Q <- function(Y, A, W, method = c("glm", "SL"),
                   Qform = NULL, SL.library = c("SL.glm"), ...) {
  method <- match.arg(method)
  W <- as.data.frame(W)

  if (method == "glm") {
    if (is.null(Qform)) Qform <- as.formula(paste("Y ~ A*(", paste(names(W), collapse = "+"), ")"))
    fit <- glm(Qform, data = data.frame(Y = Y, A = A, W), family = binomial())
    QAW <- predict(fit, newdata = data.frame(A = A, W), type = "response")
    Q1W <- predict(fit, newdata = data.frame(A = 1, W), type = "response")
    Q0W <- predict(fit, newdata = data.frame(A = 0, W), type = "response")
  } else {
    WA <- data.frame(A = A, W)
    fit <- SuperLearner::SuperLearner(Y = Y, X = WA, family = binomial(),
                                       SL.library = SL.library, ...)
    QAW <- as.numeric(predict(fit, newdata = WA)$pred)
    Q1W <- as.numeric(predict(fit, newdata = data.frame(A = 1, W))$pred)
    Q0W <- as.numeric(predict(fit, newdata = data.frame(A = 0, W))$pred)
  }

  list(QAW = QAW, Q1W = Q1W, Q0W = Q0W, fit = fit, method = method)
}

## gbounds = c(0.025, 0.975) is default truncates; 
fit_g <- function(A, W, method = c("glm", "SL"),
                   gform = NULL, SL.library = c("SL.glm"), gbounds = c(0.025, 0.975), ...) {
  method <- match.arg(method)
  W <- as.data.frame(W)

  if (method == "glm") {
    if (is.null(gform)) gform <- as.formula(paste("A ~", paste(names(W), collapse = "+")))
    fit <- glm(gform, data = data.frame(A = A, W), family = binomial())
    g1W <- predict(fit, newdata = W, type = "response")
  } else {
    fit <- SuperLearner::SuperLearner(Y = A, X = W, family = binomial(),
                                       SL.library = SL.library, ...)
    g1W <- as.numeric(predict(fit, newdata = W)$pred)
  }

  g1W <- pmin(pmax(g1W, gbounds[1]), gbounds[2])
  list(g1W = g1W, g0W = 1 - g1W, fit = fit, method = method)
}

## Predicts P(A=1|X) from either a glm or a SuperLearner g-fit object.
## predict.SuperLearner doesn't accept type="response" and returns a
## list (unlike predict.glm) -- this dispatch lets any caller (e.g. a
## simulation's own Monte Carlo truth calculation) predict from a
## fit_g() result without caring which method produced it.
predict_ghat <- function(g_fit, newdata) {
  if (inherits(g_fit, "SuperLearner")) {
    as.numeric(predict(g_fit, newdata = newdata)$pred)
  } else {
    as.numeric(predict(g_fit, newdata = newdata, type = "response"))
  }
}


## ------------------------------------------------------------------
## LAYER 3: THE USER-FACING ENTRY POINT
## ------------------------------------------------------------------

## Core estimation step given ALREADY-FITTED Q and g (as returned by
## fit_Q()/fit_g()) plus an h(x) specification. estimate_tau_h() below
## calls this after fitting Q/g from raw data; exposed separately so
## that anything re-using the same Q/g fit across multiple h's (e.g. a
## simulation study looping over h-choices under one nuisance-model
## scenario) doesn't have to refit.
estimate_tau_h_given_fits <- function(Y, A, h, Qfit, gfit, n_boot = 300, boot_strata = NULL) {
  h_vec <- if (is.function(h)) h(gfit$g1W, gfit$g0W) else h
  stopifnot(length(h_vec) == length(Y))
  strata <- if (is.null(boot_strata)) A else boot_strata

  mle_est  <- mle_h(h_vec, Qfit$Q1W, Qfit$Q0W)
  ipw_est  <- ipw_h(Y, A, h_vec, gfit$g1W)
  aipw_out <- aipw_h(Y, A, h_vec, Qfit$Q1W, Qfit$Q0W, gfit$g1W)
  tmle_out <- tmle_h(Y, A, h_vec, Qfit$QAW, Qfit$Q1W, Qfit$Q0W, gfit$g1W)
  boot_out <- tmle_bootstrap_h(Y, A, h_vec, Qfit$Q1W, Qfit$Q0W, gfit$g1W,
                                tmle_fit = tmle_out, n_boot = n_boot, strata = strata)

  list(h = h_vec,
       MLE = mle_est,
       IPW = ipw_est,
       AIPW = aipw_out$estimate, AIPW_se = aipw_out$se,
       TMLE = tmle_out$estimate, TMLE_se = tmle_out$se, TMLE_CI = tmle_out$CI,
       Bootstrap_SE = boot_out$se,
       Bootstrap_CI_normal = boot_out$CI_normal,
       Bootstrap_CI_percentile = boot_out$CI_percentile,
       Qfit = Qfit, gfit = gfit, tmle_fit = tmle_out, boot_fit = boot_out)
}

## The single entry point: given raw data (Y, A, W), an h(x)
## specification (fixed vector, or function(g1W, g0W)), and choices
## for how Q0/g0 are fit, this fits Q and g, builds h, and returns
## MLE/IPW/AIPW/TMLE point estimates plus TMLE's EIF-Wald and
## target-step bootstrap inference.
##
##   truncate = TRUE  (default) -> gbounds as given (default [0.025,0.975])
##   truncate = FALSE            -> gbounds forced to c(0,1), i.e. no truncation
##
## ... is passed through to fit_Q()/fit_g() (e.g. cvControl for SL).
estimate_tau_h <- function(Y, A, W, h,
                            Q_method = c("glm", "SL"), Qform = NULL, Q_SL.library = c("SL.glm"),
                            g_method = c("glm", "SL"), gform = NULL, g_SL.library = c("SL.glm"),
                            truncate = TRUE, gbounds = c(0.025, 0.975),
                            n_boot = 300, boot_strata = NULL, ...) {
  Q_method <- match.arg(Q_method); g_method <- match.arg(g_method)
  if (!truncate) gbounds <- c(0, 1)

  Qfit <- fit_Q(Y, A, W, method = Q_method, Qform = Qform, SL.library = Q_SL.library, ...)
  gfit <- fit_g(A, W, method = g_method, gform = gform, SL.library = g_SL.library,
                gbounds = gbounds, ...)

  estimate_tau_h_given_fits(Y, A, h, Qfit, gfit, n_boot = n_boot, boot_strata = boot_strata)
}
