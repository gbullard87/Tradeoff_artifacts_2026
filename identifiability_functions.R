
# ==============================================================================
# Reaction-norm slope-intercept null test
# ==============================================================================
#
# Purpose
# -------
# Test whether an empirical association between group-specific reaction-norm
# intercepts and slopes is more extreme than expected when the *latent* group
# intercepts and slopes are uncorrelated.
#
# The procedure is a parametric bootstrap:
#   1. Fit response ~ predictor separately within each group.
#   2. Estimate an independent random-intercept/random-slope null model by
#      maximum likelihood while retaining the observed sampling design.
#   3. Simulate many complete response vectors under that null model using the
#      exact predictor values, group sizes, imbalance, and replication in the
#      supplied data.
#   4. Refit the within-group regressions after every simulation.
#   5. Build null distributions of recovered slope-intercept correlations and/or
#      covariances.
#   6. Calculate Monte Carlo two-tailed p-values relative to the center of each
#      induced null distribution.
#
# ==============================================================================

#Tests whether the observed correlation or covariance between group-specific reaction-norm 
#intercepts and slopes is more extreme than expected from the sampling design and estimation 
#error alone. Inputs: data must be a data frame containing a numeric response, a numeric predictor 
#with at least two distinct values per retained group, and a grouping variable; response, predictor, 
#and group may be supplied as bare column names or single character strings. statistic is "correlation", 
#"covariance", or "both"; n_sims is a positive integer giving the number of parametric-bootstrap simulations; 
#center_predictor and test_uncentered_intercept are logical values; predictor_center is NULL or a single finite 
#numeric value; min_group_n is an integer of at least 2; seed is NULL or an integer; optimizer_control is a list 
#passed to optim(); keep_simulated_fits and progress are logical values; and alternative is "two.sided", "less", 
#or "greater". Rows with missing or non-finite response or predictor values are removed, and groups lacking sufficient 
#observations or predictor variation are excluded. Output: an object of class "reaction_norm_null_test" containing the 
#empirical statistic, bootstrap null distribution, Monte Carlo p-value and null-distribution summary, group-specific OLS 
#intercepts and slopes, estimated null-model parameters, predictor-centering information, optimization and data-screening 
#diagnostics, and, optionally, coefficients recovered from every simulated dataset. The object has dedicated print() and plot() methods.

reaction_norm_null_test <- function(
    data,
    response,
    predictor,
    group,
    statistic = c("correlation", "covariance", "both"),
    n_sims = 2000L,
    center_predictor = FALSE,
    predictor_center = NULL,
    min_group_n = 3L,
    seed = NULL,
    optimizer_control = list(maxit = 1000),
    keep_simulated_fits = FALSE,
    #test original intercept after centering if you want
    test_uncentered_intercept = FALSE,
    #what kind of test to use?
    alternative = c("two.sided", "less", "greater"),
    progress = interactive()
) {
  # ---------------------------------------------------------------------------
  # 0. Validate arguments and capture column names
  # ---------------------------------------------------------------------------
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  
  # Accept either bare names, e.g. response = y, or strings, e.g. response = "y".
  # We inspect the unevaluated expression so a bare column name need not exist as
  # a separate object in the calling environment.
  get_col_name <- function(arg_expr) {
    if (is.character(arg_expr) && length(arg_expr) == 1L) return(arg_expr)
    deparse(arg_expr)
  }
  
  response_name  <- get_col_name(substitute(response))
  predictor_name <- get_col_name(substitute(predictor))
  group_name     <- get_col_name(substitute(group))
  
  needed <- c(response_name, predictor_name, group_name)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0L) {
    stop("Columns not found in `data`: ", paste(missing_cols, collapse = ", "))
  }
  
  statistic <- match.arg(statistic)
  alternative <- match.arg(alternative)
  n_sims <- as.integer(n_sims)
  min_group_n <- as.integer(min_group_n)
  if (!is.finite(n_sims) || n_sims < 1L) stop("`n_sims` must be at least 1.")
  if (!is.finite(min_group_n) || min_group_n < 2L) {
    stop("`min_group_n` must be at least 2.")
  }
  if (!is.null(seed)) set.seed(seed)
  
  # Make a small internal data frame with stable names. Rows with missing or
  # non-finite response/predictor values, or missing groups, cannot be fitted.
  dat <- data.frame(
    y = data[[response_name]],
    x_original = data[[predictor_name]],
    group = data[[group_name]],
    stringsAsFactors = FALSE
  )
  
  if (!is.numeric(dat$y)) stop("The response column must be numeric.")
  if (!is.numeric(dat$x_original)) stop("The predictor column must be numeric.")
  
  keep <- is.finite(dat$y) & is.finite(dat$x_original) & !is.na(dat$group)
  n_removed_missing <- sum(!keep)
  dat <- dat[keep, , drop = FALSE]
  if (nrow(dat) == 0L) stop("No complete, finite observations remain.")
  
  # Define the predictor origin. predictor_center has precedence. This is useful
  # because slope-intercept covariance depends strongly on the predictor origin.
  if (!is.null(predictor_center)) {
    if (length(predictor_center) != 1L || !is.finite(predictor_center)) {
      stop("`predictor_center` must be one finite number.")
    }
    x_center <- as.numeric(predictor_center)
  } else if (isTRUE(center_predictor)) {
    x_center <- mean(dat$x_original)
  } else {
    x_center <- 0
  }
  dat$x <- dat$x_original - x_center
  
  recover_original_intercept <- function(intercept, slope, center_value) {
    intercept - slope * center_value
  }
  
  # Retain only groups capable of yielding a slope: enough rows and at least two
  # distinct predictor values. The same retained design is used in every bootstrap.
  split_initial <- split(dat, dat$group, drop = TRUE)
  valid_group <- vapply(
    split_initial,
    function(d) nrow(d) >= min_group_n && length(unique(d$x)) >= 2L,
    logical(1)
  )
  excluded_groups <- names(split_initial)[!valid_group]
  group_data <- split_initial[valid_group]
  if (length(group_data) < 3L) {
    stop("At least three valid groups are required to estimate an association.")
  }
  
  # Recombine and create integer group IDs. `group_data` preserves each group's
  # original row-level predictor pattern, including imbalance and replication.
  dat <- do.call(rbind, group_data)
  rownames(dat) <- NULL
  group_labels <- names(group_data)
  n_groups <- length(group_data)
  
  # ---------------------------------------------------------------------------
  # 1. Helper: fit separate OLS reaction norms and extract intercepts/slopes
  # ---------------------------------------------------------------------------
  fit_group_lines <- function(y_by_group, x_by_group, labels) {
    out <- matrix(NA_real_, nrow = length(labels), ncol = 2L,
                  dimnames = list(NULL, c("intercept", "slope")))
    for (g in seq_along(labels)) {
      Xg <- cbind(1, x_by_group[[g]])
      # qr.coef is numerically stable and avoids repeatedly constructing lm objects.
      out[g, ] <- qr.coef(qr(Xg), y_by_group[[g]])
    }
    data.frame(
      group = labels,
      intercept = out[, "intercept"],
      slope = out[, "slope"],
      stringsAsFactors = FALSE
    )
  }
  
  x_by_group <- lapply(group_data, `[[`, "x")
  y_by_group <- lapply(group_data, `[[`, "y")
  empirical_fits <- fit_group_lines(y_by_group, x_by_group, group_labels)
  
  # Defensive check. This should not fail after the group validity screening.
  finite_fit <- is.finite(empirical_fits$intercept) & is.finite(empirical_fits$slope)
  if (sum(finite_fit) < 3L) stop("Fewer than three finite group fits were obtained.")
  empirical_fits <- empirical_fits[finite_fit, , drop = FALSE]
  
  empirical_intercepts <- empirical_fits$intercept
  
  if (test_uncentered_intercept) {
    empirical_intercepts <- recover_original_intercept(
      empirical_fits$intercept,
      empirical_fits$slope,
      x_center
    )
  }
  
  empirical_correlation <- cor(
    empirical_intercepts,
    empirical_fits$slope
  )
  
  empirical_covariance <- cov(
    empirical_intercepts,
    empirical_fits$slope
  )
  
  # ---------------------------------------------------------------------------
  # 2. Fit the independent latent intercept/slope null by maximum likelihood
  # ---------------------------------------------------------------------------
  # Model for group g:
  #   y_g = X_g beta + X_g b_g + e_g
  #   b_g ~ N(0, diag(tau_intercept^2, tau_slope^2))
  #   e_g ~ N(0, sigma_e^2 I)
  #
  # Thus y_g marginally follows:
  #   N(X_g beta, X_g D X_g' + sigma_e^2 I)
  # where D is diagonal. The diagonal D explicitly imposes zero latent
  # slope-intercept covariance. We maximize the marginal likelihood by `optim`.
  
  X_by_group <- lapply(x_by_group, function(x) cbind(1, x))
  
  # Starting fixed effects come from an ordinary regression over all retained rows.
  beta_start <- as.numeric(qr.coef(qr(cbind(1, dat$x)), dat$y))
  if (any(!is.finite(beta_start))) stop("Could not obtain starting fixed effects.")
  
  # Obtain stable positive starting variances. The residual start is pooled across
  # group-specific regressions. 
  rss <- 0
  residual_df <- 0L
  for (g in seq_along(group_data)) {
    Xg <- X_by_group[[g]]
    yg <- y_by_group[[g]]
    bg <- qr.coef(qr(Xg), yg)
    rss <- rss + sum((yg - drop(Xg %*% bg))^2)
    residual_df <- residual_df + nrow(Xg) - 2L
  }
  total_var <- var(dat$y)
  variance_floor <- max(total_var, 1, na.rm = TRUE) * 1e-10
  sigma2_start <- if (residual_df > 0L) rss / residual_df else total_var / 4
  sigma2_start <- max(sigma2_start, variance_floor)
  tau_a2_start <- max(var(empirical_fits$intercept) / 2, variance_floor)
  tau_b2_start <- max(var(empirical_fits$slope) / 2, variance_floor)
  
  # Parameters are beta0, beta1, log(tau_a^2), log(tau_b^2), log(sigma_e^2).
  # Log-variances guarantee positive variance estimates during optimization.
  neg_log_likelihood <- function(par) {
    beta <- par[1:2]
    variance_values <- exp(par[3:5])
    if (any(!is.finite(variance_values))) return(.Machine$double.xmax / 100)
    D <- diag(variance_values[1:2], nrow = 2L)
    sigma2 <- variance_values[3]
    nll <- 0
    
    for (g in seq_along(X_by_group)) {
      Xg <- X_by_group[[g]]
      yg <- y_by_group[[g]]
      Vg <- Xg %*% D %*% t(Xg) + diag(sigma2, nrow(Xg))
      R <- tryCatch(chol(Vg), error = function(e) NULL)
      if (is.null(R)) return(.Machine$double.xmax / 100)
      residual <- yg - drop(Xg %*% beta)
      # If R'R = V, then ||solve(R', residual)||^2 = residual' V^-1 residual.
      z <- forwardsolve(t(R), residual)
      log_det <- 2 * sum(log(diag(R)))
      nll <- nll + 0.5 * (length(yg) * log(2 * pi) + log_det + sum(z^2))
    }
    nll
  }
  
  start <- c(beta_start,
             log(tau_a2_start), log(tau_b2_start), log(sigma2_start))
  opt <- optim(
    par = start,
    fn = neg_log_likelihood,
    method = "BFGS",
    control = optimizer_control,
    hessian = FALSE
  )
  if (opt$convergence != 0L || !is.finite(opt$value)) {
    warning("Null-model optimization did not report clean convergence: ",
            opt$message %||% paste("convergence code", opt$convergence))
  }
  
  beta_hat <- unname(opt$par[1:2])
  names(beta_hat) <- c("mean_intercept", "mean_slope")
  variances_hat <- exp(opt$par[3:5])
  names(variances_hat) <- c("intercept_variance", "slope_variance",
                            "residual_variance")
  
  # ---------------------------------------------------------------------------
  # 3. Parametric bootstrap under zero latent intercept-slope covariance
  # ---------------------------------------------------------------------------
  null_correlation <- rep(NA_real_, n_sims)
  null_covariance  <- rep(NA_real_, n_sims)
  
  if (isTRUE(keep_simulated_fits)) {
    simulated_fits <- vector("list", n_sims)
  } else {
    simulated_fits <- NULL
  }
  
  progress_points <- if (isTRUE(progress) && n_sims >= 10L) {
    unique(round(seq(1, n_sims, length.out = 11L))[-1L])
  } else integer(0)
  
  for (s in seq_len(n_sims)) {
    # Draw independent latent deviations. Independence is the biological/statistical
    # null; any association among the subsequently recovered OLS coefficients is
    # generated by their shared design and observation error.
    latent_intercept <- rnorm(
      n_groups, mean = beta_hat[1], sd = sqrt(variances_hat[1])
    )
    latent_slope <- rnorm(
      n_groups, mean = beta_hat[2], sd = sqrt(variances_hat[2])
    )
    
    sim_y <- vector("list", n_groups)
    for (g in seq_len(n_groups)) {
      xg <- x_by_group[[g]]
      sim_y[[g]] <- latent_intercept[g] + latent_slope[g] * xg +
        rnorm(length(xg), mean = 0, sd = sqrt(variances_hat[3]))
    }
    
    recovered <- fit_group_lines(sim_y, x_by_group, group_labels)
    good <- is.finite(recovered$intercept) & is.finite(recovered$slope)
    recovered <- recovered[good, , drop = FALSE]
    
    intercepts_to_test <- recovered$intercept
    
    if (test_uncentered_intercept) {
      intercepts_to_test <- recover_original_intercept(
        recovered$intercept,
        recovered$slope,
        x_center
      )
    }
    
    if (nrow(recovered) >= 3L &&
        sd(intercepts_to_test) > 0 &&
        sd(recovered$slope) > 0) {
      
      null_correlation[s] <- cor(
        intercepts_to_test,
        recovered$slope
      )
      
      null_covariance[s] <- cov(
        intercepts_to_test,
        recovered$slope
      )
    }
    if (isTRUE(keep_simulated_fits)) simulated_fits[[s]] <- recovered
    
    if (s %in% progress_points) {
      message(sprintf("Parametric bootstrap: %d%% complete",
                      round(100 * s / n_sims)))
    }
  }
  
  # ---------------------------------------------------------------------------
  # 4. Two-tailed Monte Carlo p-values
  # ---------------------------------------------------------------------------
  # The recovered null distribution can be centred away from zero because OLS
  # intercept and slope errors are mathematically associated. Therefore, extremity
  # is measured relative to the simulated null centre, not automatically zero.
  # The +1 correction prevents an estimated p-value of exactly zero.
  monte_carlo_p <- function(
    observed,
    null_values,
    alternative = c("two.sided", "less", "greater")
  ) {
    alternative <- match.arg(alternative)
    z <- null_values[is.finite(null_values)]
    if (length(z) == 0L) {
      return(NA_real_)
    }
    null_center <- mean(z)
    if (alternative == "two.sided") {
      p <- (1 + sum(
        abs(z - null_center) >=
          abs(observed - null_center)
      )) / (length(z) + 1)
    } else if (alternative == "less") {
      p <- (1 + sum(
        z <= observed
      )) / (length(z) + 1)
    } else {   # greater
      p <- (1 + sum(
        z >= observed
      )) / (length(z) + 1)
    }
    p
  }
  
  
  p_correlation <- monte_carlo_p(
    empirical_correlation,
    null_correlation,
    alternative
  )
  
  p_covariance <- monte_carlo_p(
    empirical_covariance,
    null_covariance,
    alternative
  )
  
  summary_all <- data.frame(
    statistic = c("correlation", "covariance"),
    empirical = c(empirical_correlation, empirical_covariance),
    null_mean = c(mean(null_correlation, na.rm = TRUE),
                  mean(null_covariance, na.rm = TRUE)),
    null_median = c(median(null_correlation, na.rm = TRUE),
                    median(null_covariance, na.rm = TRUE)),
    null_sd = c(sd(null_correlation, na.rm = TRUE),
                sd(null_covariance, na.rm = TRUE)),
    null_q025 = c(unname(quantile(null_correlation, 0.025, na.rm = TRUE)),
                  unname(quantile(null_covariance, 0.025, na.rm = TRUE))),
    null_q975 = c(unname(quantile(null_correlation, 0.975, na.rm = TRUE)),
                  unname(quantile(null_covariance, 0.975, na.rm = TRUE))),
    p_value = c(p_correlation, p_covariance),
    valid_simulations = c(sum(is.finite(null_correlation)),
                          sum(is.finite(null_covariance))),
    stringsAsFactors = FALSE
  )
  
  selected_rows <- switch(
    statistic,
    correlation = 1L,
    covariance = 2L,
    both = 1:2
  )
  summary_selected <- summary_all[selected_rows, , drop = FALSE]
  null_selected <- switch(
    statistic,
    correlation = data.frame(simulation = seq_len(n_sims),
                             correlation = null_correlation),
    covariance = data.frame(simulation = seq_len(n_sims),
                            covariance = null_covariance),
    both = data.frame(simulation = seq_len(n_sims),
                      correlation = null_correlation,
                      covariance = null_covariance)
  )
  
  if (test_uncentered_intercept==TRUE){empirical_fits$intercept_tested <- empirical_intercepts}
  
  result <- list(
    call = match.call(),
    statistic_requested = statistic,
    alternative = alternative,
    summary = summary_selected,
    empirical = c(correlation = empirical_correlation,
                  covariance = empirical_covariance),
    null_distribution = null_selected,
    empirical_group_fits = empirical_fits,
    null_parameters = list(
      mean_intercept = beta_hat[1],
      mean_slope = beta_hat[2],
      intercept_variance = variances_hat[1],
      slope_variance = variances_hat[2],
      latent_intercept_slope_covariance = 0,
      residual_variance = variances_hat[3]
    ),
    predictor_center = x_center,
    diagnostics = list(
      optimizer_convergence = opt$convergence,
      optimizer_message = opt$message,
      negative_log_likelihood = opt$value,
      groups_used = n_groups,
      observations_used = nrow(dat),
      rows_removed_missing_or_nonfinite = n_removed_missing,
      groups_excluded = excluded_groups
    ),
    simulated_group_fits = simulated_fits
  )
  class(result) <- "reaction_norm_null_test"
  result
}

# A small infix helper used only for an optional optimizer message.
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || !nzchar(x)) y else x


# Print a concise result while preserving all detailed components in the object.
print.reaction_norm_null_test <- function(x, digits = 4, ...) {
  cat("Reaction-norm slope-intercept null test\n")
  cat("-------------------------------------------\n")
  cat("Groups used:", x$diagnostics$groups_used, "\n")
  cat("Observations used:", x$diagnostics$observations_used, "\n")
  cat("Predictor centre:", format(x$predictor_center, digits = digits), "\n")
  cat("Optimizer convergence code:", x$diagnostics$optimizer_convergence, "\n\n")
  print(x$summary, digits = digits, row.names = FALSE)
  invisible(x)
}


# Optional base-R plotting method for the requested null distribution(s).
plot.reaction_norm_null_test <- function(x, bins = 30, ...) {
  
  requested <- x$statistic_requested
  
  stats <- if (requested == "both") {
    c("correlation", "covariance")
  } else {
    requested
  }
  
  # Combine the requested null distributions into long format.
  plot_data <- do.call(
    rbind,
    lapply(stats, function(st) {
      
      z <- x$null_distribution[[st]]
      z <- z[is.finite(z)]
      
      data.frame(
        statistic = st,
        value = z,
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(plot_data) <- NULL
  
  # Set readable facet labels and preserve their order.
  plot_data$statistic <- factor(
    plot_data$statistic,
    levels = stats,
    labels = stats
  )
  
  # Create one annotation row for each requested statistic.
  annotation_data <- do.call(
    rbind,
    lapply(stats, function(st) {
      
      z <- x$null_distribution[[st]]
      z <- z[is.finite(z)]
      
      data.frame(
        statistic_name = st,
        statistic = paste("Null", st),
        observed = x$empirical[[st]],
        null_median = median(z),
        stringsAsFactors = FALSE
      )
    })
  )
  
  annotation_data$statistic <- factor(
    annotation_data$statistic_name,
    levels = stats,
    labels = stats
  )
  
  
  
  # Construct and return the ggplot object.
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = value)
  ) +
    ggplot2::geom_histogram(
      bins = bins,
      fill = "grey80",
      color = "grey25",
      linewidth = 0.3
    ) +
    
    # Solid red line for the empirical statistic.
    ggplot2::geom_vline(
      data = annotation_data,
      ggplot2::aes(xintercept = observed),
      color = "firebrick",
      linewidth = 1
    ) +
    
    # Dashed blue line for the null median.
    ggplot2::geom_vline(
      data = annotation_data,
      ggplot2::aes(xintercept = null_median),
      color = "navy",
      linewidth = 1,
      linetype = "dashed"
    ) +
    
    # Show one panel for each statistic when statistic = "both".
    ggplot2::facet_wrap(
      ~ statistic,
      scales = "free_x",
      nrow = 1,
      strip.position = "bottom"
    ) +
    
    ggplot2::labs(
      x = NULL,
      y = "Count",
      color = NULL,
      linetype = NULL
    ) +
    
    ggplot2::theme_minimal(base_size = 12) +
    
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    )
  
  # Add a legend using invisible mapped lines. The actual vertical lines above
  # use annotation-specific x-intercepts.
  legend_data <- data.frame(
    label = c("Empirical", "Null median"),
    x = NA_real_
  )
  
  p <- p +
    ggplot2::geom_vline(
      data = legend_data,
      ggplot2::aes(
        xintercept = x,
        color = label,
        linetype = label
      ),
      linewidth = 1,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_manual(
      name = NULL,
      values = c(
        "Empirical" = "firebrick",
        "Null median" = "navy"
      )
    ) +
    ggplot2::scale_linetype_manual(
      name = NULL,
      values = c(
        "Empirical" = "solid",
        "Null median" = "dashed"
      )
    )
  
  return(p)
}


# EXAMPLE
# This synthetic example has an intentionally *uncorrelated* latent intercept and
# slope distribution. It also uses an unbalanced sampling design, demonstrating
# that the function reproduces the exact observed group/predictor structure.
# 
# if (interactive()) {
#   set.seed(42)
# 
#   n_groups <- 60
#   group_ids <- paste0("line_", seq_len(n_groups))
# 
#   # Independent latent reaction-norm parameters under the data-generating model.
#   true_intercepts <- rnorm(n_groups, mean = 3.0, sd = 0.7)
#   true_slopes <- rnorm(n_groups, mean = -0.06, sd = 0.018)
# 
#   # Build an unbalanced design: each line has observations at 3-5 predictor
#   # values, and each line-by-predictor combination has 2-5 replicates.
#   example_parts <- vector("list", n_groups)
#   for (g in seq_len(n_groups)) {
#     temperatures <- sort(sample(36:40, size = sample(3:5, 1)))
#     rows_g <- do.call(rbind, lapply(temperatures, function(temp) {
#       data.frame(
#         line = group_ids[g],
#         temperature = temp,
#         replicate = seq_len(sample(2:5, 1))
#       )
#     }))
#     rows_g$response <- true_intercepts[g] +
#       true_slopes[g] * rows_g$temperature +
#       rnorm(nrow(rows_g), mean = 0, sd = 0.12)
#     example_parts[[g]] <- rows_g
#   }
#   example_data <- do.call(rbind, example_parts)
#   rownames(example_data) <- NULL
# 
#   # Run the test. Use n_sims = 5000 or more for a final analysis. A fixed seed
#   # makes both null-model simulation and the resulting p-value reproducible.
#   example_result <- reaction_norm_null_test(
#     data = example_data,
#     response = "response",
#     predictor = "temperature",
#     group = "line",
#     statistic = "both",             # "correlation", "covariance", or "both"
#     n_sims = 1000,
#     center_predictor = TRUE,         # intercept is response at mean temperature
#     seed = 123,
#     progress = TRUE
#   )
# 
#   print(example_result)
#   plot(example_result)
# 
#   # Useful returned components:
#   head(example_result$null_distribution)
#   head(example_result$empirical_group_fits)
#   example_result$null_parameters
#   example_result$summary



















#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Mixed effects version

# Parametric-bootstrap test of random intercept-slope correlation/covariance
#
# The null model must retain random intercept and slope variation for the focal
# grouping factor but constrain their covariance to zero. The full model must
# allow that covariance. Both models must be Gaussian lme4::lmer() fits using
# the same observations in the same order.

.extract_re_statistic <- function(model, focal_group, slope_term,
                                  statistic = c("correlation", "covariance")) {
  statistic <- match.arg(statistic)
  vc_list <- lme4::VarCorr(model)

  if (!focal_group %in% names(vc_list)) {
    stop(
      "focal_group '", focal_group,
      "' was not found in VarCorr(model). Available groups: ",
      paste(names(vc_list), collapse = ", ")
    )
  }

  vc <- as.matrix(vc_list[[focal_group]])

  if (!all(c("(Intercept)", slope_term) %in% rownames(vc))) {
    stop(
      "The full/refitted model must contain a random intercept and random ",
      "slope named '", slope_term, "' for group '", focal_group, "'."
    )
  }

  if (statistic == "covariance") {
    value <- vc["(Intercept)", slope_term]
  } else {
    sd_intercept <- sqrt(vc["(Intercept)", "(Intercept)"])
    sd_slope <- sqrt(vc[slope_term, slope_term])
    denominator <- sd_intercept * sd_slope

    if (!is.finite(denominator) || denominator <= 0) {
      return(NA_real_)
    }

    value <- vc["(Intercept)", slope_term] / denominator
  }

  if (is.finite(value)) unname(value) else NA_real_
}

.extract_full_group_coefficients <- function(model, focal_group, slope_term) {
  random_effects <- lme4::ranef(model)[[focal_group]]
  fixed_effects <- lme4::fixef(model)

  if (is.null(random_effects)) {
    stop("focal_group was not found in ranef(full_model).")
  }

  required_random <- c("(Intercept)", slope_term)
  required_fixed <- c("(Intercept)", slope_term)

  if (!all(required_random %in% colnames(random_effects))) {
    stop("The focal group lacks the requested random intercept or slope.")
  }

  if (!all(required_fixed %in% names(fixed_effects))) {
    stop("The fixed effects lack the requested intercept or slope term.")
  }

  data.frame(
    group = rownames(random_effects),
    intercept = unname(
      fixed_effects["(Intercept)"] + random_effects[, "(Intercept)"]
    ),
    slope = unname(
      fixed_effects[slope_term] + random_effects[, slope_term]
    ),
    stringsAsFactors = FALSE
  )
}


#Uses a parametric bootstrap to test whether the random intercept-slope 
#correlation or covariance estimated by a linear mixed model is more extreme 
#than expected under a zero-covariance null model. Inputs: full_model and 
#null_model must be Gaussian lme4::lmer() model objects fitted to exactly the 
#same observations in the same order; the full model must estimate a correlated 
#random intercept and slope for the focal grouping factor, whereas the null model 
#must retain both random-effect variances while constraining their covariance to zero. 
#data must be the corresponding data frame with no additional or omitted rows; 
#focal_group and slope_term are single character strings naming the grouping column 
#and numeric slope-predictor column; n_sim, bins, and maxfun are positive integers; 
#statistic is "correlation" or "covariance"; seed is an integer; tags is a 
#three-element character vector; labels is a two-element character vector giving 
#predictor and response labels; plot_type is "density" or "histogram"; 
#null_x_limits is NULL or two increasing finite numeric values; exclude_singular, 
#reml, and finite_correction are logical values; singular_tol is a positive numeric 
#tolerance; and optimizer, line_color, and observed_color are character strings. 
#Output: a list containing the simulated null distribution, two-sided bootstrap 
#p-value, 95% quantiles of the null distribution, 
#the observed random-effect correlation or covariance, a three-panel patchwork figure 
#showing fitted reaction norms, the joint distribution of estimated group coefficients, 
#and the bootstrap null distribution, and the standalone null-distribution panel.

identifiability_test_mixed <- function(
    full_model,
    null_model,
    data,
    n_sim = 1000L,
    statistic = c("correlation", "covariance"),
    focal_group,
    slope_term,
    seed = 1L,
    tags = c("A", "B", "C"),
    labels = c("Predictor", "Response"),
    plot_type = c("density", "histogram"),
    bins = 40L,
    null_x_limits = NULL,
    exclude_singular = FALSE,
    singular_tol = 1e-4,
    reml = lme4::isREML(full_model),
    finite_correction = TRUE,
    optimizer = "bobyqa",
    maxfun = 200000L,
    line_color = "hotpink",
    observed_color = "blue") {

  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }

  statistic <- match.arg(statistic)
  plot_type <- match.arg(plot_type)

  if (!inherits(full_model, "lmerMod") || !inherits(null_model, "lmerMod")) {
    stop("full_model and null_model must both be fitted with lme4::lmer().")
  }
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  if (length(n_sim) != 1L || !is.finite(n_sim) || n_sim < 1) {
    stop("n_sim must be a positive integer.")
  }
  n_sim <- as.integer(n_sim)

  if (length(tags) != 3L) {
    stop("tags must contain exactly three panel labels.")
  }
  if (length(labels) != 2L) {
    stop("labels must contain predictor and response labels.")
  }
  if (!is.character(focal_group) || length(focal_group) != 1L) {
    stop("focal_group must be one character string.")
  }
  if (!is.character(slope_term) || length(slope_term) != 1L) {
    stop("slope_term must be one character string.")
  }
  if (!slope_term %in% names(data) || !is.numeric(data[[slope_term]])) {
    stop("slope_term must name a numeric column in data.")
  }
  if (!focal_group %in% names(data)) {
    stop("focal_group must name a column in data.")
  }
  if (!is.null(null_x_limits) &&
      (length(null_x_limits) != 2L || any(!is.finite(null_x_limits)) ||
       null_x_limits[1L] >= null_x_limits[2L])) {
    stop("null_x_limits must contain two increasing finite values.")
  }

  full_formula <- stats::formula(full_model)
  null_formula <- stats::formula(null_model)

  response_vars <- all.vars(full_formula[[2L]])
  null_response_vars <- all.vars(null_formula[[2L]])

  if (length(response_vars) != 1L || length(null_response_vars) != 1L ||
      response_vars[[1L]] != null_response_vars[[1L]]) {
    stop(
      "Both models must use the same single, untransformed response column."
    )
  }
  response_name <- response_vars[[1L]]

  if (!response_name %in% names(data)) {
    stop("The response column is absent from data.")
  }
  if (nrow(data) != stats::nobs(full_model) ||
      nrow(data) != stats::nobs(null_model)) {
    stop(
      "data must contain exactly the observations used by both models, in ",
      "the same order. Remove omitted rows before fitting both models."
    )
  }

  # Verify that both fitted models refer to the same response observations.
  full_response <- stats::model.response(stats::model.frame(full_model))
  null_response <- stats::model.response(stats::model.frame(null_model))
  if (!isTRUE(all.equal(as.numeric(full_response), as.numeric(null_response)))) {
    stop("full_model and null_model do not use the same response observations.")
  }

  observed_statistic <- .extract_re_statistic(
    full_model,
    focal_group,
    slope_term,
    statistic
  )
  if (!is.finite(observed_statistic)) {
    stop("The empirical full-model statistic is nonfinite.")
  }

  observed_coefficients <- .extract_full_group_coefficients(
    full_model,
    focal_group,
    slope_term
  )

  predictor_sd <- stats::sd(data[[slope_term]], na.rm = TRUE)
  if (!is.finite(predictor_sd) || predictor_sd <= 0) {
    stop("The slope predictor has zero or nonfinite standard deviation.")
  }

  set.seed(seed)
  null_distribution <- rep(NA_real_, n_sim)
  failed_count <- 0L
  singular_count <- 0L

  for (simulation_index in seq_len(n_sim)) {
    # Simulate one complete response vector under the zero-covariance null.
    simulated_response <- 
      stats::simulate(null_model, nsim = 1)

    if (is.null(simulated_response)) {
      failed_count <- failed_count + 1L
      next
    }

    simulated_data <- data
    simulated_data[[response_name]] <- as.numeric(simulated_response$sim_1)

    # Refit the covariance-allowing model so the null distribution reflects
    # the same estimator used for the empirical statistic.
    simulated_full_model <- 
        lme4::lmer(
          formula = full_formula,
          data = simulated_data,
          REML = reml,
          control = lme4::lmerControl(
            optimizer = optimizer,
            optCtrl = list(maxfun = maxfun),
            check.conv.singular = "ignore"
          )
        )
      

    if (is.null(simulated_full_model)) {
      failed_count <- failed_count + 1L
      next
    }

    singular <- lme4::isSingular(
      simulated_full_model,
      tol = singular_tol
    )

    if (singular) {
      singular_count <- singular_count + 1L
      if (isTRUE(exclude_singular)) {
        next
      }
    }
    
    null_distribution[simulation_index] <-
    .extract_re_statistic(
      simulated_full_model,
      focal_group,
      slope_term,
      statistic
    )


}
  null_distribution <- null_distribution[is.finite(null_distribution)]
  n_valid <- length(null_distribution)

  if (n_valid == 0L) {
    stop("No simulation produced a finite null statistic.")
  }
  if (n_valid < 0.5 * n_sim) {
    warning(
      "Fewer than half of simulations were usable: ",
      n_valid, "/", n_sim, "."
    )
  }

  null_median <- stats::median(null_distribution)

  # Two-sided empirical-tail p-value. This remains valid when the simulated
  # null distribution is asymmetric or not centered exactly at zero.
  if (isTRUE(finite_correction)) {
    p_lower <- (sum(null_distribution <= observed_statistic) + 1) /
      (n_valid + 1)
    p_upper <- (sum(null_distribution >= observed_statistic) + 1) /
      (n_valid + 1)
  } else {
    p_lower <- mean(null_distribution <= observed_statistic)
    p_upper <- mean(null_distribution >= observed_statistic)
  }
  p_value <- min(1, 2 * min(p_lower, p_upper))
  #95% conf int
  confint = quantile(null_distribution, probs = c(0.025, 0.975))

  # Panel A: fitted reaction norms for the focal group.
  predictor_grid <- seq(
    min(data[[slope_term]], na.rm = TRUE),
    max(data[[slope_term]], na.rm = TRUE),
    length.out = 100L
  )

  reaction_norms <- do.call(
    rbind,
    lapply(seq_len(nrow(observed_coefficients)), function(i) {
      data.frame(
        predictor = predictor_grid,
        response = observed_coefficients$intercept[i] +
          observed_coefficients$slope[i] * predictor_grid,
        group = observed_coefficients$group[i]
      )
    })
  )

  fixed_effects <- lme4::fixef(full_model)
  population_line <- data.frame(
    predictor = predictor_grid,
    response = unname(fixed_effects["(Intercept)"]) +
      unname(fixed_effects[slope_term]) * predictor_grid
  )

  p1 <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = reaction_norms,
      ggplot2::aes(x = predictor, y = response, group = group),
      color = line_color,
      alpha = 0.12,
      linewidth = 0.45
    ) +
    ggplot2::geom_line(
      data = population_line,
      ggplot2::aes(x = predictor, y = response),
      color = line_color,
      linewidth = 1.2
    ) +
    ggplot2::labs(
      x = labels[1L],
      y = labels[2L],
      tag = tags[1L]
    ) +
    ggplot2::theme_classic()

  # Panel B: observed joint density of fitted group intercepts and slopes.
  p2 <- ggplot2::ggplot(
    observed_coefficients,
    ggplot2::aes(x = intercept, y = slope)
  ) +
    ggplot2::stat_density_2d(
      ggplot2::aes(fill = ggplot2::after_stat(level)),
      geom = "polygon",
      contour = TRUE,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      color = "black",
      alpha = 0.35,
      size = 1
    ) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(
      x = "Estimated intercept",
      y = "Estimated slope",
      fill = "Density",
      tag = tags[2L]
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none")

  # Panel C: null distribution and empirical statistic.
  null_plot_data <- data.frame(value = null_distribution)

  if (plot_type == "density") {
    p3 <- ggplot2::ggplot(
      null_plot_data,
      ggplot2::aes(x = value)
    ) +
      ggplot2::geom_density(
        fill = "grey80",
        color = "black",
        alpha = 0.85,
        na.rm = TRUE
      ) +
      ggplot2::labs(y = "Density")
  } else {
    p3 <- ggplot2::ggplot(
      null_plot_data,
      ggplot2::aes(x = value)
    ) +
      ggplot2::geom_histogram(
        fill = "grey80",
        color = "black",
        bins = bins
      ) +
      ggplot2::labs(y = "Frequency")
  }

  x_axis_label <- if (statistic == "correlation") {
    "Random intercept-slope correlation"
  } else {
    "Random intercept-slope covariance"
  }

  p3 <- p3 +
    ggplot2::geom_vline(
      xintercept = observed_statistic,
      color = observed_color,
      linewidth = 1.2
    ) +
    ggplot2::geom_vline(
      xintercept = null_median,
      color = "black",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    ggplot2::labs(
      title = paste0(
        "two-sided p = ", format.pval(p_value, digits = 3)
      ),
      subtitle = paste0(
        "B = ", signif(B, 3),
        "; valid = ", n_valid, "/", n_sim,
        "; singular = ", singular_count
      ),
      x = x_axis_label,
      tag = tags[3L]
    ) +
    ggplot2::theme_classic()+
    xlim(c(-1,1))

  if (!is.null(null_x_limits)) {
    p3 <- p3 + ggplot2::coord_cartesian(xlim = null_x_limits)
  }

  combined_plot <- patchwork::wrap_plots(p1, p2, p3, nrow = 1)

  list(
    null_distribution = null_distribution,
    p_value = p_value,
    figure = combined_plot,
    null_figure = p3,
    confint = confint,
    observed_statistic = observed_statistic
  )
}

# Example -----------------------------------------------------------------
#
# full_model <- lme4::lmer(
#   response ~ predictor + (1 + predictor | group),
#   data = dat,
#   REML = TRUE
# )
#
# null_model <- lme4::lmer(
#   response ~ predictor +
#     (1 | group) +
#     (0 + predictor | group),
#   data = dat,
#   REML = TRUE
# )
#
# result <- identifiability_test_mixed(
#   full_model = full_model,
#   null_model = null_model,
#   data = dat,
#   n_sim = 2000,
#   statistic = "correlation",
#   focal_group = "group",
#   slope_term = "predictor",
#   labels = c("Predictor", "Response"),
#   seed = 1,
#   exclude_singular = FALSE
# )
#
# result$p_value
# result$B
# result$figure

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# =============================================================================
# POWER, SENSITIVITY, AND SPECIFICITY ANALYSIS
# FOR reaction_norm_null_test()
#
# Assumptions:
#   1. reaction_norm_null_test() is already defined in the environment.
#   2. Each simulated group has a true intercept and slope drawn from an MVN.
#   3. Observations follow:
#
#        y_ij = intercept_i + slope_i * x_ij + error_ij
#
#   4. The test evaluates whether latent intercept-slope correlation differs
#      from zero after accounting for estimation-induced correlation.
#
# Primary output:
#   A single power-curve panel showing the proportion of significant tests
#   against the true generating intercept-slope correlation.
# =============================================================================

#Estimates the power, false-positive rate, and specificity of reaction_norm_null_test() 
#across specified numbers of reaction norms and generating intercept-slope correlations. 
#Inputs: GROUP_SAMPLE_SIZES is an integer vector with values of at least 3; GENERATING_CORRELATIONS 
#is a numeric vector containing values strictly between −1 and 1, where zero estimates the 
#false-positive rate; MEAN_INTERCEPT and MEAN_SLOPE are finite numeric scalars; INTERCEPT_VARIANCE 
#and SLOPE_VARIANCE are positive numeric scalars; RESIDUAL_VARIANCE is a non-negative numeric scalar; 
#PREDICTOR_VALUES is a numeric vector containing at least two distinct values; and REPLICATES_PER_PREDICTOR, 
#N_OUTER_SIMULATIONS, and N_INNER_NULL_SIMULATIONS are positive integers. TEST_STATISTIC, TEST_ALTERNATIVE, 
#CENTER_PREDICTOR, PREDICTOR_CENTER, and TEST_UNCENTERED_INTERCEPT are passed to reaction_norm_null_test() 
#and must follow its required formatting. SIGNIFICANCE_LEVEL is a numeric probability; MASTER_SEED is an integer; 
#RESULTS_FILE, SUMMARY_FILE, and FIGURE_FILE are character file paths; and RESUME_EXISTING_RUN 
#and VERBOSE are logical values. The function requires reaction_norm_null_test() to be defined 
#before execution. Output: the power-curve ggplot object, returned invisibly or assigned by the user, 
#and three files written to the specified paths: an RDS file containing all simulation-level results, 
#a CSV file summarizing rejection rates with exact 95% binomial confidence intervals and diagnostic 
#rates, and a PNG figure showing the proportion of significant tests across generating 
#correlations and sample sizes.

tradeoff_power_analysis = function(
    # Number of independently sampled reaction norms per simulated experiment.
  GROUP_SAMPLE_SIZES = c(10,20,50,100,200,1000),
  # True latent intercept-slope correlations.
  # rho = 0 estimates the false-positive rate.
  GENERATING_CORRELATIONS = c(-0.8, -0.5, -0.3, 0, 0.3, 0.5, 0.8),
  # --------------------------
  # Generating MVN parameters
  # --------------------------
  # Means of the true intercept and slope distributions.
  MEAN_INTERCEPT,
  MEAN_SLOPE,
  # Variances of the true intercept and slope distributions.
  INTERCEPT_VARIANCE,
  SLOPE_VARIANCE,
  # Observation-level residual variance.
  #Impossible to know what this truly is, but we can have an idea for order of magnitude
  RESIDUAL_VARIANCE,
  # Predictor values at which each reaction norm is sampled.
  PREDICTOR_VALUES = c(),
  # Number of replicate observations at each predictor value, per group.
  REPLICATES_PER_PREDICTOR,
  # Number of independently generated datasets per N-by-rho condition.
  # The final analysis can be computationally expensive because every outer
  # simulation runs a complete parametric bootstrap.
  N_OUTER_SIMULATIONS = 20,
  # With a significance threshold of 0.05, values such as 999 or 1999 give
  # sufficiently fine Monte Carlo p-value resolution.
  #number of runs for the trade-off test
  N_INNER_NULL_SIMULATIONS = 200,
  # Test statistic to use for the primary power analysis.
  # Allowed values depend on reaction_norm_null_test().
  TEST_STATISTIC = "covariance",
  # Significance threshold.
  SIGNIFICANCE_LEVEL = 0.05,
  # Use the two-sided test for the main validation.
  TEST_ALTERNATIVE = "two.sided",
  # These arguments are passed directly to reaction_norm_null_test().
  CENTER_PREDICTOR = FALSE,
  PREDICTOR_CENTER = NULL,
  # If the predictor is centered, this controls whether the test uses the
  # intercept at x = 0 in the original predictor coordinates.
  TEST_UNCENTERED_INTERCEPT = FALSE,
  MASTER_SEED = 83621,
  # Intermediate results are saved after every N-by-rho condition.
  RESULTS_FILE = "reaction_norm_power_results.rds",
  SUMMARY_FILE = "reaction_norm_power_summary.csv",
  FIGURE_FILE  = "reaction_norm_power_panel_A.png",
  # Set TRUE to continue from an existing RESULTS_FILE.
  # Set FALSE to begin a completely new run.
  RESUME_EXISTING_RUN = FALSE,
  # Print iteration-level messages.
  VERBOSE = TRUE){
  
  
  # =============================================================================
  # 2. BASIC VALIDATION OF GLOBAL SETTINGS
  # =============================================================================
  
  if (!exists("reaction_norm_null_test", mode = "function")) {
    stop(
      "reaction_norm_null_test() was not found. ",
      "Define the function before running this script."
    )
  }
  
  if (length(PREDICTOR_VALUES) < 2L) {
    stop("At least two distinct predictor values are required.")
  }
  
  if (length(unique(PREDICTOR_VALUES)) < 2L) {
    stop("PREDICTOR_VALUES must contain at least two distinct values.")
  }
  
  if (REPLICATES_PER_PREDICTOR < 1L) {
    stop("REPLICATES_PER_PREDICTOR must be at least 1.")
  }
  
  if (any(GROUP_SAMPLE_SIZES < 3L)) {
    stop("All GROUP_SAMPLE_SIZES must be at least 3.")
  }
  
  if (any(abs(GENERATING_CORRELATIONS) >= 1)) {
    stop("Generating correlations must be strictly between -1 and 1.")
  }
  
  if (INTERCEPT_VARIANCE <= 0 ||
      SLOPE_VARIANCE <= 0 ||
      RESIDUAL_VARIANCE < 0) {
    stop(
      "Intercept and slope variances must be positive, ",
      "and residual variance must be nonnegative."
    )
  }
  
  if (N_OUTER_SIMULATIONS < 1L) {
    stop("N_OUTER_SIMULATIONS must be at least 1.")
  }
  
  if (N_INNER_NULL_SIMULATIONS < 1L) {
    stop("N_INNER_NULL_SIMULATIONS must be at least 1.")
  }
  
  
  # Your supplied function uses %||% in its optimization warning.
  # Define it if it is not already available.
  if (!exists("%||%", mode = "function")) {
    `%||%` <- function(x, y) {
      if (is.null(x) || length(x) == 0L) y else x
    }
  }
  
  
  # =============================================================================
  # 3. HELPER: DRAW CORRELATED INTERCEPTS AND SLOPES
  # =============================================================================
  
  draw_latent_parameters <- function(
    n_groups,
    mean_intercept,
    mean_slope,
    intercept_variance,
    slope_variance,
    correlation
  ) {
    covariance <- correlation *
      sqrt(intercept_variance) *
      sqrt(slope_variance)
    
    sigma <- matrix(
      c(
        intercept_variance, covariance,
        covariance,         slope_variance
      ),
      nrow = 2L,
      byrow = TRUE
    )
    
    eigenvalues <- eigen(
      sigma,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    
    if (any(eigenvalues <= 0)) {
      stop(
        "The requested intercept-slope covariance matrix is not ",
        "positive definite."
      )
    }
    
    # If R'R = Sigma and Z has independent standard-normal columns,
    # then Z %*% R has covariance Sigma.
    R <- chol(sigma)
    
    z <- matrix(
      rnorm(n_groups * 2L),
      nrow = n_groups,
      ncol = 2L
    )
    
    draws <- sweep(
      z %*% R,
      MARGIN = 2L,
      STATS = c(mean_intercept, mean_slope),
      FUN = "+"
    )
    
    colnames(draws) <- c("true_intercept", "true_slope")
    
    data.frame(
      group = seq_len(n_groups),
      true_intercept = draws[, "true_intercept"],
      true_slope = draws[, "true_slope"],
      stringsAsFactors = FALSE
    )
  }
  
  
  # =============================================================================
  # 4. HELPER: SIMULATE ONE COMPLETE REACTION-NORM DATASET
  # =============================================================================
  
  simulate_reaction_norm_dataset <- function(
    n_groups,
    generating_correlation,
    mean_intercept = MEAN_INTERCEPT,
    mean_slope = MEAN_SLOPE,
    intercept_variance = INTERCEPT_VARIANCE,
    slope_variance = SLOPE_VARIANCE,
    residual_variance = RESIDUAL_VARIANCE,
    predictor_values = PREDICTOR_VALUES,
    replicates_per_predictor = REPLICATES_PER_PREDICTOR
  ) {
    latent <- draw_latent_parameters(
      n_groups = n_groups,
      mean_intercept = mean_intercept,
      mean_slope = mean_slope,
      intercept_variance = intercept_variance,
      slope_variance = slope_variance,
      correlation = generating_correlation
    )
    
    within_group_x <- rep(
      predictor_values,
      each = replicates_per_predictor
    )
    
    group_id <- rep(
      latent$group,
      each = length(within_group_x)
    )
    
    x <- rep(
      within_group_x,
      times = n_groups
    )
    
    true_intercept <- latent$true_intercept[group_id]
    true_slope <- latent$true_slope[group_id]
    
    residual_sd <- sqrt(residual_variance)
    
    y <- true_intercept +
      true_slope * x +
      rnorm(
        length(x),
        mean = 0,
        sd = residual_sd
      )
    
    data.frame(
      group = factor(group_id),
      x = x,
      y = y,
      true_intercept = true_intercept,
      true_slope = true_slope,
      stringsAsFactors = FALSE
    )
  }
  
  
  # =============================================================================
  # 5. HELPER: RUN ONE OUTER SIMULATION
  # =============================================================================
  
  run_one_outer_simulation <- function(
    n_groups,
    generating_correlation,
    outer_simulation,
    simulation_seed
  ) {
    set.seed(simulation_seed)
    
    dat <- simulate_reaction_norm_dataset(
      n_groups = n_groups,
      generating_correlation = generating_correlation
    )
    
    # The outer dataset's realized latent correlation will not equal the requested
    # population correlation exactly, especially when n_groups is small.
    latent_by_group <- dat[
      !duplicated(dat$group),
      c("group", "true_intercept", "true_slope")
    ]
    
    realized_latent_correlation <- cor(
      latent_by_group$true_intercept,
      latent_by_group$true_slope
    )
    
    realized_latent_covariance <- cov(
      latent_by_group$true_intercept,
      latent_by_group$true_slope
    )
    
    test_result <- tryCatch(
      suppressWarnings(
        reaction_norm_null_test(
          data = dat,
          response = y,
          predictor = x,
          group = group,
          statistic = TEST_STATISTIC,
          n_sims = N_INNER_NULL_SIMULATIONS,
          center_predictor = CENTER_PREDICTOR,
          predictor_center = PREDICTOR_CENTER,
          min_group_n = 3L,
          seed = simulation_seed + 1000000L,
          optimizer_control = list(maxit = 1000),
          keep_simulated_fits = FALSE,
          test_uncentered_intercept = TEST_UNCENTERED_INTERCEPT,
          alternative = TEST_ALTERNATIVE,
          progress = FALSE
        )
      ),
      error = function(e) e
    )
    
    if (inherits(test_result, "error")) {
      return(
        data.frame(
          n_groups = n_groups,
          generating_correlation = generating_correlation,
          outer_simulation = outer_simulation,
          simulation_seed = simulation_seed,
          realized_latent_correlation = realized_latent_correlation,
          realized_latent_covariance = realized_latent_covariance,
          empirical_statistic = NA_real_,
          null_mean = NA_real_,
          null_sd = NA_real_,
          p_value = NA_real_,
          significant = NA,
          optimizer_convergence = NA_integer_,
          valid_inner_simulations = NA_integer_,
          success = FALSE,
          error_message = conditionMessage(test_result),
          stringsAsFactors = FALSE
        )
      )
    }
    
    result_row <- test_result$summary[
      test_result$summary$statistic == TEST_STATISTIC,
      ,
      drop = FALSE
    ]
    
    if (nrow(result_row) != 1L) {
      return(
        data.frame(
          n_groups = n_groups,
          generating_correlation = generating_correlation,
          outer_simulation = outer_simulation,
          simulation_seed = simulation_seed,
          realized_latent_correlation = realized_latent_correlation,
          realized_latent_covariance = realized_latent_covariance,
          empirical_statistic = NA_real_,
          null_mean = NA_real_,
          null_sd = NA_real_,
          p_value = NA_real_,
          significant = NA,
          optimizer_convergence = test_result$diagnostics$optimizer_convergence,
          valid_inner_simulations = NA_integer_,
          success = FALSE,
          error_message = "Requested statistic was not found in test summary.",
          stringsAsFactors = FALSE
        )
      )
    }
    
    p_value <- result_row$p_value
    
    data.frame(
      n_groups = n_groups,
      generating_correlation = generating_correlation,
      outer_simulation = outer_simulation,
      simulation_seed = simulation_seed,
      realized_latent_correlation = realized_latent_correlation,
      realized_latent_covariance = realized_latent_covariance,
      empirical_statistic = result_row$empirical,
      null_mean = result_row$null_mean,
      null_sd = result_row$null_sd,
      p_value = p_value,
      significant = is.finite(p_value) &&
        p_value < SIGNIFICANCE_LEVEL,
      optimizer_convergence =
        test_result$diagnostics$optimizer_convergence,
      valid_inner_simulations = result_row$valid_simulations,
      success = is.finite(p_value),
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  
  
  # =============================================================================
  # 6. CONSTRUCT THE CONDITION GRID
  # =============================================================================
  
  condition_grid <- expand.grid(
    n_groups = GROUP_SAMPLE_SIZES,
    generating_correlation = GENERATING_CORRELATIONS,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  condition_grid <- condition_grid[
    order(
      condition_grid$n_groups,
      condition_grid$generating_correlation
    ),
    ,
    drop = FALSE
  ]
  
  condition_grid$condition_id <- seq_len(nrow(condition_grid))
  
  
  # =============================================================================
  # 7. LOAD EXISTING RESULTS, IF REQUESTED
  # =============================================================================
  
  if (RESUME_EXISTING_RUN && file.exists(RESULTS_FILE)) {
    simulation_results <- readRDS(RESULTS_FILE)
    
    message(
      "Loaded ",
      nrow(simulation_results),
      " existing simulation results from ",
      RESULTS_FILE,
      "."
    )
  } else {
    simulation_results <- data.frame(
      n_groups = integer(0),
      generating_correlation = numeric(0),
      outer_simulation = integer(0),
      simulation_seed = integer(0),
      realized_latent_correlation = numeric(0),
      realized_latent_covariance = numeric(0),
      empirical_statistic = numeric(0),
      null_mean = numeric(0),
      null_sd = numeric(0),
      p_value = numeric(0),
      significant = logical(0),
      optimizer_convergence = integer(0),
      valid_inner_simulations = integer(0),
      success = logical(0),
      error_message = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  
  # =============================================================================
  # 8. RUN THE OUTER SIMULATION STUDY
  # =============================================================================
  
  set.seed(MASTER_SEED)
  
  # Assign a fixed seed to every condition-by-repetition combination.
  seed_grid <- expand.grid(
    condition_id = condition_grid$condition_id,
    outer_simulation = seq_len(N_OUTER_SIMULATIONS),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  seed_grid$simulation_seed <- sample.int(
    .Machine$integer.max,
    size = nrow(seed_grid),
    replace = FALSE
  )
  
  for (condition_index in seq_len(nrow(condition_grid))) {
    
    current_n <- condition_grid$n_groups[condition_index]
    current_rho <- condition_grid$generating_correlation[condition_index]
    current_condition_id <- condition_grid$condition_id[condition_index]
    
    already_completed <- simulation_results$outer_simulation[
      simulation_results$n_groups == current_n &
        abs(
          simulation_results$generating_correlation - current_rho
        ) < 1e-12
    ]
    
    simulations_to_run <- setdiff(
      seq_len(N_OUTER_SIMULATIONS),
      already_completed
    )
    
    if (length(simulations_to_run) == 0L) {
      if (VERBOSE) {
        message(
          "Skipping completed condition: N = ",
          current_n,
          ", rho = ",
          current_rho
        )
      }
      next
    }
    
    if (VERBOSE) {
      message(
        "\nStarting condition ",
        condition_index,
        " of ",
        nrow(condition_grid),
        ": N = ",
        current_n,
        ", rho = ",
        current_rho,
        ", remaining outer simulations = ",
        length(simulations_to_run)
      )
    }
    
    condition_results <- vector(
      "list",
      length(simulations_to_run)
    )
    
    for (k in seq_along(simulations_to_run)) {
      
      outer_index <- simulations_to_run[k]
      
      simulation_seed <- seed_grid$simulation_seed[
        seed_grid$condition_id == current_condition_id &
          seed_grid$outer_simulation == outer_index
      ]
      
      condition_results[[k]] <- run_one_outer_simulation(
        n_groups = current_n,
        generating_correlation = current_rho,
        outer_simulation = outer_index,
        simulation_seed = simulation_seed
      )
      
      if (VERBOSE &&
          (k == 1L ||
           k %% 10L == 0L ||
           k == length(simulations_to_run))) {
        message(
          "  Completed ",
          k,
          " of ",
          length(simulations_to_run),
          " remaining simulations."
        )
      }
    }
    
    condition_results <- do.call(
      rbind,
      condition_results
    )
    
    simulation_results <- rbind(
      simulation_results,
      condition_results
    )
    
    simulation_results <- simulation_results[
      order(
        simulation_results$n_groups,
        simulation_results$generating_correlation,
        simulation_results$outer_simulation
      ),
      ,
      drop = FALSE
    ]
    
    saveRDS(
      simulation_results,
      file = RESULTS_FILE
    )
    
    if (VERBOSE) {
      message(
        "Saved intermediate results to ",
        RESULTS_FILE,
        "."
      )
    }
  }
  
  
  # =============================================================================
  # 9. SUMMARIZE POWER AND FALSE-POSITIVE RATES
  # =============================================================================
  
  # Exact binomial confidence intervals are used for the estimated rejection rate.
  binomial_summary <- function(significant, success) {
    
    usable <- isTRUE(success) | success %in% TRUE
    significant <- significant[usable]
    significant <- significant[!is.na(significant)]
    
    n_valid <- length(significant)
    n_significant <- sum(significant)
    
    if (n_valid == 0L) {
      return(
        c(
          n_valid = 0,
          n_significant = NA,
          rejection_rate = NA,
          confidence_low = NA,
          confidence_high = NA
        )
      )
    }
    
    interval <- binom.test(
      x = n_significant,
      n = n_valid,
      conf.level = 0.95
    )$conf.int
    
    c(
      n_valid = n_valid,
      n_significant = n_significant,
      rejection_rate = n_significant / n_valid,
      confidence_low = interval[1],
      confidence_high = interval[2]
    )
  }
  
  summary_split <- split(
    simulation_results,
    interaction(
      simulation_results$n_groups,
      simulation_results$generating_correlation,
      drop = TRUE
    )
  )
  
  power_summary <- lapply(
    summary_split,
    function(d) {
      binomial_results <- binomial_summary(
        significant = d$significant,
        success = d$success
      )
      
      data.frame(
        n_groups = unique(d$n_groups),
        generating_correlation =
          unique(d$generating_correlation),
        n_attempted = nrow(d),
        n_valid = unname(binomial_results["n_valid"]),
        n_significant =
          unname(binomial_results["n_significant"]),
        rejection_rate =
          unname(binomial_results["rejection_rate"]),
        confidence_low =
          unname(binomial_results["confidence_low"]),
        confidence_high =
          unname(binomial_results["confidence_high"]),
        failure_rate = mean(!d$success),
        optimizer_nonconvergence_rate = mean(
          d$optimizer_convergence != 0,
          na.rm = TRUE
        ),
        mean_realized_latent_correlation = mean(
          d$realized_latent_correlation,
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  )
  
  power_summary <- do.call(
    rbind,
    power_summary
  )
  
  rownames(power_summary) <- NULL
  
  power_summary <- power_summary[
    order(
      power_summary$n_groups,
      power_summary$generating_correlation
    ),
    ,
    drop = FALSE
  ]
  
  power_summary$result_type <- ifelse(
    abs(power_summary$generating_correlation) < 1e-12,
    "False-positive rate",
    "Power"
  )
  
  write.csv(
    power_summary,
    file = SUMMARY_FILE,
    row.names = FALSE
  )
  
  print(power_summary)
  
  
  # =============================================================================
  # 10. REPORT NULL-CALIBRATION RESULTS
  # =============================================================================
  
  null_summary <- power_summary[
    abs(power_summary$generating_correlation) < 1e-12,
    ,
    drop = FALSE
  ]
  
  if (nrow(null_summary) > 0L) {
    null_summary$specificity <- 1 - null_summary$rejection_rate
    
    message("\nFalse-positive rate and specificity when generating rho = 0:")
    
    print(
      null_summary[
        ,
        c(
          "n_groups",
          "n_valid",
          "rejection_rate",
          "specificity",
          "confidence_low",
          "confidence_high",
          "failure_rate"
        )
      ]
    )
  }
  
  
  # =============================================================================
  # 11. PANEL A: POWER CURVES
  # =============================================================================
  
  # This section uses ggplot2 only for the final figure.
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "The ggplot2 package is required for the figure. ",
      "Install it using install.packages('ggplot2')."
    )
  }
  
  plot_data <- power_summary[
    is.finite(power_summary$rejection_rate),
    ,
    drop = FALSE
  ]
  
  plot_data$n_groups_label <- factor(
    paste0("N = ", plot_data$n_groups),
    levels = paste0("N = ", sort(unique(plot_data$n_groups)))
  )
  
  power_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = generating_correlation,
      y = rejection_rate,
      color = n_groups_label,
      group = n_groups_label
    )
  ) +
    ggplot2::geom_hline(
      yintercept = SIGNIFICANCE_LEVEL,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.5
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dotted",
      color = "grey55",
      linewidth = 0.5
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = confidence_low,
        ymax = confidence_high
      ),
      width = 0.025,
      linewidth = 0.55
    ) +
    ggplot2::geom_line(
      linewidth = 0.9
    ) +
    ggplot2::geom_point(
      size = 2.5
    ) +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(plot_data$generating_correlation)),
      limits = range(plot_data$generating_correlation)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      x = "Generating intercept-slope correlation",
      y = "Proportion of significant tests",
      color = "Reaction norms"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.grid = ggplot2::element_blank()
    )
  
  print(power_plot)
  
  ggplot2::ggsave(
    filename = FIGURE_FILE,
    plot = power_plot,
    width = 7.2,
    height = 5.2,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  
  message(
    "\nAnalysis complete.\n",
    "Raw simulation results: ", RESULTS_FILE, "\n",
    "Power summary: ", SUMMARY_FILE, "\n",
    "Panel A figure: ", FIGURE_FILE
  )
  return(power_plot)
}


