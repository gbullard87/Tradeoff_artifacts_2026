# Parametric-bootstrap identifiability test for a linear mixed model with
# random intercepts and slopes.
#
# The function follows this procedure:
#   1. Extract group-specific intercepts and slopes from the fitted model.
#   2. Compute observed x* = -Cov(intercept, slope) / Var(slope).
#   3. Simulate a full response vector from the fitted model.
#   4. Refit the same model to the simulated response.
#   5. Exclude failed and singular fits.
#   6. Recalculate x* for each valid simulation.
#   7. Compute a two-tailed p-value from |x*_sim| >= |x*_obs|.
#   8. Return three objects:
#        $x_star_distribution
#        $density_plot
#        $p_value

# Internal helper: extract full group-level intercepts and slopes.
.extract_group_coefficients <- function(model, focal_group,
                                        slope_term = NULL) {
  if (!inherits(model, "merMod")) {
    stop("model must be a fitted lme4 mixed-effects model.")
  }
  
  random_list <- lme4::ranef(model)
  
  if (!focal_group %in% names(random_list)) {
    stop(
      "focal_group '", focal_group,
      "' was not found in ranef(model). Available groups: ",
      paste(names(random_list), collapse = ", ")
    )
  }
  
  random_effects <- random_list[[focal_group]]
  fixed_effects <- lme4::fixef(model)
  
  if (!"(Intercept)" %in% colnames(random_effects)) {
    stop("The focal grouping factor does not have a random intercept.")
  }
  
  if (!"(Intercept)" %in% names(fixed_effects)) {
    stop("The model does not have a fixed intercept.")
  }
  
  random_slope_terms <- setdiff(colnames(random_effects), "(Intercept)")
  
  if (is.null(slope_term)) {
    if (length(random_slope_terms) != 1L) {
      stop(
        "The focal grouping factor must contain exactly one random slope, ",
        "or slope_term must be supplied. Available random slopes: ",
        paste(random_slope_terms, collapse = ", ")
      )
    }
    slope_term <- random_slope_terms[[1L]]
  }
  
  if (!slope_term %in% colnames(random_effects)) {
    stop(
      "slope_term '", slope_term,
      "' was not found among the random effects for '", focal_group, "'."
    )
  }
  
  if (!slope_term %in% names(fixed_effects)) {
    stop(
      "slope_term '", slope_term,
      "' was not found among the fixed-effect coefficients."
    )
  }
  
  data.frame(
    group = rownames(random_effects),
    intercept = unname(
      fixed_effects[["(Intercept)"]] + random_effects[["(Intercept)"]]
    ),
    slope = unname(
      fixed_effects[[slope_term]] + random_effects[[slope_term]]
    ),
    stringsAsFactors = FALSE
  )
}

# Main function -----------------------------------------------------------
#
# Args:
#   model:
#     A fitted lme4::lmer model containing a random intercept and random slope
#     for focal_group.
#
#   data:
#     The data frame used to fit model. It must contain the response named on
#     the left side of formula(model).
#
#   focal_group:
#     Character string naming the grouping factor whose full intercepts and
#     slopes define x*. Example: "MaleID".
#
#   slope_term:
#     Character string naming the slope coefficient. Example: "stage".
#     If NULL, the function requires exactly one random slope for focal_group.
#
#   nsim:
#     Number of parametric-bootstrap simulations.
#
#   seed:
#     Random-number seed.
#
#   singular_tol:
#     Tolerance passed to lme4::isSingular().
#
#   exclude_singular:
#     If TRUE, singular refits are excluded
#
#   reml:
#     Whether simulated models are refitted using REML. Defaults to the method
#     used by the supplied code.
#
#   bins:
#     Number of histogram bins.
#
#   x_limits:
#     Optional length-two vector of plot limits. Use c(-5, 5) to reproduce the
#     Great Tit panel. NULL lets ggplot choose the range.
#
#   plot_type:
#     "histogram" reproduces the supplied code. "density" produces a true
#     kernel-density plot.
#
#   finite_correction:
#     FALSE exactly matches mean(abs(sim) >= abs(observed)). TRUE applies the
#     finite-bootstrap correction (extreme + 1)/(valid + 1).
#
# Returns:
#   A list containing exactly:
#     x_star_distribution : numeric vector of valid simulated x* values
#     density_plot         : ggplot object with observed x* as a blue line
#     p_value              : two-tailed bootstrap p-value
#
identifiability_test_mixed <- function(
    model,
    data,
    focal_group,
    slope_term = NULL,
    nsim = 1000L,
    seed = 1L,
    singular_tol = 1e-4,
    exclude_singular = TRUE,
    reml = TRUE,
    bins = 40L,
    x_limits = NULL,
    plot_type = c("histogram", "density"),
    finite_correction = FALSE) {

  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required.")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  plot_type <- match.arg(plot_type)

  if (!inherits(model, "lmerMod")) {
    stop("model must be fitted with lme4::lmer().")
  }

  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }

  if (length(nsim) != 1L || !is.finite(nsim) || nsim < 1L) {
    stop("nsim must be a positive integer.")
  }
  nsim <- as.integer(nsim)

  if (nsim < 999L) {
    warning("Use at least 999 simulations for final inference.")
  }

  if (!is.null(x_limits)) {
    if (length(x_limits) != 2L || any(!is.finite(x_limits)) ||
        x_limits[1L] >= x_limits[2L]) {
      stop("x_limits must contain two increasing finite numbers.")
    }
  }

  # Identify the response column from the fitted model formula.
  response_name <- all.vars(stats::formula(model))[[1L]]

  if (!response_name %in% names(data)) {
    stop(
      "The response column '", response_name,
      "' from formula(model) is absent from data."
    )
  }

  # Observed group-specific coefficients and x*.
  observed_coefficients <- .extract_group_coefficients(
    model = model,
    focal_group = focal_group,
    slope_term = slope_term
  )

  observed_slope_variance <- stats::var(
    observed_coefficients$slope,
    na.rm = TRUE
  )

  if (!is.finite(observed_slope_variance) ||
      observed_slope_variance <= 0) {
    stop("Observed group-specific slope variance is zero or nonfinite.")
  }

  x_star_observed <- -stats::cov(
    observed_coefficients$intercept,
    observed_coefficients$slope,
    use = "complete.obs"
  ) / observed_slope_variance

  if (!is.finite(x_star_observed)) {
    stop("Observed x* is nonfinite.")
  }

  # Parametric bootstrap.
  set.seed(seed)
  x_star_simulated <- rep(NA_real_, nsim)

  for (simulation_index in seq_len(nsim)) {

    # This intentionally matches simulate(model)[[1]] from the supplied code.
    simulated_response <- tryCatch(
      stats::simulate(model, nsim = 1L)[[1L]],
      error = function(e) NULL
    )

    if (is.null(simulated_response)) {
      next
    }

    simulated_data <- data
    simulated_data[[response_name]] <- as.numeric(simulated_response)

    # Refit the original model specification to the simulated response.
    simulated_model <- tryCatch(
      suppressWarnings(
        lme4::lmer(
          formula = stats::formula(model),
          data = simulated_data,
          REML = reml,
          control = lme4::lmerControl(
            check.conv.singular = "ignore"
          )
        )
      ),
      error = function(e) NULL
    )

    if (is.null(simulated_model)) {
      next
    }

    if (isTRUE(exclude_singular) &&
        lme4::isSingular(simulated_model, tol = singular_tol)) {
      next
    }

    simulated_coefficients <- tryCatch(
      .extract_group_coefficients(
        model = simulated_model,
        focal_group = focal_group,
        slope_term = slope_term
      ),
      error = function(e) NULL
    )

    if (is.null(simulated_coefficients)) {
      next
    }

    simulated_slope_variance <- stats::var(
      simulated_coefficients$slope,
      na.rm = TRUE
    )

    if (!is.finite(simulated_slope_variance) ||
        simulated_slope_variance <= 0) {
      next
    }

    x_star_value <- -stats::cov(
      simulated_coefficients$intercept,
      simulated_coefficients$slope,
      use = "complete.obs"
    ) / simulated_slope_variance

    if (is.finite(x_star_value)) {
      x_star_simulated[[simulation_index]] <- x_star_value
    }
  }

  x_star_simulated <- x_star_simulated[is.finite(x_star_simulated)]

  if (length(x_star_simulated) == 0L) {
    stop("No simulations produced a valid x* value.")
  }

  if (length(x_star_simulated) < 0.5 * nsim) {
    warning(
      "Fewer than half of the requested simulations produced valid x* values: ",
      length(x_star_simulated), "/", nsim, "."
    )
  }

  # Two-tailed p-value, matching the supplied absolute-value procedure.
  n_extreme <- sum(abs(x_star_simulated) >= abs(x_star_observed))

  if (isTRUE(finite_correction)) {
    p_value <- (n_extreme + 1) / (length(x_star_simulated) + 1)
  } else {
    p_value <- n_extreme / length(x_star_simulated)
  }

  plot_data <- data.frame(x_star = x_star_simulated)

  if (plot_type == "histogram") {
    result_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = x_star)
    ) +
      ggplot2::geom_histogram(
        fill = "grey80",
        color = "black",
        bins = bins
      )
  } else {
    result_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = x_star)
    ) +
      ggplot2::geom_density(
        fill = "grey80",
        color = "black",
        alpha = 0.8,
        na.rm = TRUE
      )
  }

  result_plot <- result_plot +
    ggplot2::geom_vline(
      xintercept = x_star_observed,
      color = "blue",
      linewidth = 1.2
    ) +
    ggplot2::labs(
      title = paste0("p = ", format.pval(p_value, digits = 3)),
      x = "x*",
      y = if (plot_type == "histogram") "Frequency" else "Density"
    ) +
    ggplot2::theme_classic()

  if (!is.null(x_limits)) {
    # coord_cartesian avoids deleting extreme bootstrap values from the plot's
    # statistical calculation, unlike scale_x_continuous(limits = ...).
    result_plot <- result_plot +
      ggplot2::coord_cartesian(xlim = x_limits)
  }

  # Return exactly the three requested objects.
  list(
    x_star_distribution = x_star_simulated,
    density_plot = result_plot,
    p_value = p_value
  )
}


# Example matching the Great Tit analysis --------------------------------

library(lme4)
library(ggplot2)

dat <- read.csv("Great_Tits.csv")
dat$stage <- dat$NestStage2
dat$y <- dat$TMaleMinDistance

mod <- lmer(
  y ~ stage +
    (1 + stage | MaleID) +
    (1 + stage | BroodID),
  data = dat,
  REML = TRUE
)

test_output <- identifiability_test_mixed(
  model = mod,
  data = dat,
  focal_group = "MaleID",
  slope_term = "stage",
  nsim = 1000,
  seed = 1,
  singular_tol = 1e-4,
  exclude_singular = TRUE,
  reml = TRUE,
  bins = 40,
  x_limits = c(-5, 5),
  plot_type = "histogram",
  finite_correction = FALSE
)

x_star_sim <- test_output$x_star_distribution
p3 <- test_output$density_plot
p_val <- test_output$p_value

print(p3)
p_val
