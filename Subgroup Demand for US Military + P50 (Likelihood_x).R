library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(gridExtra)

# Read the data
data <- read.csv("~/Downloads/Subgroup Demand for US Military.csv")

# Calculate Odds_Against from Likelihood variable in the data
data$Odds_Against <- (100 - data$Likelihood) / data$Likelihood

# Define the demographic pairs to analyze
demographic_pairs <- list(
  list(var_1 = "Older_Age",      var_0 = "Younger_Age",     label1 = "Older Age",      label0 = "Younger Age"),
  list(var_1 = "Male",           var_0 = "Female",          label1 = "Male",           label0 = "Female"),
  list(var_1 = "North_Residence",var_0 = "South_Residence", label1 = "North Residence",label0 = "South Residence"),
  list(var_1 = "High_NG_Trust",  var_0 = "Low_NG_Trust",    label1 = "High NG Trust",  label0 = "Low NG Trust"),
  list(var_1 = "High_US_Trust",  var_0 = "Low_US_Trust",    label1 = "High US Trust",  label0 = "Low US Trust"),
  list(var_1 = "High_Religiosity",var_0 = "Low_Religiosity",label1 = "High Religiosity",label0 = "Low Religiosity")
)

# GLOBAL k: Fixed at 2 for all analyses
k_global <- 2

cat("######################################\n")
cat("GLOBAL PARAMETER\n")
cat("######################################\n")
cat("k (fixed globally) =", k_global, "\n")
cat("This assumes demand asymptotes to ~1% of Q0\n")
cat("k is treated as a property of the measurement instrument\n\n")

# Define the Koff function using Likelihood
Koff_likelihood <- function(likelihood, alpha, Qo, k_fixed) {
  # Convert likelihood back to odds_against for the model
  odds_against <- (100 - likelihood) / likelihood
  Qo * 10^(k_fixed * (exp(-alpha * Qo * odds_against) - 1))
}

# Function to calculate P50 (likelihood at 50% support)
calculate_p50 <- function(alpha, Qo, k_fixed, data_subset, target_q = 50) {
  if (is.na(k_fixed) || alpha <= 0 || Qo <= 0 || k_fixed <= 0) {
    return(NA)
  }
  
  objective_function <- function(likelihood) {
    odds_against <- (100 - likelihood) / likelihood
    q_pred <- Qo * 10^(k_fixed * (exp(-alpha * Qo * odds_against) - 1))
    q_pred - target_q
  }
  
  q_at_100 <- Qo * 10^(k_fixed * (exp(0) - 1))
  
  if (q_at_100 < target_q) {
    warning("Q0 is below target Q - P50 cannot be calculated for this group")
    return(NA)
  }
  
  lower_bound <- max(0.1, min(data_subset$Likelihood[data_subset$Likelihood > 0], na.rm = TRUE))
  upper_bound <- min(99.9, max(data_subset$Likelihood, na.rm = TRUE))
  
  f_lower <- objective_function(lower_bound)
  f_upper <- objective_function(upper_bound)
  
  if (f_lower * f_upper > 0) {
    warning("Target Q is not crossed within the observed likelihood range for this group")
    return(NA)
  }
  
  result <- tryCatch({
    uniroot(objective_function, interval = c(lower_bound, upper_bound), tol = 0.01)
  }, error = function(e) {
    warning("Error in root-finding: ", e$message)
    NULL
  })
  
  if (!is.null(result)) {
    result$root
  } else {
    NA
  }
}

# Function to fit model for a demographic subgroup
fit_demographic_group <- function(data_subset, group_label, k_fixed) {
  valid_data <- data_subset[!is.na(data_subset$value) & data_subset$Likelihood > 0, ]
  
  if (nrow(valid_data) < 3) {
    message("Not enough data for ", group_label)
    return(NULL)
  }
  
  cat("\n=== ", group_label, " ===\n", sep = "")
  
  fit <- tryCatch({
    nls(
      formula = value ~ Koff_likelihood(Likelihood, alpha, Qo, k_fixed),
      data    = valid_data,
      start   = list(alpha = 0.0000001, Qo = 100),
      algorithm = "port",
      lower   = c(alpha = 0,   Qo = 0),
      upper   = c(alpha = 0.1, Qo = 100),
      control = nls.control(maxiter = 50000)
    )
  }, error = function(e) {
    message("Error fitting model for ", group_label, ": ", e$message)
    NULL
  })
  
  if (!is.null(fit)) {
    residuals <- residuals(fit)
    tss <- sum((valid_data$value - mean(valid_data$value, na.rm = TRUE))^2, na.rm = TRUE)
    rss <- sum(residuals^2)
    r_squared <- 1 - (rss / tss)
    
    params <- coef(fit)
    p50 <- calculate_p50(params["alpha"], params["Qo"], k_fixed, valid_data, target_q = 50)
    
    cat("R-squared =", round(r_squared, 4), "\n")
    cat("alpha =", format(params["alpha"], scientific = TRUE), "\n")
    cat("Q0 =", round(params["Qo"], 2), "\n")
    
    if (!is.na(p50)) {
      cat("P50 (likelihood at 50% support) =", round(p50, 2), "%\n")
    } else {
      cat("P50 could not be calculated\n")
    }
    
    return(list(
      fit       = fit,
      r_squared = r_squared,
      p50       = p50,
      label     = group_label,
      alpha     = params["alpha"],
      Qo        = params["Qo"],
      k         = k_fixed
    ))
  }
  
  NULL
}

# Function to create comparison plot with overall curve
create_comparison_plot <- function(data, demographic_pair, fit_results_1, fit_results_0,
                                   k_fixed, overall_data, show_title = TRUE) {
  col_1 <- demographic_pair$var_1
  col_0 <- demographic_pair$var_0
  
  data_1 <- data.frame(Likelihood = data$Likelihood, value = data[[col_1]]) %>%
    filter(!is.na(value) & Likelihood > 0)
  data_0 <- data.frame(Likelihood = data$Likelihood, value = data[[col_0]]) %>%
    filter(!is.na(value) & Likelihood > 0)
  
  likelihood_range <- seq(from = max(data$Likelihood[data$Likelihood > 0]),
                          to = min(data$Likelihood[data$Likelihood > 0]), 
                          length.out = 100)
  plot_data <- data.frame(Likelihood = likelihood_range)
  
  # Overall curve
  if (!is.null(overall_data$fit)) {
    params_overall <- coef(overall_data$fit)
    plot_data$overall_pred <- Koff_likelihood(
      likelihood_range,
      alpha   = params_overall["alpha"],
      Qo      = params_overall["Qo"],
      k_fixed = k_fixed
    )
  }
  
  # Subgroup curves
  if (!is.null(fit_results_1$fit)) {
    params_1 <- coef(fit_results_1$fit)
    plot_data$group1_pred <- Koff_likelihood(
      likelihood_range,
      alpha   = params_1["alpha"],
      Qo      = params_1["Qo"],
      k_fixed = k_fixed
    )
  }
  
  if (!is.null(fit_results_0$fit)) {
    params_0 <- coef(fit_results_0$fit)
    plot_data$group0_pred <- Koff_likelihood(
      likelihood_range,
      alpha   = params_0["alpha"],
      Qo      = params_0["Qo"],
      k_fixed = k_fixed
    )
  }
  
  p <- ggplot()
  
  # Overall curve (background)
  if (!is.null(overall_data$fit)) {
    p <- p +
      geom_line(
        data  = plot_data,
        aes(x = Likelihood, y = overall_pred),
        color    = "black",
        linewidth = 0.6,
        alpha    = 1,
        linetype = "dashed"
      )
  }
  
  # Named color mapping for subgroups
  subgroup_colors <- setNames(
    c("slateblue", "yellow4"),
    c(demographic_pair$label1, demographic_pair$label0)
  )
  
  # Points
  if (nrow(data_1) > 0) {
    p <- p + geom_point(
      data  = data_1,
      aes(x = Likelihood, y = value, color = demographic_pair$label1),
      shape = 16,
      alpha = 0.6,
      size  = 1.5
    )
  }
  
  if (nrow(data_0) > 0) {
    p <- p + geom_point(
      data  = data_0,
      aes(x = Likelihood, y = value, color = demographic_pair$label0),
      shape = 17,
      alpha = 0.6,
      size  = 1.5
    )
  }
  
  # Curves
  if (!is.null(fit_results_1$fit)) {
    p <- p + geom_line(
      data  = plot_data,
      aes(x = Likelihood, y = group1_pred, color = demographic_pair$label1),
      linewidth = 0.8
    )
  }
  
  if (!is.null(fit_results_0$fit)) {
    p <- p + geom_line(
      data  = plot_data,
      aes(x = Likelihood, y = group0_pred, color = demographic_pair$label0),
      linewidth = 0.8
    )
  }
  
  # Horizontal 50 percent line
  p <- p + geom_hline(
    yintercept = 50,
    linetype   = "dotted",
    color      = "grey50",
    alpha      = 0.7
  )
  
  # P50 vertical lines, using names to index colors
  if (!is.null(fit_results_1$p50) && !is.na(fit_results_1$p50)) {
    p <- p + geom_vline(
      xintercept = fit_results_1$p50,
      linetype   = "dashed",
      color      = subgroup_colors[demographic_pair$label1],
      alpha      = 0.4
    )
  }
  
  if (!is.null(fit_results_0$p50) && !is.na(fit_results_0$p50)) {
    p <- p + geom_vline(
      xintercept = fit_results_0$p50,
      linetype   = "dashed",
      color      = subgroup_colors[demographic_pair$label0],
      alpha      = 0.4
    )
  }
  
  # Panel titles
  plot_title <- if (show_title) {
    if (grepl("Age", demographic_pair$label1)) "Age"
    else if (grepl("Male|Female", demographic_pair$label1)) "Gender"
    else if (grepl("Residence", demographic_pair$label1)) "Residence"
    else if (grepl("NG Trust", demographic_pair$label1)) "Trust in Nigerian Government"
    else if (grepl("US Trust", demographic_pair$label1)) "Trust in US Government"
    else if (grepl("Religiosity", demographic_pair$label1)) "Religiosity"
    else ""
  } else ""
  
  p <- p +
    scale_x_continuous(
      trans = "reverse",
      breaks = seq(0, 100, 20),
      limits = c(100, 0)
    ) +
    scale_color_manual(values = subgroup_colors, name = "") +
    labs(
      title = plot_title,
      x     = "Likelihood of successful suppression \nof regional violence (%)",
      y     = "% supporting US military \nintervention"
    ) +
    theme_minimal() +
    theme(
      panel.grid    = element_blank(),
      axis.line     = element_line(color = "black"),
      plot.title    = element_text(hjust = 0.5, size = 10, face = "bold"),
      axis.title    = element_text(size = 8),
      axis.text     = element_text(size = 6),
      legend.position  = "bottom",
      legend.text      = element_text(size = 8),
      legend.key.size  = unit(0.4, "cm"),
      plot.margin      = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
    )
  
  p
}

# Create the "All" fit first (for overlay)
cat("\n\n######################################\n")
cat("Creating 'All' plot\n")
cat("######################################\n")

valid_data_all <- data[!is.na(data$All) & data$Likelihood > 0, ]

fit_all <- tryCatch({
  nls(
    formula = All ~ Koff_likelihood(Likelihood, alpha, Qo, k_global),
    data    = valid_data_all,
    start   = list(alpha = 0.0000001, Qo = 100),
    algorithm = "port",
    lower   = c(alpha = 0, Qo = 0),
    upper   = c(alpha = 0.1, Qo = 100),
    control = nls.control(maxiter = 50000)
  )
}, error = function(e) NULL)

overall_fit_data <- NULL

if (!is.null(fit_all)) {
  params_all <- coef(fit_all)
  p50_all <- calculate_p50(params_all["alpha"], params_all["Qo"], k_global,
                           valid_data_all, target_q = 50)
  
  overall_fit_data <- list(
    fit   = fit_all,
    alpha = params_all["alpha"],
    Qo    = params_all["Qo"],
    p50   = p50_all
  )
}

# Main analysis loop for demographic pairs
all_results <- list()

for (i in seq_along(demographic_pairs)) {
  demo_pair <- demographic_pairs[[i]]
  
  cat("\n\n######################################\n")
  cat("Analyzing:", demo_pair$label1, "vs", demo_pair$label0, "\n")
  cat("######################################\n")
  
  col_1 <- demo_pair$var_1
  col_0 <- demo_pair$var_0
  
  if (!col_1 %in% names(data) || !col_0 %in% names(data)) {
    cat("Columns", col_1, "or", col_0, "not found in data. Skipping.\n")
    next
  }
  
  data_1 <- data.frame(Likelihood = data$Likelihood, value = data[[col_1]])
  data_0 <- data.frame(Likelihood = data$Likelihood, value = data[[col_0]])
  
  fit_1 <- fit_demographic_group(data_1, demo_pair$label1, k_global)
  fit_0 <- fit_demographic_group(data_0, demo_pair$label0, k_global)
  
  if (!is.null(fit_1) || !is.null(fit_0)) {
    p <- create_comparison_plot(
      data         = data,
      demographic_pair = demo_pair,
      fit_results_1    = fit_1,
      fit_results_0    = fit_0,
      k_fixed          = k_global,
      overall_data     = overall_fit_data,
      show_title       = TRUE
    )
    
    comparison_id <- paste(demo_pair$var_1, demo_pair$var_0, sep = "_vs_")
    all_results[[comparison_id]] <- list(
      group1 = fit_1,
      group0 = fit_0,
      plot   = p
    )
  }
}

# Create grid with all demographic plots
grid_plots <- list()

for (result_name in names(all_results)) {
  grid_plots[[length(grid_plots) + 1]] <- all_results[[result_name]]$plot
}

if (length(grid_plots) > 0) {
  grid_combined <- grid.arrange(grobs = grid_plots, ncol = 2, nrow = 3)
  cat("\n\nGrid plot with", length(grid_plots), "panels saved and displayed\n")
}

# Create summary table of P50 values
p50_summary <- data.frame(
  Comparison = character(),
  Group      = character(),
  P50        = numeric(),
  Alpha      = numeric(),
  Q0         = numeric(),
  k          = numeric(),
  R_squared  = numeric(),
  stringsAsFactors = FALSE
)

for (result_name in names(all_results)) {
  result <- all_results[[result_name]]
  
  if (!is.null(result$group1)) {
    p50_summary <- rbind(
      p50_summary,
      data.frame(
        Comparison = result_name,
        Group      = result$group1$label,
        P50        = ifelse(is.na(result$group1$p50), NA, result$group1$p50),
        Alpha      = result$group1$alpha,
        Q0         = result$group1$Qo,
        k          = result$group1$k,
        R_squared  = result$group1$r_squared
      )
    )
  }
  
  if (!is.null(result$group0)) {
    p50_summary <- rbind(
      p50_summary,
      data.frame(
        Comparison = result_name,
        Group      = result$group0$label,
        P50        = ifelse(is.na(result$group0$p50), NA, result$group0$p50),
        Alpha      = result$group0$alpha,
        Q0         = result$group0$Qo,
        k          = result$group0$k,
        R_squared  = result$group0$r_squared
      )
    )
  }
}

cat("\n\n######################################\n")
cat("P50 SUMMARY TABLE\n")
cat("######################################\n")
cat("Note: k = ", k_global, " (fixed globally for all groups)\n", sep = "")
cat("All P50 values are directly comparable on the same scale\n")
cat("P50 is now reported as likelihood percentage at 50% support\n\n")
print(p50_summary, row.names = FALSE)