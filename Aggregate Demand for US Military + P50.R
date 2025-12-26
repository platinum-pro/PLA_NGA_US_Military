library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

# Read the data
data <- read.csv("~/Downloads/Subgroup Demand for US Military.csv")

# Convert All column to numeric if needed
data$All <- as.numeric(data$All)

# GLOBAL k: Fixed at 2 for all analyses
k_global <- 2

cat("######################################\n")
cat("GLOBAL PARAMETER\n")
cat("######################################\n")
cat("k (fixed globally) =", k_global, "\n")
cat("This assumes demand asymptotes to ~1% of Q0\n")
cat("k is treated as a property of the measurement instrument\n\n")

# Define the Koff function with FIXED k
# Note: Now using Odds_Against directly (no price conversion needed)
Koff <- function(odds_against, alpha, Qo, k_fixed) {
  Qo * 10^(k_fixed * (exp(-alpha * Qo * odds_against) - 1))
}

# Function to convert odds against to likelihood percentage
odds_to_likelihood <- function(odds_against) {
  # odds_against = (1 - p) / p
  # Solving for p: p = 1 / (1 + odds_against)
  likelihood <- 100 / (1 + odds_against)
  return(likelihood)
}

# Function to calculate P50 (odds against at 50% support)
calculate_p50 <- function(alpha, Qo, k_fixed, data, target_q = 50) {
  # Add check for valid inputs
  if(alpha <= 0 || Qo <= 0 || k_fixed <= 0) {
    return(NA)
  }
  
  # Define the function we want to find the root of: Q(Odds) - target_q = 0
  objective_function <- function(odds) {
    q_pred <- Qo * 10^(k_fixed * (exp(-alpha * Qo * odds) - 1))
    return(q_pred - target_q)
  }
  
  # Check if Q0 is above or below the target
  q_at_zero <- Qo * 10^(k_fixed * (exp(0) - 1))  # Q at odds = 0
  
  if(q_at_zero < target_q) {
    warning("Q0 is below target Q - P50 cannot be calculated")
    return(NA)
  }
  
  # Define search interval (from minimum to maximum observed odds_against)
  lower_bound <- min(data$Odds_Against[data$Odds_Against > 0], na.rm = TRUE)
  upper_bound <- max(data$Odds_Against, na.rm = TRUE)
  
  # Check if the function crosses the target within our range
  f_lower <- objective_function(lower_bound)
  f_upper <- objective_function(upper_bound)
  
  if(f_lower * f_upper > 0) {
    warning("Target Q is not crossed within the observed odds range")
    return(NA)
  }
  
  # Use uniroot to find P50
  result <- tryCatch({
    uniroot(objective_function, 
            interval = c(lower_bound, upper_bound),
            tol = 0.01)
  }, error = function(e) {
    warning("Error in root-finding: ", e$message)
    return(NULL)
  })
  
  if(!is.null(result)) {
    return(result$root)
  } else {
    return(NA)
  }
}

# Function to fit model and calculate R-squared
fit_and_calculate_rsq <- function(data, k_fixed) {
  # Remove NA values
  valid_data <- data[!is.na(data$All) & data$Odds_Against > 0, ]
  
  fit <- tryCatch({
    nls(formula = All ~ Koff(Odds_Against, alpha, Qo, k_fixed),
        data = valid_data,
        start = list(alpha = 0.0000001, Qo = 100),
        algorithm = "port",
        lower = c(alpha = 0, Qo = 0),
        upper = c(alpha = 0.1, Qo = 100),
        control = nls.control(maxiter = 50000))
  }, error = function(e) {
    message("Error fitting model: ", e$message)
    return(NULL)
  })
  
  if (!is.null(fit)) {
    residuals <- residuals(fit)
    tss <- sum((valid_data$All - mean(valid_data$All, na.rm = TRUE))^2, na.rm = TRUE)
    rss <- sum(residuals^2)
    r_squared <- 1 - (rss / tss)
    
    # Get parameters
    params <- coef(fit)
    
    # Calculate P50 using the fixed k
    p50 <- calculate_p50(params["alpha"], params["Qo"], k_fixed, valid_data, target_q = 50)
    
    cat("R-squared =", round(r_squared, 4), "\n")
    cat("k (fixed) =", k_fixed, "\n")
    cat("alpha =", format(params["alpha"], scientific = TRUE), "\n")
    cat("Q0 =", round(params["Qo"], 2), "\n")
    
    if(!is.na(p50)) {
      likelihood_at_p50 <- odds_to_likelihood(p50)
      cat("P50 (odds against at 50% support) =", round(p50, 3), "\n")
      cat("P50 (as likelihood of success) =", round(likelihood_at_p50, 2), "%\n")
    } else {
      cat("P50 could not be calculated\n")
    }
    
    return(list(fit = fit, r_squared = r_squared, p50 = p50, k = k_fixed))
  }
  
  return(NULL)
}

# Create plot function
create_plot <- function(data, fit_results, k_fixed) {
  # Remove zero or negative odds_against values
  valid_data <- data[data$Odds_Against > 0, ]
  odds_range <- 10^seq(log10(min(valid_data$Odds_Against)), 
                       log10(max(valid_data$Odds_Against)), 
                       length.out = 100)
  plot_data <- data.frame(Odds_Against = odds_range)
  
  if (!is.null(fit_results$fit)) {
    # Generate predictions using the fixed k
    params <- coef(fit_results$fit)
    plot_data$All_pred <- Koff(odds_range, 
                               alpha = params["alpha"], 
                               Qo = params["Qo"], 
                               k_fixed = k_fixed)
  }
  
  p <- ggplot() +
    geom_point(data = valid_data, aes(x = Odds_Against, y = All), shape = 16) +
    geom_line(data = plot_data, aes(x = Odds_Against, y = All_pred))
  
  # Add horizontal line at Q = 50%
  p <- p + geom_hline(yintercept = 50, linetype = "dotted", 
                      color = "grey50", alpha = 0.7)
  
  # Add vertical line for P50 if it exists
  if (!is.null(fit_results$p50) && !is.na(fit_results$p50)) {
    if(fit_results$p50 >= min(valid_data$Odds_Against) && 
       fit_results$p50 <= max(valid_data$Odds_Against)) {
      p <- p + 
        geom_vline(xintercept = fit_results$p50, linetype = "dashed", 
                   color = "blue", alpha = 0.7) +
        annotate("text", x = fit_results$p50, 
                 y = max(valid_data$All, na.rm = TRUE) * 0.95,
                 label = paste0("P50 = ", round(fit_results$p50, 2)),
                 hjust = -0.1, vjust = 1)
    }
  }
  
  p <- p +
    scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000),
                  labels = scales::comma,
                  limits = range(valid_data$Odds_Against)) +
    labs(title = "",
         x = "Odds Against Successful Outcome",
         y = "Proportion of respondents supporting \n US military intervention (%)") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          plot.margin = margin(t = 10, r = 30, b = 10, l = 10, unit = "pt")) +
    annotation_logticks(sides = "b")
  
  return(p)
}

# Fit model and calculate R-squared using global k
fit_results <- fit_and_calculate_rsq(data, k_global)

# Print parameter values
if (!is.null(fit_results$fit)) {
  cat("\nParameters:\n")
  print(summary(fit_results$fit)$parameters)
}

# Create and display plot
p <- create_plot(data, fit_results, k_global)
print(p)