library(ggplot2)
library(tidyr)
library(scales)
library(dplyr)
library(gridExtra)

# Read your data
df_long <- read.csv("~/Downloads/Demand for US Military Intervention - Interviews.csv")

# Create mapping for participant names
participant_names <- setNames(
  paste("Participant", 1:12),
  paste0("P", 1:12)
)

# Define colors for each participant (extending to 12 colors)
participant_colors <- setNames(
  c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "brown",
    "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB"),
  paste0("P", 1:12)
)

# Plot function with line and color
plot_with_line <- function(subject_data, subject_code) {
  subject_name <- participant_names[subject_code]
  subject_color <- participant_colors[subject_code]
  
  ggplot(subject_data, aes(x = Odds_Against, y = Willingness)) +
    geom_point(size = 2, color = subject_color) +
    geom_line(color = subject_color) +
    scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                  labels = function(x) scales::comma(x, accuracy = 0.01)) +
    scale_y_continuous(breaks = c(0, 1), limits = c(0, 1)) +
    labs(title = subject_name,
         x = "Odds Against Success",
         y = "Willingness to support\nUS intervention") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8)) +
    annotation_logticks(sides = "b")
}

# Create plots for all participants
plot_order <- paste0("P", 1:12)
plots <- lapply(plot_order, function(subject_code) {
  data <- df_long[df_long$ID == subject_code, ]
  plot_with_line(data, subject_code)
})

# Display plots (4 columns x 3 rows to fit all 12)
grid.arrange(grobs = plots, ncol = 4)