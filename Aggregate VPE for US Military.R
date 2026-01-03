library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

# Read the data
data <- read.csv("~/Downloads/Subgroup Demand for US Military.csv")

# Create bar chart of VPE vs All
# Order VPE by decreasing likelihood (matching your list order)
vpe_order <- c("Certain", "Almost Certain", "Very Likely", "Probable", 
               "About Even Chance", "Possible", "Unlikely", "Very Unlikely", 
               "Remote", "Impossible")

# Convert VPE to factor with specified order
data$VPE <- factor(data$VPE, levels = vpe_order)

bar_plot <- ggplot(data, aes(x = VPE, y = All, fill = VPE)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "",
       x = "Likelihood of Successful Outcomes from US Military Intervention",
       y = "Proportion of respondents supporting \n US military intervention (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

print(bar_plot)