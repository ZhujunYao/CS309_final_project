# Load required libraries
library(ggplot2)

data<-read.csv("C:/Users/20228/Downloads/Algorithm_data.csv")

# Define colors based on running time
colors <- data$Time

# Define dot size based on neglogMSE
dot_sizes <- -log(data$MSE)

# Create dot plot
ggplot(data, aes(x = R.squared, y = dot_sizes, size = dot_sizes, color = Time)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(name = "Running Time") +
  scale_size_continuous(range = c(1, 10)) +
  labs(x = "R.squared", y = "-log(MSE)", title = "Summary of Algorithm Performance") +
  theme_minimal() +
  theme(legend.position = "right",  # Legend position
        panel.border = element_rect(color = "black", fill = NA, size = 1),  # Panel border
        plot.title = element_text(hjust = 0.5))  # Title alignment
##########################################
# Load required libraries
library(dplyr)
library(ggplot2)


# Summary of the data
summary_data <- data %>%
  group_by(Algorithm) %>%
  summarize(
    Mean_MSE = mean(MSE),
    Mean_R_squared = mean(R.squared),
    Mean_Time = mean(Time)
  )

# Define colors based on running time
colors <- summary_data$Mean_Time

# Define dot size based on neglogMSE
dot_sizes <- -log(summary_data$Mean_MSE)

# Create dot plot
ggplot(summary_data, aes(x = Mean_R_squared, y = Algorithm, color = Mean_Time)) +
#  geom_point(size = 3, alpha = 0.7) +  # Sample dots
  geom_point(aes(size = dot_sizes), alpha = 0.7) +  # MSE dots
  scale_color_gradient(name = "Running Time") +
  labs(x = "Mean R-squared", y = "Algorithm", title = "Summary of Algorithm Performance") +
  scale_size_continuous(range = c(1, 10), guide = FALSE) +  # Disable size legend
  theme_minimal() +
  theme(legend.position = "right",  # Legend position
        panel.border = element_rect(color = "black", fill = NA, size = 1),  # Panel border
        plot.title = element_text(hjust = 0.5))  # Title alignment
# Create dot plot
ggplot(summary_data, aes(x = Mean_R_squared, y = Algorithm, color = Mean_Time)) +
 # geom_point(size = 3, alpha = 0.7) +  # Sample dots
  geom_point(aes(size = dot_sizes), alpha = 0.7) +  # MSE dots
  scale_color_gradient(name = "Running Time") +
  labs(x = "Mean R-squared", y = "Algorithm", title = "Summary of Algorithm Performance") +
  scale_size_continuous(name = "-log(MSE)", range = c(1, 10)) +  # Update size legend
  theme_minimal() +
  theme(legend.position = "right",  # Legend position
        panel.border = element_rect(color = "black", fill = NA, size = 1),  # Panel border
        plot.title = element_text(hjust = 0.5))  # Title alignment

#########################################################
# Load required library
library(dplyr)

# Sample data (replace with your actual data)
data <- data.frame(
  Algorithm = c('Algorithm A', 'Algorithm A', 'Algorithm A', 'Algorithm B', 'Algorithm B', 'Algorithm B'),
  MSE = c(0.1, 0.2, 0.3, 0.15, 0.25, 0.35),
  R_squared = c(0.8, 0.7, 0.6, 0.85, 0.75, 0.65),
  Time = c(10, 12, 15, 8, 11, 14)
)

# Summary of the data
summary_data <- data %>%
  group_by(Algorithm) %>%
  summarize(
    Mean_MSE = mean(MSE),
    Min_MSE = min(MSE),
    Mean_R_squared = mean(R.squared),
    Max_R_squared = max(R.squared),
    Mean_Time = mean(Time),
    Min_Time = min(Time)
  )
summary_data <- as.data.frame(lapply(summary_data, function(x) if(is.numeric(x)) round(x, 4) else x))

# View the summary data
print(summary_data)
setwd("C:/Users/20228/Downloads/")
write.csv(summary_data,file = "summary_algorithm_mean_max_most_recent.csv")
