# Install necessary packages
install.packages("lubridate")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("qcc")

# Load the required libraries
library(dplyr)
library(ggplot2)
library(lubridate)
library(qcc)

# Load R data file containing processed data
misc = load("Data/Processed/misc.Rdata")

# Initialize the station codes dataset
station_codes <- station_codes_agg

# Read the aggregated Origin-Destination (OD) matrix data
OD <- read.csv("Data/Processed/IPF/OD_Trenord_IPF_agg_LM_1_e.csv")

# Remove rows with NA values
OD <- na.omit(OD)

# Create a list to store weekly OD matrices
OD_list = list()

# Extract column names representing weeks
weeks <- colnames(OD)[-c(1, 2)]
weeks <- substring(weeks, 6)  # Remove the first 5 characters (if necessary)

# Populate the OD_list with matrices for each week
for (w in weeks) {
  OD_list[[w]] = OD_week(OD, w, station_codes)  # Function to create weekly OD matrices
}

# Convert each week's data from data frame to matrix
for (w in weeks) {
  OD_list[[w]] = as.matrix(OD_list[[w]])
}

# Calculate the length of the OD_list
t = length(OD_list)

# Calculate the Frobenius norm for each matrix in the OD_list
norm_frobenius_1 <- sapply(OD_list, function(mat) sqrt(sum(mat^2)))

# Compute and visualize the autocorrelation of Frobenius norms
acf_result <- acf(norm_frobenius_1, main = "Autocorrelation of Frobenius Norm in OD Matrices")
print(acf_result)

# Calculate the mean OD matrix
media_OD <- Reduce("+", OD_list) / t

# Compute the norm of the mean matrix
norm_media = sqrt(sum(media_OD^2))

# Calculate the difference between the Frobenius norms and the mean norm
norm_difference = sapply(norm_frobenius_1, function(mat) (mat - norm_media))

# Visualize the autocorrelation of the norm difference
acf_result <- acf(norm_difference, main = "Autocorrelation of Frobenius Norm Differences")
print(acf_result)

# Generate a quantile-quantile (QQ) plot and perform Shapiro-Wilk test for normality
qqnorm(norm_difference)
qqline(norm_difference, col = "red")
shapiro.test(norm_difference)

# Plot the density of the norm differences
plot(density(norm_difference))

# Calculate mean and standard deviation of the norm differences
media_norm1 <- mean(norm_difference)
dev_std_norm1 <- sd(norm_difference)

# Standardize the norm differences
norm_difference = as.vector(norm_difference)
stand = (norm_difference - media_norm1) / dev_std_norm1

# Parameters for CUSUM Chart
h <- 5  # Alert limit
k <- 0.5 # Significant change
head.start <- 0  # Initial value of CUSUM
decision.interval <- h  # Set a decision interval

# Create a CUSUM chart
cusum_chart <- cusum(stand, 1, center = 0, std.dev = 1, head.start = 0, decision.interval = 5, se.shift = 1, "standardized norm differences", plot = TRUE)

# Data frame for plotting moving averages
dati = data.frame(
  week = 1:51,
  valore = norm_difference  # Simulated data
)

# Calculate moving averages and control limits
window_size <- 3
dati <- dati %>% 
  mutate(media_mobile = zoo::rollmean(valore, k = window_size, fill = NA, align = "right"))

# Compute standard deviation for limits
deviazione_standard <- sd(dati$valore, na.rm = TRUE)
dati <- dati %>%
  mutate(
    limite_superiore = media_mobile + 3 * deviazione_standard / sqrt(window_size),
    limite_inferiore = media_mobile - 3 * deviazione_standard / sqrt(window_size),
    limite_superiore1 = media_mobile + deviazione_standard / sqrt(window_size),
    limite_inferiore1 = media_mobile - deviazione_standard / sqrt(window_size),
    limite_superiore2 = media_mobile + 2 * deviazione_standard / sqrt(window_size),
    limite_inferiore2 = media_mobile - 2 * deviazione_standard / sqrt(window_size)
  )
# Adjust limits for NA values in moving average
dati[is.na(dati$media_mobile), "limite_superiore"] = 3 * deviazione_standard / sqrt(dati[is.na(dati$media_mobile), "week"])
dati[is.na(dati$media_mobile), "limite_inferiore"] = -3 * deviazione_standard / sqrt(dati[is.na(dati$media_mobile), "week"])
dati[is.na(dati$media_mobile), "limite_superiore1"] = deviazione_standard / sqrt(dati[is.na(dati$media_mobile), "week"])
dati[is.na(dati$media_mobile), "limite_inferiore1"] = -deviazione_standard / sqrt(dati[is.na(dati$media_mobile), "week"])
dati[is.na(dati$media_mobile), "limite_superiore2"] = 2 * deviazione_standard / sqrt(dati[is.na(dati$media_mobile), "week"])
dati[is.na(dati$media_mobile), "limite_inferiore2"] = -2 * deviazione_standard / sqrt(dati[is.na(dati$media_mobile), "week"])

# Plot the control chart for moving average
plot(dati$valore, type = "o", col = "blue", pch = 16,
     main = "Moving Average Control Chart for Milano-Varese OD matrices",
     xlab = "Week", ylab = "Difference of Frobenius norms",
     ylim = c(min(dati$limite_inferiore), max(dati$limite_superiore)))

# Add limit lines to the plot
lines(dati$limite_superiore1, col = "green", lty = 2)
lines(dati$limite_inferiore1, col = "green", lty = 2)
lines(dati$limite_superiore2, col = "orange", lty = 2)
lines(dati$limite_inferiore2, col = "orange", lty = 2)
lines(dati$limite_superiore, col = "red", lty = 2)
lines(dati$limite_inferiore, col = "red", lty = 2)

# Add a horizontal line for a reference value
abline(h = 86.17416, col = "black", lty = 2)

# Add a legend to the plot
legend("bottomright", legend = c("Difference Frobenius norms", "Zone C", "Warning limits", "Control limits"),
       col = c("blue", "green", "orange", "red"), lty = c(1, 2, 3, 4), pch = c(8, NA))

# Check for violations in CUSUM chart
cat("Violations:", cusum_chart$violations, "\n")

# Calculate upper and lower limits for control based on mean and standard deviation
limite_superiore1 <- media_norm1 + 2 * dev_std_norm1
limite_inferiore1 <- media_norm1 - 2 * dev_std_norm1
limite_superiore2 <- media_norm1 + 3 * dev_std_norm1
limite_inferiore2 <- media_norm1 - 3 * dev_std_norm1
limite_superiore0 <- media_norm1 + dev_std_norm1
limite_inferiore0 <- media_norm1 - dev_std_norm1

# Calculate the control chart for norm differences
cap_ratio = (limite_superiore2 - limite_inferiore2) / (6 * dev_std_norm1)

# Plot the norm differences
plot(norm_difference, type = "o", col = "blue", pch = 16,
     main = "Control Chart for Milano-Varese OD matrices",
     xlab = "Week", ylab = "Difference of Frobenius norms",
     ylim = c(min(norm_difference, limite_inferiore2), max(norm_difference, limite_superiore1, limite_superiore2)))

# Add control limits to the chart
abline(h = limite_superiore1, col = "orange", lty = 2)
abline(h = limite_inferiore1, col = "orange", lty = 2)
abline(h = limite_superiore2, col = "red", lty = 2)
abline(h = limite_inferiore2, col = "red", lty = 2)
abline(h = limite_superiore0, col = "green", lty = 2)
abline(h = limite_inferiore0, col = "green", lty = 2)

# Add a horizontal line for a reference value
abline(h = 86.17416, col = "black", lty = 2)
legend("bottomright", legend = c("Difference Frobenius norms", "Zone C", "Warning limits", "Control limits"),
       col = c("blue", "green", "orange", "red"), lty = c(1, 2, 3, 4), pch = c(8, NA))

# Calculate Frobenius norm for each matrix in OD_list
norm_frobenius <- sapply(OD_list, function(mat) sqrt(sum((mat - media_OD)^2)))

# Visualize the autocorrelation of Frobenius norms
acf_result <- acf(norm_frobenius, main = "Autocorrelation of Frobenius Norm in OD Matrices")
print(acf_result)

# Create a QQ plot for the Frobenius norms
qqnorm(norm_frobenius)
qqline(norm_frobenius, col = "red")

# Perform a Shapiro-Wilk test for normality on the Frobenius norms
shapiro.test(norm_frobenius)

# Plot the density of the Frobenius norms
plot(density(norm_frobenius))

# Calculate the maximum norm for OD matrices
norm_max <- sapply(OD_list, function(mat) max(abs(mat - media_OD)))  # Maximum norm

# Calculate total variation for OD matrices
total_variation <- sapply(OD_list, function(mat) sum(abs(mat - media_OD)))  # Total variation

# Calculate mean and standard deviation of Frobenius norms
media_norm <- mean(norm_frobenius)
dev_std_norm <- sd(norm_frobenius)

# Define control limits for Frobenius norms
limite_superiore <- media_norm + 1 * dev_std_norm
limite_superiore2 <- media_norm + 2 * dev_std_norm
limite_superiore3 <- media_norm + 3 * dev_std_norm

# Calculate mean and standard deviation of max norms
media_max <- mean(norm_max)
dev_std_max <- sd(norm_max)
limite_superiore_max <- media_max + 2 * dev_std_max

# Calculate mean and standard deviation of total variation
media_total_variation <- mean(total_variation)
dev_std_total_variation <- sd(total_variation)
limite_superiore_total_variation <- media_total_variation + 2 * dev_std_total_variation
limite_superiore_total_variation2 <- media_total_variation + 3 * dev_std_total_variation

# Set up the plotting area for multiple plots
par(mfrow = c(3, 1))

# Plot Frobenius norms
plot(norm_frobenius, type = "o", col = "blue", pch = 16,
     main = "Control Chart for Milano-Varese OD matrices",
     xlab = "Week", ylab = "Frobenius norm",
     ylim = c(min(norm_frobenius), max(norm_frobenius, limite_superiore, limite_superiore2, limite_superiore3)))

# Add horizontal lines for control limits
abline(h = limite_superiore, col = "green", lty = 2)
abline(h = limite_superiore2, col = "orange", lty = 2)
abline(h = limite_superiore3, col = "red", lty = 2)
abline(h = media_norm, col = "black", lty = 2)

# Add a legend to the Frobenius norms plot
legend("topright", legend = c("Frobenius norm", "Zone C", "Warning limit", "Control limit"),
       col = c("blue", "green", "orange", "red"), lty = c(1, 2, 3, 4), pch = c(8, NA))

# Plot the maximum norms
plot(norm_max, type = "o", col = "green", pch = 16,
     main = "Control Chart for Milano-Varese OD matrices - Max Norm",
     xlab = "Week", ylab = "Max norm",
     ylim = c(min(norm_max), max(norm_max, limite_superiore_max)))

# Add control limit line for maximum norms
abline(h = limite_superiore_max, col = "red", lty = 2)

# Add a legend to the max norm plot
legend("topright", legend = c("Max norm", "Control limit"),
       col = c("green", "red"), lty = c(1, 2), pch = c(16, NA))

# Plot the total variation
plot(total_variation, type = "o", col = "purple", pch = 16,
     main = "Control Chart for Milano-Varese OD matrices - Total Variation",
     xlab = "Week", ylab = "Total variation",
     ylim = c(min(total_variation), max(total_variation, limite_superiore_total_variation)))
# Add control limit line for total variation
abline(h = limite_superiore_total_variation, col = "red", lty = 2)

# Add a legend to the total variation plot
legend("topright", legend = c("Total variation", "Control limit"),
       col = c("purple", "red"), lty = c(1, 2), pch = c(16, NA))

frobenius_norm <- function(mat, media) {
  sqrt(sum((mat - media)^2))
}

# Function to calculate the norm difference
norm_diff <- function(mat, media) {
  sqrt(sum((mat^2))) - sqrt(sum((media^2)))
}

# Initialize vectors to store results for adaptive control chart
window_size <- 3
norm_frobenius <- numeric(t - window_size + 1)
limiti_superiori1 <- numeric(t - window_size + 1)
limiti_superiori2 <- numeric(t - window_size + 1)
limiti_superiori3 <- numeric(t - window_size + 1)
limiti_inferiori1 <- numeric(t - window_size + 1)
limiti_inferiori2 <- numeric(t - window_size + 1)
limiti_inferiori3 <- numeric(t - window_size + 1)

# Loop through data to calculate norms and control limits
for (i in 1:(t - window_size + 1)) {
  # Extract the current window of data
  window_data <- OD_list[i:(i + window_size - 1)]
  
  # Calculate the mean of the current window
  media_window <- Reduce("+", window_data) / length(window_data)
  
  # Calculate the Frobenius norm difference for the current matrix
  norm_frobenius[i] <- norm_diff(OD_list[[i]], media_window)
  
  # Calculate norms within the current window for further analysis
  norms_window <- sapply(window_data, norm_diff, media = media_window)
  media_norm_window <- mean(norms_window)
  dev_std_norm_window <- sd(norms_window)
  
  # Calculate adaptive control limits for the current window
  limiti_superiori1[i] <- media_norm_window + 1 * dev_std_norm_window
  limiti_superiori2[i] <- media_norm_window + 2 * dev_std_norm_window
  limiti_superiori3[i] <- media_norm_window + 3 * dev_std_norm_window
  limiti_inferiori1[i] <- media_norm_window - dev_std_norm_window
  limiti_inferiori2[i] <- media_norm_window - 2 * dev_std_norm_window
  limiti_inferiori3[i] <- media_norm_window - 3 * dev_std_norm_window
}

# Plot the adaptive control chart for Frobenius norms
plot(norm_frobenius, type = "o", col = "blue", pch = 16,
     main = "Adaptive Control Chart for Milano-Varese OD matrices",
     xlab = "Week", ylab = "Frobenius norm",
     ylim = range(c(norm_frobenius, limiti_inferiori3, limiti_superiori3)))

# Add adaptive control limits to the plot
lines(limiti_superiori1, col = "green", lty = 2)
lines(limiti_superiori2, col = "orange", lty = 2)
lines(limiti_superiori3, col = "red", lty = 2)
lines(limiti_inferiori1, col = "green", lty = 2)
lines(limiti_inferiori2, col = "orange", lty = 2)
lines(limiti_inferiori3, col = "red", lty = 2)

# Add a legend to the adaptive control chart
legend("topright", legend = c("Frobenius norm", "Zone C", "Warning limits", "Control limits"),
       col = c("blue", "green", "orange", "red"), lty = c(1, 2, 3, 4), pch = c(8, NA))

# Parameters for the Exponentially Weighted Moving Average (EWMA) chart
lambda <- 0.2  # Smoothing parameter
threshold <- 3  # Threshold for detecting shifts

# Create EWMA chart for norm differences
ewma_chart <- ewma(norm_difference, lambda = lambda, center = mean(norm_difference), 
                   std.dev = sd(norm_difference), chart = TRUE)