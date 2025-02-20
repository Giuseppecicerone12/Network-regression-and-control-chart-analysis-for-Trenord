# Install and load necessary libraries
install.packages("lubridate")
library(lubridate)
library(SuperCENT)

# Import the dataset
raw_df = read.csv("Data/Raw/apc.csv")

# Convert the 'data' column to date format
raw_df$data <- dmy(raw_df$data)

# Convert time-related columns to date-time format
for (x in c("partenza", "arrivo", "ora_ingresso", "ora_uscita")) {
  raw_df[[x]] <- dmy_hms(raw_df[[x]])
}

# Process the removed stop flag
raw_df$fermata_soppressa <- raw_df$fermata_soppressa == "true"
# Remove records with negative departure delays
raw_df = raw_df[raw_df$ritardo_uscita >= 0,]
# Set negative arrival delays to zero
raw_df[raw_df$ritardo_ingresso < 0, "ritardo_ingresso"] = 0

# Plot the density of arrival delays
plot(density(raw_df$ritardo_ingresso))
# Calculate total delay and week number
raw_df$ritardo = raw_df$ritardo_uscita - raw_df$ritardo_ingresso 
raw_df$settimana = week(raw_df$data)
plot(density(raw_df$ritardo),xlim=c(-1000,1000))
# Calculate percentiles and mean delays
percentili <- quantile(raw_df$ritardo, probs = c(0.25, 0.5, 0.75))
print(percentili)
media = mean(raw_df$ritardo)
dev = sd(raw_df$ritardo)

# Calculate noise to add to the delays
dev_rumore = 0.01 * dev
set.seed(1)
rumore = rnorm(length(raw_df$ritardo), 0, dev_rumore)
raw_df$ritardo = raw_df$ritardo + rumore

# Load station codes
misc = load("Data/Processed/misc.Rdata")
station_codes <- station_codes_agg

# Aggregate data by stations and weeks
raw_df[raw_df$stazione %in% urban_area, "stazione"] <- "MI_AGG"
raw_df <- raw_df |> group_by(stazione, settimana) |> summarise(ritardo_totale = sum(ritardo))
plot(density(raw_df$ritardo_totale),xlim=c(-30000,30000))
# Create a list of aggregated data for the weeks
Y = list()
for (w in raw_df$settimana) {
  Y[[w]] = raw_df[raw_df$settimana == w, c("stazione", "ritardo_totale")]
}

# Import the OD (origin-destination) dataset
OD <- read.csv("Data/Processed/IPF/OD_Trenord_IPF_agg_LM_1_e.csv")
OD <- na.omit(OD)

# Process weeks for the OD dataset
weeks <- colnames(OD)[-c(1, 2)]
weeks <- substring(weeks, 6)
weeks <- as.Date(unlist(lapply(weeks, function(x) {
  if (str_sub(x, -2) == "00") {
    return(as.character(as.Date(paste(substr(x, 1, 4), 1, 1, sep = "-"))))
  }
  return(as.character(as.Date(paste(x, 1, sep = "_"), "%Y_%W_%w")))
})))

# Create a list of OD matrices for the weeks
A_list = list()
for (w in weeks) {
  A_list[[w]] = OD_week(OD, w, station_codes)
}

# Clean the data for the weeks
for (w in 1:length(weeks)) {
  Y[[w]] = Y[[w]][Y[[w]]$stazione %in% rownames(A_list[[weeks[1]]]),]
  Y[[w]] = Y[[w]][order(match(Y[[w]]$stazione, rownames(A_list[[weeks[1]]]))),]
}

# Prepare the design matrix for analyses
X <- rep(1, nrow(Y[[1]]))
X = matrix(X)
A = list()
for (w in weeks) {
  A[[w]] = as.matrix(A_list[[w]])
}

# Remove specific weeks and prepare for modeling
Y = Y[-c(52, 53)]
new_weeks = weeks[-c( 52, 53)]
lin_mod = list()
# Create a list to store results
ret = list()
new_weeks = weeks[-c(35, 36, 52, 53)]

# Initialize linear models
lin_mod = list()
for (w in 1:length(Y)) {
  # Transform the total delays using a negative log transformation
  Y[[w]]$ritardo_totale = neglogTrans$transform(Y[[w]]$ritardo_totale)
  # Plot the density of the transformed delays
  plot(density(Y[[w]]$ritardo_totale))
}

# List to store results for the SuperCENT model
sup=list()
for (w in 1:length(Y)) {
  ret[[w]] <- two_stage(A[[new_weeks[w]]], X, Y[[w]]$ritardo_totale, lrange = 2^4, gap = 2, folds = 3)
}

# Cross-validation for the SuperCENT model
for (w in 1:length(Y)) {
  ret[[w]] <- cv.supercent(A[[new_weeks[w]]], X, Y[[w]]$ritardo_totale, l = nrow(Y[[1]]) * (ret[[w]]$epsy)^2 / (ret[[w]]$epsa), lrange = 2^4, gap = 2, folds = 4)
}

# Second stage modeling
for (w in 1:length(Y)) {
  ret[[w]] <- two_stage(A[[new_weeks[w]]], X, Y[[w]]$ritardo_totale, lrange = 2^4, gap = 2, folds = 3)
}

# Re-run SuperCENT model
for (w in 1:length(Y)) {
  sup[[w]] <- supercent(A[[new_weeks[w]]], X, Y[[w]]$ritardo_totale, l = nrow(Y[[1]]) * (ret[[w]]$epsy)^2 / (ret[[w]]$epsa)^2)
}

# Save the results of the SuperCENT model
save(sup, file = "centralities_oracle.RData")

# Preparing for analysis of the effects
betas_u=list()
betas_v=list()
U = list()
V=list()
pvalues_u = list()
pvalues_v = list()
lambda = list()
diff_cent = list()
diff = list()

# Loop through each result to extract parameters
for (w in 1:length(sup)) {
  betas_u[[w]] = sup[[w]]$beta[2]
  betas_v[[w]] = sup[[w]]$beta[3]
  U[[w]] = sup[[w]]$u
  V[[w]] = sup[[w]]$v
  pvalues_u[[w]] = confint(sup[[w]])$p[2]
  pvalues_v[[w]] = confint(sup[[w]])$p[3]
  lambda[[w]] = sup[[w]]$l
  diff_cent[[w]] = sqrt(sum((U[[w]] - V[[w]])^2))  # Calculate the difference in centrality
  diff[[w]] = (U[[w]] - V[[w]])  # Difference between U and V
}

# Calculate errors and RMSE for model assessment
error = list()
rmse = list()
for (w in 1:(length(Y) - 1)) {
  matrice = cbind(X, U[[w]], V[[w]])
  error[[w + 1]] = Y[[w + 1]]$ritardo_totale - matrice %*% sup[[w]]$beta
  rmse[[w + 1]] = sqrt(mean(error[[w + 1]]^2))  # Compute root mean square error
}

# Save the linear model results
save(lin_mod, file = "linmod_neglog.RData")

# Prepare data for plotting centralities
st_names = rownames(A_list[[weeks[1]]])
names = load_stations_names()
names = names[names$Code %in% st_names,]
names = names[match(st_names, names$Code),]

# Create a data frame for SVD results
svd_data <- data.frame(
  "Station" = names$Name,
  "Hub" = U[[28]],
  "Authority" = V[[28]]
)
svd_data = svd_data[-1,]  # Remove the first row if not needed

# Load ggplot2 and ggrepel for plotting
library(ggplot2)
library(ggrepel)
# Plot hub and authority centrality
ggplot(svd_data, aes(x = Hub, y = Authority, label = Station)) +
  geom_point(color = "blue") +
  geom_text_repel(color = "red") + 
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") + 
  labs(title = "Values of Hub and Authority Centrality in Week 28") +
  theme_minimal()

# Create a data frame for hub centrality across stations and weeks
hub_data <- expand.grid(Week = weeks, Station = colnames(A[[1]])) |>
  add_column(hub_centrality = NA)

# Fill the hub centrality data frame
for (i in 1:length(colnames(A[[1]]))) {
  for (w in 1:length(ret)) {
    hub_data[hub_data$Week == weeks[w] & hub_data$Station == colnames(A[[1]])[i], "hub_centrality"] = U[[w]][i]
  }
}

hub_data = na.omit(hub_data)  # Remove NA values
hb_df <- hub_data |> 
  dplyr::select(Week, Station, hub_centrality) |> 
  pivot_wider(names_from = Week, values_from = hub_centrality)

# Rename columns for clarity
w_tot <- length(colnames(hb_df))
colnames(hb_df)[2:w_tot] <- paste0("Week_", colnames(hb_df)[2:w_tot])
row.names(hb_df) <- hb_df$Station

# Transpose the data frame for spline analysis
x <- t(hb_df[, 2:w_tot])
t <- rep(1:(w_tot - 1))
degree <- 3  # Spline degree
m <- degree + 1  # Spline order

# Select number of basis functions and compute GCV
nb <- as.data.frame(hb_df[, 1])
s_tot <- length(station_codes)
x = na.omit(x)

# Initialize basis count
for (i in 1:s_tot) {
  nb[i, 2] <- 0
}

# Compute GCV for various basis counts
for (j in 1:s_tot) {
  x1 <- x[, j]
  nbasis <- 4:20
  gcv <- numeric(length(nbasis))
  
  for (i in 1:length(nbasis)) {
    basis <- create.bspline.basis(rangeval = c(1, w_tot - 1), nbasis[i], m)
    gcv[i] <- smooth.basis(t, x1, basis)$gcv
  }
  
  # Plot GCV results
  par(mfrow = c(1, 1))
  plot(nbasis, gcv)
  nb[j, 2] <- nbasis[which.min(gcv)]  # Store the optimal basis count
}

# Print the mean optimal number of basis functions
print(mean(nb[, 2]))

# Round to select the number of basis functions
nbasis <- round(mean(nb[, 2]))

# Create B-spline basis
basis <- create.bspline.basis(rangeval = c(1, (w_tot - 1)), nbasis = nbasis, norder = m)
Xsp <- smooth.basis(argvals = t, y = x, fdParobj = basis)
Xsp0bis <- eval.fd(t, Xsp$fd)

# Prepare the data frame for the smoothed hub centrality
hb_func_df <- as.data.frame(Xsp0bis)
colnames(hb_func_df) <- hb_df$Station
hb_func_df$Week <- rownames(x)

# Format the Week column
hb_func_df$Week <- as.Date(unlist(lapply(hb_func_df$Week, function(x) {
  if (str_sub(x, -2) == "00") {
    return(as.character(as.Date(paste(substr(x, 1, 4), 1, 1, sep = "-"))))
  }
  return(as.character(as.Date(paste(x, 1, sep = "_"), "%Y_%W_%w")))
})))

# Convert to long format for visualization
hb_func_df <- hb_func_df |> pivot_longer(cols = colnames(hb_func_df)[1:s_tot], names_to = "Station", values_to = "y")

# Create functional data object for hub centrality
fd_hub <- Data2fd(y = t(as.matrix(hb_df[-1, -1])), argvals = t, basisobj = basis)
plot.fd(fd_hub, titles = colnames(A[[1]]), xlab = "Week", ylab = "Hub Centralities")
lines(mean.fd(fd_hub), lwd = 3)  #
# Highlight specific curves to label
curve_indices <- c(8, 18, 25, 31, 33, 34, 35, 36)  # Indices of curves to label
curve_labels <- c("Bergamo", "Verdello Daimine", "Treviglio", "Rovato", "Brescia", "Desenzano Del Garda-Sirmione", "Peschiera Del Garda", "Verona Porta Nuova")  # Corresponding labels

# Extract data to position the labels
x_vals <- seq(min(fd_hub$basis$rangeval), max(fd_hub$basis$rangeval), length.out = 100)  # X interval
y_vals <- eval.fd(x_vals, fd_hub)[, curve_indices]  # Values of the curves to label

# Position the labels near the maximum or average values of the curves
for (i in seq_along(curve_indices)) {
  text(x_vals[length(x_vals) / 2], y_vals[nrow(y_vals) / 2, i], labels = curve_labels[i], pos = 4, col = i)
}

# Create functional data object for Milan hub centrality
fd_hub_Mi <- Data2fd(y = t(as.matrix(hb_df[1, -1])), argvals = t, basisobj = basis)
plot.fd(fd_hub_Mi, xlab = "Week", ylab = "Milan Hub Centrality Over Weeks")

# Calculate authority centralities
aut_data <- expand.grid(Week = weeks, Station = colnames(A[[1]])) |>
  add_column(aut_centrality = NA)

# Fill the authority centrality data frame
for (i in 1:length(colnames(A[[1]]))) {
  for (w in 1:length(ret)) {
    aut_data[aut_data$Week == weeks[w] & aut_data$Station == colnames(A[[1]])[i], "aut_centrality"] = V[[w]][i]
  }
}

aut_data = na.omit(aut_data)  # Remove NA values
aut_df <- aut_data |>
  dplyr::select(Week, Station, aut_centrality) |>
  pivot_wider(names_from = Week, values_from = aut_centrality)

# Rename columns for clarity
w_tot <- length(colnames(aut_df))
colnames(aut_df)[2:w_tot] <- paste0("Week_", colnames(aut_df)[2:w_tot])
row.names(aut_df) <- aut_df$Station

# Transpose the data frame for spline analysis
x <- t(aut_df[, 2:w_tot])
t <- rep(1:(w_tot - 1))
degree <- 3  # Spline degree
m <- degree + 1  # Spline order

# Select number of basis functions and compute GCV for authority centrality
nb <- as.data.frame(aut_df[, 1])
s_tot <- length(station_codes)
x = na.omit(x)

# Initialize basis count for authority centrality
for (i in 1:s_tot) {
  nb[i, 2] <- 0
}

# Compute GCV for various basis counts
for (j in 1:s_tot) {
  x1 <- x[, j]
  nbasis <- 4:20
  gcv <- numeric(length(nbasis))
  
  for (i in 1:length(nbasis)) {
    basis <- create.bspline.basis(rangeval = c(1, w_tot - 1), nbasis[i], m)
    gcv[i] <- smooth.basis(t, x1, basis)$gcv
  }
  
  # Plot GCV results
  par(mfrow = c(1, 1))
  plot(nbasis, gcv)
  nb[j, 2] <- nbasis[which.min(gcv)]  # Store the optimal basis count
}

# Print the mean optimal number of basis functions for authority centrality
print(mean(nb[, 2]))

# Round to select the number of basis functions
nbasis <- round(mean(nb[, 2]))

# Create B-spline basis for authority centrality
basis <- create.bspline.basis(rangeval = c(1, (w_tot - 1)), nbasis = nbasis, norder = m)
Xsp <- smooth.basis(argvals = t, y = x, fdParobj = basis)
Xsp0bis <- eval.fd(t, Xsp$fd)

# Prepare the data frame for the smoothed authority centrality
aut_func_df <- as.data.frame(Xsp0bis)
colnames(aut_func_df) <- aut_df$Station
aut_func_df$Week <- rownames(x)
# Format the Week column for authority centrality
aut_func_df$Week <- as.Date(unlist(lapply(aut_func_df$Week, function(x) {
  if (str_sub(x, -2) == "00") {
    return(as.character(as.Date(paste(substr(x, 1, 4), 1, 1, sep = "-"))))
  }
  return(as.character(as.Date(paste(x, 1, sep = "_"), "%Y_%W_%w")))
})))

# Convert to long format for visualization of authority centrality
aut_func_df <- aut_func_df |> pivot_longer(cols = colnames(aut_func_df)[1:s_tot], names_to = "Station", values_to = "y")

# Create functional data object for authority centrality
fd_aut <- Data2fd(y = t(as.matrix(aut_df[-1, -1])), argvals = t, basisobj = basis)

# Plot authority centralities
plot.fd(fd_aut, titles = colnames(A[[1]]), xlab = "Week", ylab = "Authority Centralities")
lines(mean.fd(fd_aut), lwd = 3)  # Add a line for the mean authority centrality

# Highlight specific curves to label for authority centrality
curve_indices <- c(8, 18, 25, 31, 33, 34, 35, 36)  # Indices of curves to label
curve_labels <- c("Bergamo", "Verdello Daimine", "Treviglio", "Rovato", "Brescia", "Desenzano Del Garda-Sirmione", "Peschiera Del Garda", "Verona Porta Nuova")  # Corresponding labels

# Extract data to position the labels for authority centrality
x_vals <- seq(min(fd_aut$basis$rangeval), max(fd_aut$basis$rangeval), length.out = 100)  # X interval
y_vals <- eval.fd(x_vals, fd_aut)[, curve_indices]  # Values of the curves to label

# Position the labels near the maximum or average values of the curves for authority centrality
for (i in seq_along(curve_indices)) {
  text(x_vals[length(x_vals) / 2], y_vals[nrow(y_vals) / 2, i], labels = curve_labels[i], pos = 4, col = i)
}

# Create functional data object for Milan authority centrality
fd_aut_Mi <- Data2fd(y = t(as.matrix(aut_df[1, -1])), argvals = t, basisobj = basis)
plot.fd(fd_aut_Mi, xlab = "Week", ylab = "Milan Authority Centrality Over Weeks")




