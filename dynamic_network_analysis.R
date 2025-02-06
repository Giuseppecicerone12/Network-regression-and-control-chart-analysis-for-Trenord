## Analyze the Trenord dynamic OD matrices using network analysis and FDA techniques (as per https://github.com/GretaGalliani/dynamic-OD-estimation-railway-network)
## Prerequisities: run pipeline (https://github.com/adri-pi/od-estimation) 
## Input:
##  strategy: gapfilling strategy identifier (defaults to LM_1_e)
##  aggregateUrban: if we should aggregate urban area (defaults to TRUE)
##  subdir: subdirectory for specific filters, defaults to "" (no filter)
##  savePlots: if we should save plots instead of showing them (defaults to TRUE)
## Output: plot is saved in folder as per line above
analyze <- function(strategy = "LM_1_e", aggregateUrban = TRUE, subdir = "", savePlots = TRUE) {
  ## Load stations codes
  load("Data/Processed/misc.Rdata")
  if (aggregateUrban) {
    station_codes <- station_codes_agg
    strategy <- paste0("agg_", strategy)
  }
  
  print(station_codes)
  ## Load network
  OD <- read.csv(paste0("Data/Processed/IPF/", subdir, "/OD_Trenord_IPF_", strategy, ".csv"))
  daytag <- read_excel("Data/daytag.xlsx")
  daytag <- daytag[ format(daytag$Data, "%Y") == "2023",]
  daytag <- daytag |> mutate(Week = as.character(format(as.Date(Data), "%Y_%W")))
  daytag <- daytag |> group_by(Week) |> summarise(Festivity= sum(Festivity),Holiday = sum(Holiday),Strike=sum(Strike))
  daytag$event = rep(NA,nrow(daytag))
  daytag$event <- apply(daytag[, c("Festivity", "Holiday", "Strike")], 1, function(x) {
    colnames(daytag[,c("Festivity", "Holiday", "Strike")])[which.max(x)]
  })
  daytag[daytag$Festivity+daytag$Holiday+daytag$Strike == 0,"event"] = "Normal week"
  daytag$Week <- as.Date(unlist(lapply(daytag$Week, function(x) {
    if (str_sub(x, -2) == "00") {
      return(as.character(as.Date(paste(substr(x, 1, 4), 1, 1, sep = "-"))))
    }
    return(as.character(as.Date(paste(x, 1, sep = "_"), "%Y_%W_%w")))
  })))
  daytag$Week <- as.Date(daytag$Week)
  #names(OD) <- make.names(names(OD)) 
  #colnames(OD)[1] <- "Start"
  OD <- na.omit(OD)
  if (!dir.exists(paste0("Data/Processed/Dynamic_network_analysis/", subdir))) dir.create(paste0("Data/Processed/Dynamic_network_analysis/", subdir), recursive = TRUE)

  ####  1. GLOBAL METRICS ----
  #### 1.1. Global comparison between matrices through RMSE-----
  weeks <- colnames(OD)[-c(1, 2)]
  weeks <- substring(weeks, 6)
  #print(weeks)
  ## Preparing dataset
  Diff_norm <- data.frame(
    "Week" = weeks,
    "RMSE" = rep(NA, length(weeks))
  )

  for (i in 2:length(weeks)) {
    ## Select needed weeks for the comparison
    w_prev <- weeks[i - 1]
    w_cur <- weeks[i]
    
    ## Selecting the OD matrices for the two needed weeks in the whole dataset
    OD_prev <- OD_week(OD, w_prev, station_codes)
    OD_cur <- OD_week(OD, w_cur, station_codes)
    ## Computing the RMSE
    Diff_norm[Diff_norm$Week == w_cur, "RMSE"] <- sqrt(sum((OD_cur - OD_prev)^2) / length(OD_cur))
  }
  #print(OD_prev)
  #print(OD_cur)
  Diff_norm$Week <- as.Date(unlist(lapply(Diff_norm$Week, function(x) {
    if (str_sub(x, -2) == "00") {
      return(as.character(as.Date(paste(substr(x, 1, 4), 1, 1, sep = "-"))))
    }
    return(as.character(as.Date(paste(x, 1, sep = "_"), "%Y_%W_%w")))
  })))
Diff_norm <- left_join(Diff_norm, daytag, by = "Week")
Diff_norm <- na.omit(Diff_norm)
event_colors <- c("Festivity" = "green", "Holiday" = "blue", "Strikes" = "red","Normal_week" = "black")
ggplot(Diff_norm, aes(x = Week, y = RMSE, group = 1)) +
  geom_line(color = "black", size = 1) +  # Linea nera
  geom_point(aes(color = Event_type), size = 4) +  # Punti colorati
  scale_color_manual(values = event_colors) +  # Applica i colori definiti
  labs(title = "Estimation of dynamic OD matrices in a railway transportation network",
       x = "2023",
       y = "RMSE",
       color = "Events' type") +  # Titoli degli assi e legenda
  theme_minimal() +  # Tema minimale
  theme(plot.title = element_text(hjust = 0.5))
  ## Saving the result
  write.csv(Diff_norm, paste0("Data/Processed/Dynamic_network_analysis/", subdir, "global_MSE_", strategy, ".csv"), row.names = FALSE)

  #### 1.2. Mean Strength -----
  strength <- expand.grid(Week = weeks, Station = station_codes) |>
    add_column(Strength_in = NA, Strength_out = NA, Strength_all = NA)

  for (w in weeks) {
    ## Get the OD matrix
    OD_cur <- OD_week(OD, w, station_codes)

    ## Convert to a weighted directed network object
    net <- graph_from_adjacency_matrix(as.matrix(OD_cur), mode = "directed", weighted = TRUE)
    
    ## Get local strength in/out/all
    strength[strength$Week == w, "Strength_all"] <- strength(net, vids = station_codes, mode = "all")
    strength[strength$Week == w, "Strength_in"] <- strength(net, vids = station_codes, mode = "in")
    strength[strength$Week == w, "Strength_out"] <- strength(net, vids = station_codes, mode = "out")
  }

  ## Fixing date format
  strength$Week <- as.Date(unlist(lapply(strength$Week, function(x) {
    if (str_sub(x, -2) == "00") {
      return(as.character(as.Date(paste(substr(x, 1, 4), 1, 1, sep = "-"))))
    }
    return(as.character(as.Date(paste(x, 1, sep = "_"), "%Y_%W_%w")))
  })))

  strength <- strength |>
    group_by(Station) |>
    mutate(
      Total_all = sum(Strength_all),
      Total_in = sum(Strength_in),
      Total_out = sum(Strength_out)
    ) |>
    ungroup() |>
    ## Divide each Strenghts value by the corresponding Total
    mutate(
      Strength_all_norm = Strength_all / Total_all,
      Strength_in_norm = Strength_in / Total_in,
      Strength_out_norm = Strength_out / Total_out
    ) |>
    dplyr::select(-c(Total_all, Total_in, Total_out))
  write.csv(strength, paste0("Data/Processed/Dynamic_network_analysis/", subdir, "strength_", strategy, ".csv"), row.names = FALSE)

  ## I compute the mean strenght in the network
  mean_strength <- strength |>
    dplyr::select(-c(Strength_all_norm, Strength_in_norm, Strength_out_norm)) |>
    group_by(Week) |>
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
  mean_strength<- left_join(mean_strength, daytag, by = "Week")
  mean_strength <- na.omit(mean_strength)
  ## Saving the result
  write.csv(mean_strength, paste0("Data/Processed/Dynamic_network_analysis/", subdir, "mean_strength_", strategy, ".csv"), row.names = FALSE)


  #### 2. LOCAL METRICS  ----
  ## Strength can also be viewed as a local metric

  #### 2.1. FDA ON STRENGTH VALUES ----
  ### GOAL: Analysis on total strength using FDA techniques
  ## I consider strength, normalized dividing by the total strength for each station through the period considered

  #### 2.1.1. FDA - SMOOTHING ----
  ## Transform dataset in the format needed to apply FDA
  st_df <- strength |>
    dplyr::select(Week, Station, Strength_all_norm) |>
    pivot_wider(names_from = Week, values_from = Strength_all_norm)
  w_tot <- length(colnames(st_df))
  colnames(st_df)[2:w_tot] <- paste0("Week_", colnames(st_df)[2:w_tot])

  ## SMOOTHING - Cubic splines
  row.names(st_df) <- st_df$Station
  x <- t(st_df[, 2:w_tot])
  t <- rep(1:(w_tot - 1))
  degree <- 3 ## spline degree
  m <- degree + 1 ## spline order

  ## Select nbasis -> compute gcv
  nb <- as.data.frame(st_df[, 1])
  s_tot <- length(station_codes)
  for (i in 1:s_tot) {
    nb[i, 2] <- 0
  }

  for (j in 1:s_tot) {
    x1 <- x[, j]
    nbasis <- 4:20
    gcv <- numeric(length(nbasis))
    for (i in 1:length(nbasis)) {
      basis <- create.bspline.basis(rangeval = c(1, w_tot - 1), nbasis[i], m)
      gcv[i] <- smooth.basis(t, x1, basis)$gcv
    }
    par(mfrow = c(1, 1))
    plot(nbasis, gcv)
    nb[j, 2] <- nbasis[which.min(gcv)]
  }

  print(mean(nb[, 2]))
  ## I round this quantity and select nbasis
  nbasis <- round(mean(nb[, 2]))

  ## Create b-spline basis
  basis <- create.bspline.basis(rangeval = c(1, (w_tot - 1)), nbasis = nbasis, norder = m)
  Xsp <- smooth.basis(argvals = t, y = x, fdParobj = basis)
  Xsp0bis <- eval.fd(t, Xsp$fd)

  ## Adjusting functional data
  st_func_df <- as.data.frame(Xsp0bis)
  colnames(st_func_df) <- st_df$Station
  st_func_df$Week <- weeks
  ## Fixing date format
  st_func_df$Week <- as.Date(unlist(lapply(st_func_df$Week, function(x) {
    if (str_sub(x, -2) == "00") {
      return(as.character(as.Date(paste(substr(x, 1, 4), 1, 1, sep = "-"))))
    }
    return(as.character(as.Date(paste(x, 1, sep = "_"), "%Y_%W_%w")))
  })))
  st_func_df <- st_func_df |> pivot_longer(cols = colnames(st_func_df)[1:s_tot], names_to = "Station", values_to = "y")
  fd <- Data2fd(argvals = x, y = t, basisobj = basis)
  pca_results <- pca.fd(fd,nharm = 5)
  pca <- plot.pca.fd(pca_results)
  print(pca_results$varprop)
  plot(pca_results$scores[, 1], pca_results$scores[, 2],
       xlab = "PC1", ylab = "PC2", 
       main = "Scores along PC1 e PC2",
       pch = 19, col = "blue")
  
  # Aggiungere i nomi delle osservazioni
  text(pca_results$scores[, 1], pca_results$scores[, 2], 
       labels = st_names$Name, 
       pos = 4, # Posizione del testo (4 = a destra del punto)
       cex = 0.2, # Dimensione del testo
       col = "red") # Colore del testo
  pc_data <- data.frame(
    "Station" = st_names$Name,
    "PC1" = pca_results$scores[, 1],
    "PC2" = pca_results$scores[, 2]
  )
  library(ggplot2)
  library(ggrepel)
  ggplot(pc_data, aes(x = PC1, y = PC2, label = Station)) +
    geom_point(color = "blue") +
    geom_text_repel(color = "red") +  # Usa geom_text_repel per evitare sovrapposizioni
    labs(title = "Score along PC1 e PC2") +
    theme_minimal()
  ## Saving the result
  write.csv(st_func_df, paste0("Data/Processed/Dynamic_network_analysis/", subdir, "local_strength_function_", strategy, ".csv"), row.names = FALSE)
  if (savePlots) png(filename = paste0("Data/Processed/Dynamic_network_analysis/", subdir, "FPCA_", strategy, ".png"), width = 800, height = 600)
  ## Plot the functional mean
  if (savePlots) png(filename = paste0("Data/Processed/Dynamic_network_analysis/", subdir, "functional_mean_", strategy, ".png"), width = 800, height = 600)
  fd_strength <- Data2fd(y = t(as.matrix(st_df[, -1])), argvals = t, basisobj = basis)
  plot.fd(fd_strength, titles = st_df$Station)
  lines(mean.fd(fd_strength), lwd = 3)
  if (savePlots) dev.off()

  ## Plot the functional covariance
  eval <- eval.fd(t, fd_strength)
  if (savePlots) png(filename = paste0("Data/Processed/Dynamic_network_analysis/", subdir, "functional_covariance_", strategy, ".png"), width = 800, height = 600)
  image.plot(t, t, (cov(t(eval))[1:(w_tot - 1), ]))
  if (savePlots) dev.off()
  ## It makes sense: weeks near each other are more similar than weeks far in time

  #### 2.1.2. FDA - OUTLIERS ----
  ## Consider data
  grid <- seq(1, (w_tot - 1))
  fData <- fData(grid, t(eval))
  plot(fData)
  if (savePlots) png(filename = paste0("Data/Processed/Dynamic_network_analysis/", subdir, "mobility_data_smoothed_", strategy, ".png"), width = 800, height = 600)
  matplot(Xsp0bis, lty = 1, type = "l", main = "Mobility data smoothed - Provinces", ylab = "Densities of strength", xlab = "Weeks")
  if (savePlots) dev.off()

  ## a. Functional boxplot
  if (savePlots) png(filename = paste0("Data/Processed/Dynamic_network_analysis/", subdir, "functional_boxplot_", strategy, ".png"), width = 800, height = 600)
  fb <- fbplot(fData, xlab = "Week", ylab = "Mobility densities", main = "Functional Boxplot")
  fb
  if (savePlots) dev.off()
  ## b. Outliergram
  if (savePlots) png(filename = paste0("Data/Processed/Dynamic_network_analysis/", subdir, "outliergram_", strategy, ".png"), width = 800, height = 600)
  pr <- outliergram(fData)
  if (savePlots) dev.off()

  ## The outliers are
  outliers <- st_df$Station[fb$ID_outliers]
  print(outliers)
  names = load_stations_names()
  #names = names[names$Code %in% station_codes,"Name"]
  names = names[names$Code %in% outliers,"Name" ]
  print(names)
  outliers <- as.character(st_df$Station[fb$ID_outliers])
  station_codes <- as.character(station_codes)
  palette_colors <- rainbow(length(unique(station_codes[station_codes %in% outliers])))
  my_colors <- rep("black", length(station_codes))
  
  print(palette_colors)
  k=1
  for(i in 1:length(station_codes)){
    if(station_codes[i] %in% outliers){
      my_colors[i] = palette_colors[k]
      k<-k+1
    }
  }
  #print(my_colors)   
  #print(k)    
  print(length(my_colors))
  
  
  print(ncol(Xsp0bis))
    
  fb_color <- matplot(Xsp0bis, lty = 1, type = "l", main = "Mobility data smoothed - Stations", ylab = "Densities of strength", xlab = "Weeks",col=my_colors)
  legend("topright", legend = names, col = palette_colors, pch = 15, title = "Outliers",cex = 0.7, box.lwd = 0.05)
  fb_color
  ## Outliers interpreted by station
  st_func_df$out_group <- "Normal"
  st_func_df[st_func_df$Station %in% outliers, "out_group"] <- "Outlier"
  
  write.csv(st_func_df, paste0("Data/Processed/Dynamic_network_analysis/", subdir, "functional_outliers_", strategy, ".csv"), row.names = F)
}
