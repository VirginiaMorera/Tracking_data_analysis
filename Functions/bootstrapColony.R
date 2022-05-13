# Developed for the manuscript *Detecting recurrent sources of variability in animal tracking studies* by Virginia Morera-Pujol *et al.* submitted to Ecological Applications in 2021.

### bootstrapColony function ####
# tracks: Tracking data. DF or SPDF, containing Latitude, Longitude, Species and Population
# simulatedDistributions: DF or SPDF containing the simulated distribution of the studied 
#   species, with at least the columns Longitude, Latitude, Population and Species. Can be directly the output of the simulateDistribution() function
# tripID: unique identifier for each trip 
# udLev: Numeric, Utilisation Distribution Level. Level of the KDE to be selected. For example, if 50 is used, the overlap is calculated between the 
#   smallest areas for which the probability to find the animal is equal to 0.50 (but see explanation for argument "conditional". From the adehabitatHR::kerneloverlap() function.
# scale is the smoothing factor to be used in the Kernel Density Estimation (in Km)
# grid is a number giving the size of the grid on which the UD should be estimated. 
# iterations:  Numeric, number of iterations for the bootstrap

bootstrapColony <- function(tracks, simulatedDistributions, tripId, udLev = 50, scale, grid = 500, iterations = 50){
   
   # package dependencies
   if (!requireNamespace("sp", quietly = TRUE)) {
    stop("Package \"sp\" needed for  function to work. Please install.",
         call. = FALSE)  }
   if (!requireNamespace("geosphere", quietly = TRUE)) {
    stop("Package \"geosphere\" needed for  function to work. Please install.",
         call. = FALSE)  }
   if (!requireNamespace("rgdal", quietly = TRUE)) {
    stop("Package \"rgdal\" needed for  function to work. Please install.",
         call. = FALSE)  }
   if (!requireNamespace("adehabitatHR", quietly = TRUE)) {
      stop("Package \"adehabitatHR\" needed for  function to work. Please install.",
           call. = FALSE)  }  require(foreach)
   if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("Package \"doParallel\" needed for  function to work. Please install.",
         call. = FALSE)  }
   if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package \"parallel\" needed for  function to work. Please install.",
         call. = FALSE)  }
   if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("Package \"foreach\" needed for  function to work. Please install.",
         call. = FALSE)  }
  
  # initial checks
  tracks$tripId <- tracks[,tripId]
  if (!"Latitude" %in% names(tracks)) stop("Latitude field does not exist")
  if (!"Longitude" %in% names(tracks)) stop("Longitude field does not exist")
  if (!"tripId" %in% names(tracks)) stop("tripId field does not exist")
  if (!"Population" %in% names(tracks)) stop("Population field does not exist")
  if (!"Species" %in% names(tracks)) stop("Sp field does not exist")
  
  # select species corresponding to the Tracking data
  simulatedDistributions <- simulatedDistributions[simulatedDistributions$Species == unique(tracks$Species),]
  
  # Converts simulatedDistributions to spatial dataframe and projects it (and if already is SPDF check projection)
  if (class(simulatedDistributions) != "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
  {
    ## filter DF to the minimum fields that are needed
    CleansimulatedDistributions <- simulatedDistributions %>%
      dplyr::select(Species, Population, Longitude, Latitude)
    mid_point <- data.frame(centroid(cbind(CleansimulatedDistributions$Longitude, CleansimulatedDistributions$Latitude)))
    
    ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
    if (min(CleansimulatedDistributions$Longitude) < -170 &  max(CleansimulatedDistributions$Longitude) > 170) {
      longs = ifelse(CleansimulatedDistributions$Longitude < 0,CleansimulatedDistributions$Longitude + 360,CleansimulatedDistributions$Longitude)
      mid_point$lon <- ifelse(median(longs) > 180,median(longs) - 360,median(longs))}
    
    simulatedDistributions.Wgs <- SpatialPoints(data.frame(CleansimulatedDistributions$Longitude, CleansimulatedDistributions$Latitude), proj4string = CRS("+proj=longlat + datum=wgs84"))
    proj.Sim <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep = ""))
    simulatedDistributions.Projected <- spTransform(simulatedDistributions.Wgs, CRS = proj.Sim)
    simulatedDistributionsSpatial <- SpatialPointsDataFrame(simulatedDistributions.Projected, data = CleansimulatedDistributions)
    simulatedDistributionsSpatial@data <- simulatedDistributionsSpatial@data %>% dplyr::select(Species, Population)
    simulatedDistributions.Wgs <- NULL
    simulatedDistributions.Projected <- NULL
    
  } else {## if data are already in a SpatialPointsDataFrame then check for projection 
    if (is.na(simulatedDistributions@proj4string)) stop("The proj4string slot of simulatedDistributions can't be empty. Check and assign projection and re-run")
    if (!is.projected(simulatedDistributions)) { # if it's not projected (lonlat) project to laea around midpoint
      mid_point <- data.frame(centroid(cbind(simulatedDistributions@coords[,1], simulatedDistributions@coords[,2])))
      proj.Sim <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep = ""))
      ### MB  This part prevents projection problems around the DATELINE 
      if (min(simulatedDistributions@coords[,1]) < -170 &  max(simulatedDistributions@coords[1,]) > 170) {
        longs = ifelse(simulatedDistributions@coords[,1] < 0, simulatedDistributions@coords[,1] + 360, simulatedDistributions@coords[,1])
        mid_point$lon <- ifelse(median(longs) > 180, median(longs) - 360, median(longs))}
      simulatedDistributionsSpatial <- spTransform(simulatedDistributions, proj.Sim)
      simulatedDistributionsSpatial@data <- simulatedDistributionsSpatial@data %>% dplyr::select(Species, Population)
    } else {  ## if projected, "unproject" to WGS to project laea around midpoint
      simulatedDistributions.WGS <- spTransform(simulatedDistributions, CRS("+proj=longlat + datum=wgs84"))
      mid_point <- data.frame(centroid(cbind(simulatedDistributions.WGS@coords[,1], simulatedDistributions.WGS@coords[,2])))
      
      ### MB  This part prevents projection problems around the DATELINE 
      if (min(simulatedDistributions.WGS@coords[,1]) < -170 &  max(simulatedDistributions.WGS@coords[1,]) > 170) {
        longs = ifelse(simulatedDistributions@coords[,1] < 0, simulatedDistributions@coords[,1] + 360, simulatedDistributions@coords[,1])
        mid_point$lon <- ifelse(median(longs) > 180, median(longs) - 360, median(longs))}
      
      proj.Sim <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep = ""))
      simulatedDistributionsSpatial <- spTransform(simulatedDistributions.WGS, CRS = proj.Sim)
      simulatedDistributionsSpatial@data <- simulatedDistributionsSpatial@data %>% dplyr::select(Species, Population)
    }
  }
  
  # Converts tracks to spatial dataframe and projects it (and if already is SPDF check projection)
  if (class(tracks) != "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
  {
    ## filter DF to the minimum fields that are needed
    CleanTracks <- tracks %>%
      dplyr::select(Species, Population, tripId, Latitude, Longitude)
    
    Tracks.Wgs <- SpatialPoints(data.frame(CleanTracks$Longitude, CleanTracks$Latitude), proj4string = CRS("+proj=longlat + datum=wgs84"))
    Tracks.Projected <- spTransform(Tracks.Wgs, CRS = simulatedDistributionsSpatial@proj4string)
    TracksSpatial <- SpatialPointsDataFrame(Tracks.Projected, data = CleanTracks)
    TracksSpatial@data <- TracksSpatial@data %>% dplyr::select(Species, Population, tripId)
    Tracks.Wgs <- NULL
    Tracks.Projected <- NULL
    
  } else {## if data are already in a SpatialPointsDataFrame then check for projection
    if (is.na(tracks@proj4string)) stop("The proj4string slot of tracks can't be empty. Check and assign projection and re-run")
    TracksSpatial <- spTransform(tracks, CRS = simulatedDistributionsSpatial@proj4string)
    TracksSpatial@data <- TracksSpatial@data %>% dplyr::select(Species, Population, tripId)
  }
    
  
  # remove tripId tracks with < 6 points as they can't be used to calculate kernel
  UIDs <- names(which(table(TracksSpatial@data[, tripId]) > 6))             
  TracksSpatial <- TracksSpatial[TracksSpatial@data[, tripId] %in% UIDs, ]
  TracksSpatial@data[ ,tripId] <- droplevels(as.factor(TracksSpatial@data[ ,tripId]))
  
  TracksSpatial$X <- TracksSpatial@coords[,1]
  TracksSpatial$Y <- TracksSpatial@coords[,2]
  
  UIDs <- as.character(unique(TracksSpatial$tripId))
  Ntrips <- length(UIDs)
  Nloop <- seq(1,(Ntrips - 1),ifelse(Ntrips > 100,10,1))
  DoubleLoop <- data.frame(SampleSize = rep(Nloop,each = iterations), Iteration = rep(seq(1:iterations),length(Nloop)))
  LoopNr <- seq(1:dim(DoubleLoop)[1])	
  # udLev <- UDLev
  
  #setup parallel backend to use 4 processors
  cl <- makeCluster(detectCores(), outfile = "")
  registerDoParallel(cl)
  Result <- data.frame()
  
  Result <- foreach(LoopN = LoopNr, .combine = rbind, .packages = c("sp","adehabitatHR","geosphere","rgdal")) %dopar% {
    # output_list <- list()
    # for (j in seq_along(LoopNr))  {
    # j = 20
    # LoopN <- LoopNr[j]
    N <- DoubleLoop$SampleSize[LoopN]
    i <- DoubleLoop$Iteration[LoopN]
    # Coverage <- NULL
    # Inclusion <- NULL
    # History <- NULL
    
    Output <- data.frame(SampleSize = N, InclusionMean = 0,Iteration = i)
    # set.seed(123)
    RanNum <- sample(UIDs, N, replace = F)
    sink("D:/selected_ids_colony_bootstrap.txt", append = TRUE)
    cat(RanNum, "\n", "\n")
    sink()
    SelectedCoords <- TracksSpatial[TracksSpatial@data[,tripId] %in% RanNum,]
    # Ext <- (min(SelectedCoords@coords[,1]) + 3 * diff(range(SelectedCoords@coords[,1])))
    # if(Ext < (scale * 1000 * 2)) {
    # BExt <- ceiling((scale * 1000 * 3)/(diff(range(SelectedCoords@coords[,1])))) #} else {BExt <- 3}
    KDE.Surface <- kernelUD(SelectedCoords, h = scale*1000, grid = 500,  same4all = FALSE)
    try(KDE.UD <- getverticeshr(KDE.Surface, percent = udLev))
    if (isTRUE(class(KDE.UD) == "try-error")) {
      sink("errors.txt", append = T)
      cat(paste("Failed in iteration", i, "with sample", RanNum, sep = " "))
      sink()} else {
        KDE.UD@proj4string <- simulatedDistributionsSpatial@proj4string
        Overlain <- over(simulatedDistributionsSpatial, KDE.UD)$area
        Output$InclusionMean <- length(Overlain[!is.na(Overlain)])/nrow(simulatedDistributionsSpatial@data)
        return(Output)
        # output_list[j] <- list(Output)
      }
  }
    
  ## stop the cluster
  stopCluster(cl)
  closeAllConnections() 
  
  par(mfrow = c(1,1), mai = c(1,1,1,1))
  #Result <- Output[1:nrow(Output) - 1,]
  Result$Population <- unique(TracksSpatial$Population)
  M1 <- try(nls(Result$InclusionMean ~ (a*Result$SampleSize)/(1+b*Result$SampleSize), data = Result, start = list(a = 1,b = 0.1)), silent = TRUE)
  if (class(M1) != "try-error") {     ### run this only if nls was successful
    Result$pred <- predict(M1)
    P2 <- aggregate(pred ~ SampleSize, Result, FUN = mean)
    P2$sd <- aggregate(InclusionMean ~ SampleSize, Result, FUN = sd)[,2]
    plot(InclusionMean ~ SampleSize, data = Result, 
         pch = 16, cex = 0.2, col = "darkgray", ylim = c(0, 1),  
         ylab = "Inclusion", xlab = "Sample Size", main = paste(unique(TracksSpatial$Population), "udLev", udLev, sep = "_"))
    yTemp <- c((P2[,2] + P2[,3]), rev(P2[,2] - P2[,3]))
    xTemp <- c(P2[,1], rev(P2[,1]))
    polygon(x = xTemp, y = yTemp, col = "gray93", border = F)
    points(InclusionMean ~ SampleSize, data = Result, pch = 16, cex = 0.2, col = "darkgray")
    lines(P2, lty = 1,lwd = 2)
    Asymptote <- (summary(M1)$coefficients[1]/summary(M1)$coefficients[2])
    RepresentativeValue <- max(P2$pred)/Asymptote*100
    Result$RepresentativeValue <- RepresentativeValue
    print(RepresentativeValue)
    text(x = 2.5, y = 0.9,paste(round(RepresentativeValue,2), "%", sep = ""), cex = 2, col = "gray45", adj = 0)
  } else{RepresentativeValue <- mean(Result$InclusionMean[Result$SampleSize == max(Result$SampleSize)])   ### if nls is unsuccessful then use mean output for largest sample size
         Result$RepresentativeValue <- (RepresentativeValue/(udLev/100))*100
         Result$pred <- NA}    ## added by Jono Handley to convert to same scale as nls output
  
  Result$Asymptote <- Asymptote
  write.table(Result, paste(unique(TracksSpatial$Population), "UDLEv", udLev, "bootout_temp.csv", sep = "_"), row.names = F, sep = ",")
  return(Result)
}
