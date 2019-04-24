### bootstrap colony effect ####

# Tracks: Tracking data. DF or SPDF, containing Latitude, Longitude, Species and Population
# SimulatedDistributions: DF or SPDF containing the simulated distribution of the studied 
#   species, with at least the columns Longitude, Latitude, Population and Species. Can be directly the output of the simulateDistribution() function
# tripID: unique identifier for each trip 
# UDLev: Numeric, Utilisation Distribution Level. Level of the KDE to be selected. For example, if 50 is used, the overlap is calculated between the 
#   smallest areas for which the probability to find the animal is equal to 0.50 (but see explanation for argument "conditional". From the adehabitatHR::kerneloverlap() function.
# Scale is the smoothing factor to be used in the Kernel Density Estimation (in Km)
# Grid is a number giving the size of the grid on which the UD should be estimated. 
# Iterations:  Numeric, number of iterations for the bootstrap

bootstrapColony <- function(Tracks, SimulatedDistributions, tripID, UDLev = 50, Scale, Grid = 500, Iterations = 50){
  require(sp)
  require(geosphere)
  require(rgdal)
  require(adehabitatHR)
  require(foreach)
  require(doParallel)
  require(parallel)
  
  Tracks$tripID <- Tracks[,tripID]
  if (!"Latitude" %in% names(Tracks)) stop("Latitude field does not exist")
  if (!"Longitude" %in% names(Tracks)) stop("Longitude field does not exist")
  if (!"tripID" %in% names(Tracks)) stop("tripID field does not exist")
  if (!"Population" %in% names(Tracks)) stop("Population field does not exist")
  if (!"Species" %in% names(Tracks)) stop("Sp field does not exist")
  
  # select species corresponding to the Tracking data
  SimulatedDistributions <- SimulatedDistributions[SimulatedDistributions$Species == unique(Tracks$Species),]
  
  # Converts SimulatedDistributions to spatial dataframe and projects it (and if already is SPDF check projection)
  if (class(SimulatedDistributions) != "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
  {
    ## filter DF to the minimum fields that are needed
    CleanSimulatedDistributions <- SimulatedDistributions %>%
      dplyr::select(Species, Population, Longitude, Latitude)
    mid_point <- data.frame(centroid(cbind(CleanSimulatedDistributions$Longitude, CleanSimulatedDistributions$Latitude)))
    
    ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
    if (min(CleanSimulatedDistributions$Longitude) < -170 &  max(CleanSimulatedDistributions$Longitude) > 170) {
      longs = ifelse(CleanSimulatedDistributions$Longitude < 0,CleanSimulatedDistributions$Longitude + 360,CleanSimulatedDistributions$Longitude)
      mid_point$lon <- ifelse(median(longs) > 180,median(longs) - 360,median(longs))}
    
    SimulatedDistributions.Wgs <- SpatialPoints(data.frame(CleanSimulatedDistributions$Longitude, CleanSimulatedDistributions$Latitude), proj4string = CRS("+proj=longlat + datum=wgs84"))
    proj.Sim <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep = ""))
    SimulatedDistributions.Projected <- spTransform(SimulatedDistributions.Wgs, CRS = proj.Sim)
    SimulatedDistributionsSpatial <- SpatialPointsDataFrame(SimulatedDistributions.Projected, data = CleanSimulatedDistributions)
    SimulatedDistributionsSpatial@data <- SimulatedDistributionsSpatial@data %>% dplyr::select(Species, Population)
    SimulatedDistributions.Wgs <- NULL
    SimulatedDistributions.Projected <- NULL
    
  } else {## if data are already in a SpatialPointsDataFrame then check for projection 
    if (is.na(SimulatedDistributions@proj4string)) stop("The proj4string slot of SimulatedDistributions can't be empty. Check and assign projection and re-run")
    if (!is.projected(SimulatedDistributions)) { # if it's not projected (lonlat) project to laea around midpoint
      mid_point <- data.frame(centroid(cbind(SimulatedDistributions@coords[,1], SimulatedDistributions@coords[,2])))
      proj.Sim <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep = ""))
      ### MB  This part prevents projection problems around the DATELINE 
      if (min(SimulatedDistributions@coords[,1]) < -170 &  max(SimulatedDistributions@coords[1,]) > 170) {
        longs = ifelse(SimulatedDistributions@coords[,1] < 0, SimulatedDistributions@coords[,1] + 360, SimulatedDistributions@coords[,1])
        mid_point$lon <- ifelse(median(longs) > 180, median(longs) - 360, median(longs))}
      SimulatedDistributionsSpatial <- spTransform(SimulatedDistributions, proj.Sim)
      SimulatedDistributionsSpatial@data <- SimulatedDistributionsSpatial@data %>% dplyr::select(Species, Population)
    } else {  ## if projected, "unproject" to WGS to project laea around midpoint
      SimulatedDistributions.WGS <- spTransform(SimulatedDistributions, CRS("+proj=longlat + datum=wgs84"))
      mid_point <- data.frame(centroid(cbind(SimulatedDistributions.WGS@coords[,1], SimulatedDistributions.WGS@coords[,2])))
      
      ### MB  This part prevents projection problems around the DATELINE 
      if (min(SimulatedDistributions.WGS@coords[,1]) < -170 &  max(SimulatedDistributions.WGS@coords[1,]) > 170) {
        longs = ifelse(SimulatedDistributions@coords[,1] < 0, SimulatedDistributions@coords[,1] + 360, SimulatedDistributions@coords[,1])
        mid_point$lon <- ifelse(median(longs) > 180, median(longs) - 360, median(longs))}
      
      proj.Sim <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep = ""))
      SimulatedDistributionsSpatial <- spTransform(SimulatedDistributions.WGS, CRS = proj.Sim)
      SimulatedDistributionsSpatial@data <- SimulatedDistributionsSpatial@data %>% dplyr::select(Species, Population)
    }
  }
  
  # Converts Tracks to spatial dataframe and projects it (and if already is SPDF check projection)
  if (class(Tracks) != "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
  {
    ## filter DF to the minimum fields that are needed
    CleanTracks <- Tracks %>%
      dplyr::select(Species, Population, tripID, Latitude, Longitude)
    
    Tracks.Wgs <- SpatialPoints(data.frame(CleanTracks$Longitude, CleanTracks$Latitude), proj4string = CRS("+proj=longlat + datum=wgs84"))
    Tracks.Projected <- spTransform(Tracks.Wgs, CRS = SimulatedDistributionsSpatial@proj4string)
    TracksSpatial <- SpatialPointsDataFrame(Tracks.Projected, data = CleanTracks)
    TracksSpatial@data <- TracksSpatial@data %>% dplyr::select(Species, Population, tripID)
    Tracks.Wgs <- NULL
    Tracks.Projected <- NULL
    
  } else {## if data are already in a SpatialPointsDataFrame then check for projection
    if (is.na(Tracks@proj4string)) stop("The proj4string slot of Tracks can't be empty. Check and assign projection and re-run")
    TracksSpatial <- spTransform(Tracks, CRS = SimulatedDistributionsSpatial@proj4string)
    TracksSpatial@data <- TracksSpatial@data %>% dplyr::select(Species, Population, tripID)
  }
    
  
  # remove tripID tracks with < 6 points as they can't be used to calculate kernel
  UIDs <- names(which(table(TracksSpatial@data[, tripID]) > 6))             
  TracksSpatial <- TracksSpatial[TracksSpatial@data[, tripID] %in% UIDs, ]
  TracksSpatial@data[ ,tripID] <- droplevels(as.factor(TracksSpatial@data[ ,tripID]))
  
  TracksSpatial$X <- TracksSpatial@coords[,1]
  TracksSpatial$Y <- TracksSpatial@coords[,2]
  
  UIDs <- as.character(unique(TracksSpatial$tripID))
  Ntrips <- length(UIDs)
  Nloop <- seq(1,(Ntrips - 1),ifelse(Ntrips > 100,10,1))
  DoubleLoop <- data.frame(SampleSize = rep(Nloop,each = Iterations), Iteration = rep(seq(1:Iterations),length(Nloop)))
  LoopNr <- seq(1:dim(DoubleLoop)[1])	
  UDLev <- UDLev
  
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
    SelectedCoords <- TracksSpatial[TracksSpatial@data[,tripID] %in% RanNum,]
    # Ext <- (min(SelectedCoords@coords[,1]) + 3 * diff(range(SelectedCoords@coords[,1])))
    # if(Ext < (Scale * 1000 * 2)) {
    # BExt <- ceiling((Scale * 1000 * 3)/(diff(range(SelectedCoords@coords[,1])))) #} else {BExt <- 3}
    KDE.Surface <- kernelUD(SelectedCoords, h = Scale*1000, grid = 500,  same4all = FALSE)
    try(KDE.UD <- getverticeshr(KDE.Surface, percent = UDLev))
    if (isTRUE(class(KDE.UD) == "try-error")) {
      sink("errors.txt", append = T)
      cat(paste("Failed in iteration", i, "with sample", RanNum, sep = " "))
      sink()} else {
        KDE.UD@proj4string <- SimulatedDistributionsSpatial@proj4string
        Overlain <- over(SimulatedDistributionsSpatial, KDE.UD)$area
        Output$InclusionMean <- length(Overlain[!is.na(Overlain)])/nrow(SimulatedDistributionsSpatial@data)
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
         ylab = "Inclusion", xlab = "Sample Size", main = paste(unique(TracksSpatial$Population), "UDLev", UDLev, sep = "_"))
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
         Result$RepresentativeValue <- (RepresentativeValue/(UDLev/100))*100
         Result$pred <- NA}    ## added by Jono Handley to convert to same scale as nls output
  
  Result$Asymptote <- Asymptote
  write.table(Result, paste(unique(TracksSpatial$Population), "UDLEv", UDLev, "bootout_temp.csv", sep = "_"), row.names = F, sep = ",")
  return(Result)
}
