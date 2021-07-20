# Developed for the manuscript *Detecting recurrent sources of variability in animal tracking studies* by Virginia Morera-Pujol *et al.* submitted to Ecological Applications in 2021.

#### simulate distribution ####

# tracks must be a data frame or SPDF with at least Longitude and Latitude (not projected), Sp and Colony fields
# populationInfo must be a data frame with at least Colony, Sp and Pairs
# scale is the smoothing factor to be used in the Kernel Density Estimation (in Km)
# grid is a number giving the size of the grid on which the UD should be estimated. 
# multiFactor is the number by which we want to multiply each population (i.e. simulate that number of positions from each animal). The higher the number, the larger the resulting point pattern


simulateDistribution <- function(tracks, populationInfo, scale, grid = 500, multiFactor = 5){
  
  # tracks <- colony_boot_list[[2]] #this can be used for testing if something goes wrong
  # package dependencies
  if (!requireNamespace("adehabitatHR", quietly = TRUE)) {
    stop("Package \"adehabitatHR\" needed for  function to work. Please install.",
         call. = FALSE)  }
  if (!requireNamespace("sp", quietly = TRUE)) {
    stop("Package \"sp\" needed for  function to work. Please install.",
         call. = FALSE)  }
  if (!requireNamespace("spatstat", quietly = TRUE)) {
    stop("Package \"spatstat\" needed for  function to work. Please install.",
         call. = FALSE)  }
  if (!requireNamespace("geosphere", quietly = TRUE)) {
    stop("Package \"geosphere\" needed for  function to work. Please install.",
         call. = FALSE)  }
  
  # initial checks
  if (!"Latitude" %in% names(tracks)) stop("Latitude field does not exist")
  if (!"Longitude" %in% names(tracks)) stop("Longitude field does not exist")
  if (!"Species" %in% names(tracks)) stop("Species field does not exist")
  if (!"Population" %in% names(tracks)) stop("Population field does not exist")
  
  # Convert tracks to spatial dataframe and project them, or check projection if already SPDF (and if already is SPDF it accepts this)
  if (class(tracks) != "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
  {
    ## filter DF to the minimum fields that are needed
    CleanTracks <- tracks %>%
      dplyr::select(Species, Population, Latitude, Longitude)
    mid_point <- data.frame(centroid(cbind(CleanTracks$Longitude, CleanTracks$Latitude)))
    
    ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
    if (min(CleanTracks$Longitude) < -170 &  max(CleanTracks$Longitude) > 170) {
      longs = ifelse(CleanTracks$Longitude < 0,CleanTracks$Longitude + 360,CleanTracks$Longitude)
      mid_point$lon <- ifelse(median(longs) > 180,median(longs) - 360,median(longs))}
    
    Tracks.Wgs <- SpatialPoints(data.frame(CleanTracks$Longitude, CleanTracks$Latitude), proj4string = CRS("+proj=longlat + datum=wgs84"))
    proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep = ""))
    Tracks.Projected <- spTransform(Tracks.Wgs, CRS = proj.UTM )
    TracksSpatial <- SpatialPointsDataFrame(Tracks.Projected, data = CleanTracks)
    TracksSpatial@data <- TracksSpatial@data %>% dplyr::select(Species, Population, Latitude, Longitude)
    Tracks.Wgs <- NULL
    Tracks.Projected <- NULL
    
  }else {## if data are already in a SpatialPointsDataFrame then check for projection
    
    if (is.na(Tracks@proj4string)) stop("proj4string slot can't be NA. Assign the correct CRS object")
    
    if (is.projected(Tracks)) {
      TracksSpatial@data <- TracksSpatial@data %>% dplyr::select(Species, Population, Latitude, Longitude)
    }else {## project data to UTM if not projected
      mid_point <- data.frame(centroid(cbind(Tracks@data$Longitude, Tracks@data$Latitude)))
      
      ### MB  This part prevents projection problems around the DATELINE 
      if (min(Tracks@data$Longitude) < -170 &  max(Tracks@data$Longitude) > 170) {
        longs = ifelse(Tracks@data$Longitude < 0, Tracks@data$Longitude + 360, Tracks@data$Longitude)
        mid_point$lon <- ifelse(median(longs) > 180, median(longs) - 360, median(longs))}
      
      proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep = ""))
      TracksSpatial <- spTransform(Tracks, CRS = proj.UTM)
      TracksSpatial@data <- TracksSpatial@data %>% dplyr::select(GroupVar, tripID, Latitude, Longitude)
    }
  }
  
  map <- rworldmap::getMap(resolution = "coarse")
  map <- spTransform(map, TracksSpatial@proj4string)
  # generate kernel
  Kernel.est <- kernelUD(TracksSpatial,  h = scale*1000, grid = grid)
  
  # convert to pixel image
  r <- raster(as(Kernel.est, "SpatialPixelsDataFrame"))
  # raster.as.im function from Jeffrey Evans answer here: https://bit.ly/2TI0FXB
  kernel.im <- as.im(r)
  
  # select colony size info
  SPopulation <- as.character(unique(tracks$Population))
  SSpecies <- as.character(unique(tracks$Species))
  Pop.size <- populationInfo[populationInfo$Species == SSpecies & populationInfo$Population == SPopulation,]$Pairs*2
  
  # we're going to simulate a nº of points equal to the pop size * multiFactor
  SimulateN <- Pop.size*multiFactor
  
  # this simulates the points as ppp
  SimPoints <- rpoint(SimulateN, kernel.im)
  
  # convert to dataframe, and from there to Spatial points
  SimPoints.df <- as.data.frame(SimPoints)
  SimPoints.sp <- SimPoints.df
  coordinates(SimPoints.sp) <- ~ x+y
  
  # plot to see everything has worked
  par(mfrow = c(1,2))
  plot(kernel.im, main = SPopulation, xlim = SimPoints.sp@bbox[1,], ylim = SimPoints.sp@bbox[2,])
  plot(map, add = T, border = "white")
  # plot(kernel.im, main = SPopulation)
  plot(SimPoints.sp, pch = 20, col = "#ff000030", cex = 0.3)
  plot(map, add = T, border = "black")
  
  par(mfrow = c(1,1))
  
  # prepare output
  SimPoints.sp <- SpatialPointsDataFrame(coords = SimPoints.sp@coords, 
                                         data = data.frame(Population = rep(SPopulation, length(SimPoints.sp)), 
                                                           Species = rep(SSpecies, length(SimPoints.sp))),
                                         proj4string = proj.UTM)
  SimPoints.sp <- spTransform(SimPoints.sp, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  return(SimPoints.sp)
}
