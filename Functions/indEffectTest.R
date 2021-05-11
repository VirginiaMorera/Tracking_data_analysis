indEffectTest <- function(
  tracks, tripID, groupVar, plotIT = TRUE, 
  debias = TRUE) {
  
  #tracks: dataframe containing, at least,  latitude, longitude, trip ID, group ID and Date Time fields
  #      : datetime field must be called Date_Time and be in ymd_HMS format
  #      : latitude field must be called Latitude
  #      : longitude field must be called Longitude
  #      : the names of the trip and group identifier variables are specified in the arguments "tripID" and "groupID"
  
  # plotIT: logical, if TRUE, a density plot of the within- and between-group median overlaps will be produced. 
  # debias: logical, if TRUE debiased estimators of the BA will be produced. 
   
  
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("Package \"purrr\" needed for  function to work. Please install.",
         call. = FALSE)  }
  
  if (!requireNamespace("ctmm", quietly = TRUE)) {
    stop("Package \"ctmm\" needed for  function to work. Please install.",
         call. = FALSE)  }
  
  if (!(tripID) %in% names(tracks)) stop("Within-group field does not exist")
  if (!groupVar %in% names(tracks)) stop("Group field does not exist")
  
  tracks %>% 
    dplyr::select({{groupVar}}, {{tripID}}, .data$Date_Time, .data$Latitude, .data$Longitude) -> tracks
  
  # remove tripID tracks with < 6 points as can't be used in KDE ------------
  UIDs <- names(which(table(tracks[, tripID]) > 5))
  tracks <- tracks[tracks[, tripID] %in% UIDs, ]
  tracks[ ,tripID] <- droplevels(as.factor(tracks[ ,tripID]))
  
  # create vector with value of groupVar for each trip ------------------------
  tracks %>% 
    select({{tripID}}, {{groupVar}}) %>% 
    distinct() %>% 
    select({{groupVar}}) %>% 
    deframe() -> gid
  
  # turn tracks into telemetry object  ----------------------------------------
  tracks %>% 
    dplyr::group_by_(tripID) %>% 
    distinct(Date_Time, .keep_all = T) %>% 
    arrange(Date_Time) %>% 
    ungroup() %>% 
    as.data.frame() %>% 
    select(individual.local.identifier = {{tripID}}, datetime = .data$Date_Time, 
           .data$Longitude, .data$Latitude) -> tracks_tele
  
  
  tracks_tele <- as.telemetry(tracks_tele)
  
  
  FIT.list <- list()
  
  before <- Sys.time()
  for(i in 1:length(tracks_tele)){
    one <- tracks_tele[[i]]
    # guesstimate parameters  ----------------------------------------
    GUESS <- ctmm.guess(one, interactive=F)
    # fit competing models  ----------------------------------------
    FITS <- ctmm.select(one, GUESS,  trace=TRUE, verbose=TRUE, cores=-1)
    print(paste("The selected model for", names(tracks_tele)[i], "is", names(FITS)[1], sep = " "))
    # select best model  ----------------------------------------
    FIT.list[[i]] <- FITS[[1]]
    names(FIT.list)[i] <- names(tracks_tele)[i]
    }
  
  FIT.list <- purrr::compact(FIT.list)
  # calculate overlap between tracks ------------------------------------------
  ov <- ctmm::overlap(FIT.list, level = 0.95, debias = debias)
  
  # separate the layers of the median and CI -------------------------
  X05 <- ov[,,2]
  Xlo <- ov[,,1]
  Xhi <- ov[,,3]
  
  # remove lower triangle of each (repeated values since it's a symetric matrix) ----------------
  X05[lower.tri(X05, diag = T)] <- NA
  Xlo[lower.tri(Xlo, diag = T)] <- NA
  Xhi[lower.tri(Xhi, diag = T)] <- NA
  
  
  # assign value of BirdID that we rescue from the simulated_dataset df to rows and columns ------
  rownames(X05) <- colnames(X05) <- gid
  rownames(Xlo) <- colnames(Xlo) <- gid
  rownames(Xhi) <- colnames(Xhi) <- gid
  
  # median
  WI05 <- NULL
  BW05 <- NULL
  
  # BW lower and higher CI bounds
  BWlo <- NULL
  BWhi <- NULL
  
  for (i in seq_along(rownames(X05))) {
    # i = 2
    
    # medians
    x1 <- X05[i,] 
    x2 <- x1[which(names(x1) == rownames(X05)[i])]
    x3 <- x1[which(names(x1) != rownames(X05)[i])]
    
    WI05 <- c(WI05, x2)
    BW05 <- c(BW05, x3)
    
    # BW intervals
    y1 <- Xlo[i,] 
    y2 <- y1[which(names(y1) != rownames(Xlo)[i])]
    BWlo <- c(BWlo, y2)
    
    z1 <- Xhi[i,] 
    z2 <- z1[which(names(z1) != rownames(Xhi)[i])]
    BWhi <- c(BWhi, z2)
    
  }
  
  # remove NAs
  WI05 <- WI05[!is.na(WI05)]
  BW05 <- BW05[!is.na(BW05)]
  BWlo <- BWlo[!is.na(BWlo)]
  BWhi <- BWhi[!is.na(BWhi)]
  
  # organize values in a dataframe for plotting
  Overlaps <- data.frame(
    Overlap = c(WI05, BW05),
    Type = c(rep("Within", length(WI05)), rep("Between", length(BW05)))
  )
  
  
 
  # visualise BW and WI medians
  plot <- ggplot(data = Overlaps) + 
      geom_density(aes(x = Overlap, col = Type, fill = Type), alpha = 0.1) + 
      ggtitle("Density of within individual vs. between individual median overlaps") +
      theme_bw()
  
  if(plotIT == TRUE) {print(plot)}
  
  WI_median <- quantile(WI05, 0.5)
  
  BWint <- data.frame(Low = BWlo,
                      High = BWhi)
  
  BWint %>% 
    mutate(Contains = "NA") %>% 
    mutate(Contains = dplyr::if_else((Low > WI_median), "Below", 
                                     dplyr::if_else(Low < WI_median & High > WI_median, "Within", 
                                                    dplyr::if_else(High < WI_median, "Above", Contains)))) -> BWint
  
  proportions <- as.data.frame(prop.table(table(BWint$Contains)))
  
  names(proportions) <- c("WI_Median_Position", "Percentage")
  proportions$Percentage <- round(proportions$Percentage*100, 1)
  
  (proportions_wide <- spread(proportions, WI_Median_Position, Percentage))
  
  # organise output
  Result <- list()
  Result[1] <- list(ov)
  Result[2] <- list(plot)
  Result[3] <- list(proportions_wide)
  names(Result) <- c("OverlapArray", "Plot", "Proportions")

  cat("In ", proportions_wide[1, "Above"], "% of 'between group' overlaps the 'within group' median falls above the 95% CI", "\n")
  cat("In ", proportions_wide[1, "Within"], "% of 'between group' overlaps the 'within group' median falls inside the 95% CI", "\n")
  cat("In ", proportions_wide[1, "Below"], "% of 'between group' overlaps the within group median falls below the 95% CI", "\n")

  return(Result)
}
