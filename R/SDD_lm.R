#' Spatial Double Differences (any direction)
#'
#' Spatial double difference (SDD) model for cross-sectional data estimated using the lm function on transformed data. Please see: Druckenmiller & Hsiang (2018).
#' @param spatial_df SpatialPolygonsDataFrame that includes your observational unit, dependent variable, and independent variable(s). This object should contain all observational units in your geography (e.g. counties in the US), not just those with complete cases.
#' @param n_channel numeric that specifies the number sampling channels you would like to use when calculating SFD. This number should be chosen such that the height of your sampling channels is approximately equal to the height of your observational units
#' @param obs_var obseravational unit 
#' @param dependent_var dependent variable 
#' @param independent_vars independent variable(s)
#' @param rotation numeric that specifies the angle at which SFD will be taken (-90 to 90 degrees), default=0. 
#' @param needs_balanced logical; if TRUE, the regression will be conducted on balanced data
#' @param balanced_df data.frame containing only those observations included in the balanced data. This data frame must include your observation variable
#' @param plot logical; if TRUE, it will plot the adjacent observations used to calculate SFD in your sample
#' @keywords SDD
#' @export
#' @examples
#' SDD_lm()

SDD_lm <- function(spatial_df, n_channel,
  obs_var, dependent_var, independent_vars,
  rotation=0, needs_balanced = FALSE, balanced_df = NA, plot = TRUE){
  print(noquote("Computing SDD..."))
  library(magrittr, quietly = "true")
  library(sp, quietly = "true")
  # Check that variables are of correct type
  if(!class(spatial_df)=="SpatialPolygonsDataFrame") stop('spatial_df must be a SpatialPolygonsDataFrame')
  if(!class(n_channel)=="numeric") stop('n_channel must be a numeric')
  if(!class(obs_var)=="character") stop('obs_var must be a character')
  if(!class(dependent_var)=="character") stop('dependent_var must be a character')
  if(!class(independent_vars)=="character") stop('independent_vars must be a character or a list of characters (e.g. independent_vars=c("age", "gender"))')
  if(!class(needs_balanced)=="logical") stop('needs_balanced must be a logical')
  if(needs_balanced==TRUE & !class(balanced_df)=="data.frame") stop('balanced_df must be a data.frame containing the obs_var')
  if(!class(plot)=="logical") stop('plot must be a logical')
  if(!class(rotation)=="numeric") stop('rotation must be a numeric')
  # Add centroids and define addtional variables
  cent <- as.data.frame(rgeos::gCentroid(spatial_df, byid=TRUE))
  spatial_df$cent.long <- cent[,1]
  spatial_df$cent.lat <- cent[,2]
  nrow <- n_channel
  ncol <- nrow*2
  values <- nrow*ncol 
  nvars <- ncol(spatial_df)
  pb1 <- txtProgressBar(min = 1, max = (nrow-1), initial = 1, style = 3)
  # Rotate map by specificed degrees
  rotateMap = function(object, angle) {
    boxpts = sp::SpatialPoints(t(sp::bbox(object)), proj4string = CRS(proj4string(object)))
    boxLL = sp::bbox(sp::spTransform(boxpts, CRS("+init=epsg:4326")))
    llc = apply(boxLL, 1, mean)
    prj = paste0("+proj=omerc +lat_0=", llc[2], " +lonc=", llc[1], " +alpha=",
      angle, " +gamma=0.0 +k=1.000000 +x_0=0.000 +y_0=0.000 +ellps=WGS84 +units=m ")
    # return as a CRS:
    CRS(prj)
  }
  if (rotation != 0){
    new_map <- rotateMap(spatial_df, rotation)
    new_map <- sp::spTransform(spatial_df, new_map)
  }
  if (rotation == 0){
    new_map <- spatial_df
  }
  # Overlay map with grid
  z <- raster::raster(new_map, nrow=nrow, ncol=ncol, vals=1:values)
  z <- as(z, 'SpatialPolygonsDataFrame')
  # Construct a row
  k <- 0
  zrow <- z[((k*ncol)+1):((k*ncol)+ncol),]
  zrow <- rgeos::gBuffer(zrow, byid=TRUE, width=0)
  # Take intersection of that row and observations
  i <- new_map*zrow
  # Arrange counties in order (west to east)
  i <- as.data.frame(i)
  i <- dplyr::arrange(i, layer)
  # Remove duplicates
  sample <- i[!duplicated(eval(parse(text=paste("i", "$", obs_var, sep = "")))),]
  # Add variables to keep track of order
  rownames(sample) <- NULL
  order <- rownames(sample)
  sample <- cbind(order=order, sample)
  sample$order <- as.numeric(as.character(sample[,1]))
  sample$row <- k
  # Now we do this for each row using a loop
  for (k in 1:(nrow-1)){
    setTxtProgressBar(pb1,k)
    # Create row
    zrow <- z[((k*ncol)+1):((k*ncol)+ncol),]
    zrow <- rgeos::gBuffer(zrow, byid=TRUE, width=0)
    # Take intersection of that row and observations
    i <- new_map*zrow
    i <- as.data.frame(i)
    # Drop any counties already sampled
    already.sampled <- unique(eval(parse(text=paste("sample", "$", obs_var, sep = ""))))
    i <- i[!(eval(parse(text=paste("i", "$", obs_var, sep = ""))) %in% already.sampled),]
    i <- i[!duplicated(eval(parse(text=paste("i", "$", obs_var, sep = "")))),]
    # Arrange counties in order (west to east)
    to_add <- dplyr::arrange(i, layer)
    # Add variables to keep track of order
    rownames(to_add) <- NULL
    order <- rownames(to_add)
    to_add <- cbind(order=order, to_add)
    to_add$order <- as.numeric(as.character(to_add[,1]))
    if (nrow(to_add) == 0) to_add[1,] <- matrix(NA,1,(nvars+2))
    to_add$row <- k
    # Add the sample from this row to the sample so far
    sample <- rbind(sample, to_add)
    not.missing <- complete.cases(sample[,eval(obs_var)])
    sample <- sample[not.missing,]
  }
  # Make sure the sample is in order
  sample <- dplyr::arrange(sample, row, order)
  # Isolate variables of interest
  vars <- sample %>% dplyr::select_(.dots = independent_vars) %>% as.matrix
  vars <- cbind(as.data.frame(eval(parse(text=paste("sample", "$", dependent_var, sep = "")))), vars)
  vars <- as.data.frame(vars)
  colnames(vars)[1] <- eval(dependent_var)
  kvars <- ncol(vars)
  ksample <- ncol(sample)
  rowvar <- ncol(sample)
  # Take spatial first differences
  to_subtract <- rbind(NA, as.data.frame(vars[1:nrow(vars)-1,1]))
  sample <- cbind(sample, (vars[,1]-to_subtract))
  colnames(sample)[(ksample + 1)] <- paste("d.", dependent_var, sep = "")
  for (r in 2:kvars){
    to_subtract <- rbind(NA, as.data.frame(vars[1:nrow(vars)-1,r]))
    sample <- cbind(sample, (vars[,r]-to_subtract))
    colnames(sample)[(ksample + r)] <- paste("d.", independent_vars[(r-1)], sep="")
  }
  to_subtract <- rbind(NA, as.data.frame(sample[1:nrow(sample)-1,ksample-1]))
  sample <- cbind(sample, (sample[,ksample-1]-to_subtract))
  colnames(sample)[(ksample + kvars + 1)] <- "d.layer"
  sample <- cbind(sample, rbind(NA, as.data.frame(eval(parse(text=paste("sample", "$", obs_var, sep = "")))[1:nrow(sample)-1])))
  colnames(sample)[ncol(sample)] <- "neighbor"
  # Remove first observation of each row (not a valid difference)
  sample <- subset(sample, order > 1)
  # Remove any observation that has big longitude jump (e.g. over a great lake)
  sample <- subset(sample, d.layer < 3)
  # If requested, include only observations in the balanced panel
  if (needs_balanced == TRUE){
    in_balanced <- unique(eval(parse(text=paste("balanced_df$", obs_var, sep = ""))))
    sample <- sample[(eval(parse(text=paste("sample$", obs_var, sep = ""))) %in% in_balanced),]
    sample <- sample[(sample$neighbor %in% in_balanced),]
  }
  # If requested, produce plot
  if (plot == TRUE){
    print(noquote("Plotting adjacent observations used for SFD..."))
    plot(spatial_df)
    title("Adjacent Observations used for SFD")
    bound <- nrow(sample)
    for (n in 1:bound){
      pair <- sample[n,]
      adj <- pair$neighbor
      match <- new_map[(eval(parse(text=paste("new_map$", obs_var, sep = ""))) %in% adj),]
      pair <- pair[,which(names(pair) %in% c("cent.long", "cent.lat"))]
      match <- match[,which(names(match) %in% c("cent.long", "cent.lat"))]
      match <- as.data.frame(match)
      pair <- rbind(pair, match)
      lines(pair$cent.long, pair$cent.lat, lwd = 1.5, col = "blue")
    }
  }
  # Take double differences
  # Isolate variables of interest
  explanatory_vars <- colnames(sample)[(ksample + 2)]
  if (kvars > 2){
    for (t in 3:kvars){
      add_on <- colnames(sample)[(ksample + t)]
      explanatory_vars <- rbind(explanatory_vars, add_on)
    }
  }
  explanatory_vars <- as.vector(explanatory_vars)
  vars <- sample %>% dplyr::select_(.dots = explanatory_vars) %>% as.matrix
  vars <- cbind(as.data.frame(eval(parse(text=paste("sample$d.", dependent_var, sep = "")))), vars)
  vars <- as.data.frame(vars)
  colnames(vars)[1] <- paste("d.", dependent_var, sep = "")
  kvars <- ncol(vars)
  ksample <- ncol(sample)
  # Take spatial double differences
  to_subtract <- rbind(NA, as.data.frame(vars[1:nrow(vars)-1,1]))
  sample <- cbind(sample, (vars[,1]-to_subtract))
  colnames(sample)[(ksample + 1)] <- paste("dd.", dependent_var, sep = "")
  for (r in 2:kvars){
    to_subtract <- rbind(NA, as.data.frame(vars[1:nrow(vars)-1,r]))
    sample <- cbind(sample, (vars[,r]-to_subtract))
    colnames(sample)[(ksample + r)] <- paste("dd.", independent_vars[(r-1)], sep="")
  }
  number <- nrow(sample)-2
  to_subtract <- rbind(NA, NA, as.data.frame(sample[1:number, rowvar]))
  sample <- cbind(sample, (sample[,rowvar]-to_subtract))
  colnames(sample)[ncol(sample)] <- "d.row"
  sample$d.row <- sample$row - lag(sample$row)
  sample$d.order <- sample$order - lag(sample$order)
  # Remove first observation of each row (not a valid difference)
  sample <- subset(sample, d.row < 1)
  sample <- subset(sample, d.order < 3)
  # Run regression
  explanatory_vars <- colnames(sample)[(ksample + 2)]
  if (kvars > 2){
    for (t in 3:kvars){
      add_on <- colnames(sample)[(ksample + t)]
      explanatory_vars <- rbind(explanatory_vars, add_on)
    }
  }
  explanatory_vars <- as.vector(explanatory_vars)
  explanatory_list <- paste(explanatory_vars, collapse = " + ")
  regression_results <- lm(
    as.formula(paste(paste(paste("dd.", dependent_var, sep = ""), "~", explanatory_list))),
    data = sample)
  return(regression_results)
}