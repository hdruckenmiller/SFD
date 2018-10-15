#' Recommends number of sampling channels
#'
#' Recommends number of sampling channels to use when taking spatial first differences. Please see: Druckenmiller & Hsiang (2018).
#' @param spatial_df SpatialPolygonsDataFrame that includes your observational unit. This object should contain all observational units in your geography (e.g. counties in the US), not just those with complete cases.
#' @keywords SFD
#' @export
#' @examples
#' n_channels()

n_channels <- function(spatial_df){
  library(magrittr, quietly = "true")
  library(sp, quietly = "true")
  avg_dist <- mean(sqrt(raster::area(spatial_df))/1000)
  total_length <- (max(spatial_df$latitude) - min(spatial_df$latitude))*(111)
  n <- total_length/avg_dist
  return(n)}
