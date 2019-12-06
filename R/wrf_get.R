# determine wrf cell indices for particular vector of WGS84 coordinates
#
# coords should be a data.frame with longitude, latitude as the two variables
wrf_cells <- function(coords) {
  # determine if matrix of wgs84 coords is valid and in the wrf grid
  # planned error messages:
  #   invalid coordinate(s): one or more longitude or latitudes is outside of allowable range
  #   not in WRF gird: one or more points falls outside of the downscaled grid (could just be a warning too that returns data for valid points)
  #

  is.data.frame(coords)
  abs(coords[, 1]) <= 180 & abs(coords[, 2]) <= 90

  coords_sp <- sp::SpatialPoints(coords, sp::CRS("+init=epsg:4326"))
  coords_sp <- sp::spTransform(coords_sp, sp::CRS(wrf_proj_str))
  # switch row/col because `wrf_mask` is transpose of grid in data files
  sw_cols <- function(rc) {
    cbind(rc[, 2], rc[, 1])
  }
  cell <- raster::cellFromXY(wrf_mask, coords_sp)
  row_col <- raster::rowColFromCell(cell, object = wrf_mask)
  sw_cols(row_col)
}

#' Access the "Historical and Projected Dynamically Downscaled Climate Data
#' for the State of Alaska and surrounding regions" datasets produced by SNAP
#'
#' @param nc_fns A netcdf file name (path) containing the downscaled data,
#' or vector of such filenames
#' @param coords A two_variable data.frame of WGS84 coordinates
#' @examples
#' #
#' # coordinates for Nome, AK
#' #
#' nome <- ak_coords[1, ]
#' fn <- "t2min_daily_wrf_GFDL-CM3_historical_1979.nc"
#' wrf_get(fn, "t2min", nome)
#' @export
wrf_get <- function(nc_fns, coords, shift = NULL) {
  var_eq <- function() {
    vars <- lapply(nc_fns, function(x) {
      nc <- ncdf4::nc_open(x)
      var <- nc$var$t2min$name
      ncdf4::nc_close(nc)
      var
    })
    vars <- unique(unlist(vars))
    if(length(vars) == 1) vars else FALSE
  }

  mk_tvec <- function(nc) {
    tstr <- strsplit(nc[[14]][[1]]$dim[[3]]$units, " ")
    tn <- nc[[14]][[1]]$varsize[3]
    tunit <- tstr[[1]][1]
    torig <- tstr[[1]][3]
    if(tunit == "days") {
      tvec <- as.Date(0:(tn - 1), origin = torig)
    } else {
      tseq <- seq(0, (tn - 1) * 3600, by = 3600)
      tvec <- as.POSIXlt(tseq, tz = "GMT", origin = torig)
    }
  }

  wrap_ncvar_get <- function(nc_fn, wrf_ijs) {
    ncvar_get_ij <- function(wrf_ij) {
      ncdf4::ncvar_get(nc, start = c(wrf_ij, 1), count = c(1, 1, -1))
    }
    n <- length(wrf_ijs)
    nc <- ncdf4::nc_open(nc_fn)
    nc_df <- lapply(wrf_ijs, ncvar_get_ij)
    nc_df <- as.data.frame(matrix(unlist(nc_df), ncol = n))
    rownames(nc_df) <- mk_tvec(nc)
    ncdf4::nc_close(nc)
    nc_df
  }

  # really need to split nc_fns on "/" and compare uniqueness of actual filenames not including paths
  if(length(nc_fns) > length(unique(nc_fns))) {
    warning("Removing duplicated file paths")
    nc_fns <- unique(nc_fns)
  }

  var <- var_eq()
  if(var == FALSE) stop()
  wrf_ijs <- wrf_cells(coords)
  # shift ijs
  if(!is.null(shift)){
    wrf_ijs[, 1] <- wrf_ijs[, 1] + shift[1]
    wrf_ijs[, 2] <- wrf_ijs[, 2] + shift[2]
  }
  wrf_ijs <- split(, row(coords))
  df <- do.call("rbind", lapply(nc_fns, wrap_ncvar_get, wrf_ijs))
  cols <- eval(var)
  n <- ncol(df)
  if(n > 1) {
    cols <- paste0(rep(eval(var), n), ".", 1:n)
  }
  names(df) <- cols
  df
}
