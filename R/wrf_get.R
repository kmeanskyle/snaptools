# determine wrf cell indices for particular vector of WGS84 coordinates
#
# coords should be a data.frame with longitude, latitude as the two variables
wrf_cells <- function(coords) {
  # planned error messages:
  #   not in WRF gird: one or more points falls outside of the downscaled grid (could just be a warning too that returns data for valid points)
  #

  if(!is.data.frame(coords))
    stop("Please supply coordinates as a data.frame")

  if(!all(abs(coords[, 1]) <= 180 & abs(coords[, 2]) <= 90))
    stop("Invalid coordinate values, latitude should be between -90 and 90, longitude between -180 and 180")

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
#' @param coords A two-variable data.frame of WGS84 coordinates
#' @param shift A vector of horizontal (1st element) and vertical (2nd element) grid cell positions to shift before querying WRF output
#' @param rc A logical arg indicating whether supplied coords are row/column values
#' @examples
#' #
#' # coordinates for Nome, AK
#' #
#' nome <- ak_coords[1, ]
#' fn <- "t2min_daily_wrf_GFDL-CM3_historical_1979.nc"
#' wrf_get(fn, "t2min", nome)
#' @export
wrf_get <- function(nc_fns,
                    coords,
                    shift = NULL,
                    rc = FALSE) {

  verify_coords <- function() {
    if(!is.data.frame(rc_df))
      stop("Please supply coordinates or row/column values as a data.frame")
  }

  var_eq <- function() {
    vars <- lapply(nc_fns, function(x) {
      nc <- ncdf4::nc_open(x)
      var <- nc$var[[1]]$name
      ncdf4::nc_close(nc)
      var
    })
    vars <- unique(unlist(vars))
    if(length(vars) == 1) vars else
      stop("More than one unique WRF variable detected")
  }

  verify_rc <- function(coords) {
    valid <- 1:262
    if(!all(unlist(rc_df) %in% valid))
      stop("Invalid row/columns supplied. Valid values are integers in 1:262")
    coords
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
  verify_coords()

  wrf_ijs <- if(rc == FALSE) wrf_cells(coords) else verify_rc(coords)

  # shift ijs
  if(!is.null(shift)){
    wrf_ijs[, 1] <- wrf_ijs[, 1] + shift[1]
    wrf_ijs[, 2] <- wrf_ijs[, 2] - shift[2]
  }
  wrf_ijs <- split(wrf_ijs, row(coords))
  df <- do.call("rbind", lapply(nc_fns, wrap_ncvar_get, wrf_ijs))
  cols <- eval(var)
  n <- ncol(df)
  if(n > 1) {
    cols <- paste0(rep(eval(var), n), ".", 1:n)
  }
  names(df) <- cols
  df
}
