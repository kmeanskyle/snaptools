#-- WRF Land Mask -------------------------------------------------------------
# make a land mask from an ancillary data file and a sample dataset from SNAP
# ancillary file: geo_em.d01.nc
# sample data file: t2min_daily_wrf_GFDL-CM3_historical_1979.nc
mk_wrf_mask <- function(anc_fn, sam_fn) {
  nc <- ncdf4::nc_open(anc_fn)
  mask <- ncvar_get(nc, "LANDMASK")
  ncdf4::nc_close(nc)
  # get spatial info not present in ancillary file
  nc <- ncdf4::nc_open(sam_fn)
  proj_str <- ncdf4::ncatt_get(nc, "xc", "proj_parameters")$value
  xc <- nc$dim$xc$vals
  yc <- nc$dim$yc$vals
  ncdf4::nc_close(nc)
  # make the mask
  r <- (xc[2] - xc[1])/2
  ext <- c(min(xc) - r, max(xc) + r, min(yc) - r, max(yc) + r)
  list(
    raster::raster(
      apply(t(mask), 2, rev),
      xmn = ext[1], xmx = ext[2], ymn = ext[3], ymx = ext[4]
    ),
    proj_str,
    list(xc = xc, yc = yc)
  )
}

library(ncdf4)

fn1 <- "data-raw/geo_em.d01.nc"
fn2 <- "data-raw/t2min_daily_wrf_GFDL-CM3_historical_1979.nc"

wrf_mask <- mk_wrf_mask(fn1, fn2)
wrf_xcyc <- wrf_mask[[3]]
wrf_proj_str <- wrf_mask[[2]]
wrf_mask <- wrf_mask[[1]]

#------------------------------------------------------------------------------



#-- WRF Truthing Data ---------------------------------------------------------
# make internal datasets from correct grid cells for Fairbanks and Nome
# sample data file: t2min_daily_wrf_GFDL-CM3_historical_1979.nc
mk_wrf_data <- function(sam_fn, coords) {
  # determine expected grid cell indices
  xy <- sp::SpatialPoints(coords, sp::CRS("+init=epsg:4326"))
  xy <- sp::spTransform(xy, sp::CRS(wrf_proj_str))
  xy <- unlist(as.data.frame(xy))

  # Set up grid_cps and grid_arr such that row number is the mapping between
  #   centerpoint coordinates and position in the output grid.
  grid_cps <- expand.grid(wrf_xcyc$xc, wrf_xcyc$yc)
  grid_arr <- matrix(1:68644, nrow = 262)
  # Then, find row number of coordinates nearest to "coords", get row/col
  eudist <- function(cp, xy) {
    sqrt((cp[1] - xy[1])^2 + (cp[2] - xy[2])^2)
  }
  dists <- apply(grid_cps, 1, eudist, xy)
  arr_ind <- which(
    grid_arr == which(dists == min(dists)),
    arr.ind = TRUE
  )

  nc <- ncdf4::nc_open(sam_fn)
  ncdat <- ncdf4::ncvar_get(
    nc, start = c(arr_ind[1], arr_ind[2], 1),
    count = c(1, 1, -1)
  )
  ncdf4::nc_close(nc)
  ncdat
}

library(ncdf4)

sam_fn <- "data-raw/t2min_daily_wrf_GFDL-CM3_historical_1979.nc"

truth <- list()
truth$nome <- mk_wrf_data(sam_fn, ak_coords[1, ])
truth$fairbanks <- mk_wrf_data(sam_fn, ak_coords[2, ])
truth$utqiagvik <- mk_wrf_data(sam_fn, ak_coords[3, ])
truth$juneau <- mk_wrf_data(sam_fn, ak_coords[4, ])

# save the following for internal use:
# land mask for grid
# proj4 string for grid
# truthing datasets for Nome, Fairbanks, Utqiagvik,
#   and Juneau (CM3 t2min, 1979)
usethis::use_data(
  wrf_mask,
  wrf_proj_str,
  truth,
  internal = TRUE, overwrite = TRUE
)

#------------------------------------------------------------------------------
