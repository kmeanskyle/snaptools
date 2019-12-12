# Goal of this script is to compare observations with data extracted
#   by snaptools functions

# prep daily observed data
prep_daily <- function(fn) {
  # only need these vars
  vars <- c("DATE","TMIN", "TMIN_ATTRIBUTES")
  # better names
  bnames <- c("date", "tmin", "tmin_attr")
  nome <- fread(fn, select = vars, col.names = bnames)

  # convert to correct type and units (m and C) and
  #   subset to matching time frame
  begin <- ymd("1979-01-01")
  end <- ymd("1979-12-31")
  nome[, ':=' (date = ymd(date), tmin = as.numeric(tmin)/10)]
  nome <- nome[date >= begin & date <= end, ]
}

# calculate euclidean distance between two vectors
eudist <- function(x, y) {
  sqrt(sum((x - y)^2))
}

library(data.table)
library(lubridate)
devtools::load_all()

fn1 <- "data-raw/Nome_daily.csv"
fn2 <- "data-raw/ _daily.csv"
nome <- prep_daily(fn1)

nc_fn <- "data-raw/t2min_daily_wrf_GFDL-CM3_historical_1979.nc"
nome_coords <- ak_coords[1, ]
# list of nome and surrounding cells
wrf_lst <- list(
  nome = wrf_get(nc_fn, nome_coords),
  nome_1N = wrf_get(nc_fn, nome_coords, c(0, 1)),
  nome_1E = wrf_get(nc_fn, nome_coords, c(1, 0)),
  nome_1S = wrf_get(nc_fn, nome_coords, c(0, -1)),
  nome_1W = wrf_get(nc_fn, nome_coords, c(-1, 0))
)

out <- lapply(wrf_lst, eudist, nome$tmin)

