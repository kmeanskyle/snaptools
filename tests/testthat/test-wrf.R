context("wrf")
library(snaptools)

fn <- "../../data-raw/t2min_daily_wrf_GFDL-CM3_historical_1979.nc"

strp <- function(x) unname(as.array(unlist(x)))

test_that("get_wrf extracts correct data", {
  expect_equal(strp(wrf_get(fn, ak_coords[1, ])), truth$nome)
  expect_equal(strp(wrf_get(fn, ak_coords[2, ])), truth$fairbanks)
  expect_equal(strp(wrf_get(fn, ak_coords[3, ])), truth$utqiagvik)
  expect_equal(strp(wrf_get(fn, ak_coords[4, ])), truth$juneau)
})
