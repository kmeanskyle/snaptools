# save some sample coordinates for external use

# coordinates for select Alaskan communities (Google)
ak_coords <- data.frame(
  lon = c(-165.4046, -147.7164, -156.7886, -134.4197),
  lat = c(64.5011, 64.8378, 71.2906, 58.3019)
)
rownames(ak_coords) <- c("Nome", "Fairbanks", "Utqiagvik", "Juneau")

usethis::use_data(ak_coords, overwrite = TRUE)
