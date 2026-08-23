# summarize_checkins.R
# Reduces loc-brightkite_totalCheckins.txt (~4.5M rows) to a small per-user
# homes file (~1 MB) that the CoLaS pipeline accepts directly.
#
# Run from the folder that contains the check-ins file, e.g.:
#   cd /Users/mariospapamichalis/Downloads
#   Rscript summarize_checkins.R
#
# Output: brightkite_homes.csv.gz  with columns
#   user, med_lat, med_lon, n_valid, n_checkins
# (median home over valid coordinates; n_valid = valid check-ins,
#  n_checkins = all check-in rows -- matching the pipeline's semantics).

infile  <- "loc-brightkite_totalCheckins.txt"
outfile <- "brightkite_homes.csv.gz"

if (requireNamespace("data.table", quietly = TRUE)) {
  library(data.table)
  ck <- fread(infile, header = FALSE, sep = "\t",
              col.names = c("user", "time", "lat", "lon", "loc_id"),
              select = c(1, 3, 4),
              colClasses = list(integer = 1, numeric = c(3, 4)))
  setnames(ck, c("user", "lat", "lon"))
  cnt_all <- ck[, .(n_checkins = .N), by = user]
  ok <- ck[is.finite(lat) & is.finite(lon) &
           !(lat == 0 & lon == 0) &
           lat >= -90 & lat <= 90 & lon >= -180 & lon <= 180]
  homes <- ok[, .(med_lat = median(lat), med_lon = median(lon),
                  n_valid = .N), by = user]
  homes <- merge(homes, cnt_all, by = "user")
  setorder(homes, user)
  fwrite(homes, outfile)   # .gz suffix -> gzip automatically
  cat(sprintf("wrote %s: %d users\n", outfile, nrow(homes)))
} else {
  # base-R fallback (slower, a few minutes; needs ~2 GB RAM)
  ck <- read.table(infile, sep = "\t", header = FALSE, quote = "",
                   comment.char = "",
                   colClasses = c("integer", "NULL", "numeric",
                                  "numeric", "NULL"))
  names(ck) <- c("user", "lat", "lon")
  n_checkins <- table(ck$user)
  ok <- ck[is.finite(ck$lat) & is.finite(ck$lon) &
           !(ck$lat == 0 & ck$lon == 0) &
           ck$lat >= -90 & ck$lat <= 90 &
           ck$lon >= -180 & ck$lon <= 180, ]
  agg <- aggregate(cbind(lat, lon) ~ user, data = ok, FUN = median)
  n_valid <- table(ok$user)
  agg$n_valid <- as.integer(n_valid[as.character(agg$user)])
  agg$n_checkins <- as.integer(n_checkins[as.character(agg$user)])
  names(agg) <- c("user", "med_lat", "med_lon", "n_valid", "n_checkins")
  agg <- agg[order(agg$user), ]
  con <- gzfile(outfile, "w")
  write.csv(agg, con, row.names = FALSE)
  close(con)
  cat(sprintf("wrote %s: %d users\n", outfile, nrow(agg)))
}
