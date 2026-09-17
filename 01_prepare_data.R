# ─────────────────────────────────────────────────────────────────────────────
# 01 — Download and prepare PIFSC WinCruz cetacean sighting records
#
# Produces three objects consumed by 02_seasonal_map.R:
#   usable   sightings that can support density estimation
#   hi_eez   Hawaiian Islands EEZ polygon
#   lookup   species -> PhyloPic silhouette uuid
# ─────────────────────────────────────────────────────────────────────────────
library(dplyr); library(readr); library(sf); library(here)

dir.create(here("data", "raw"),     recursive = TRUE, showWarnings = FALSE)
dir.create(here("data", "derived"), recursive = TRUE, showWarnings = FALSE)

# ── download ─────────────────────────────────────────────────────────────────
TGZ <- here("data", "raw", "wincruz.tgz")
SRC <- "https://oceanwatch.pifsc.noaa.gov/xfer/PIFSC_PIRO_bulk_data_download_InPort_18141.tgz"

if (!file.exists(TGZ)) {
  download.file(SRC, TGZ, mode = "wb")
  untar(TGZ, exdir = here("data", "raw", "wincruz"))
}

csv <- list.files(here("data", "raw", "wincruz"), "\\.csv$",
                  recursive = TRUE, full.names = TRUE)[1]
stopifnot(length(csv) == 1, file.exists(csv))

# Everything as character, converted deliberately below. The file uses the
# string "n/a" for missing values in numeric-looking columns, and SpeciesCode
# mixes numeric cetacean codes with alphabetic ones for turtles and monk seal.
d <- read_csv(csv, col_types = cols(.default = col_character()))

# ── classify ─────────────────────────────────────────────────────────────────
YEAR_RANGE <- c(1810, 2020)   # unused here; kept for symmetry with other projects

# Categories that are records of "a cetacean" but not of a species. They cannot
# carry a density estimate, so they are excluded from the usable set.
UNIDENTIFIED <- c(
  "Unid. Dolphin", "Unid. small delphinid", "Unid. medium delphinid",
  "Unid. large delphinid", "Unid. Beaked whale", "Mesoplodon beaked whale",
  "Unid. Rorqual", "Unid. small whale", "Unid. large whale", "Unid. Whale",
  "Unid. Cetacean", "Sei/Bryde's whale", "Pygmy/dwarf sperm whale"
)

d3 <- d |>
  mutate(
    date   = lubridate::mdy(Date),
    year   = lubridate::year(date),
    month  = lubridate::month(date),
    best   = suppressWarnings(as.numeric(na_if(TotalBest, "n/a"))),
    lon    = suppressWarnings(as.numeric(na_if(ShipBeginLong, "n/a"))),
    lat    = suppressWarnings(as.numeric(na_if(ShipBeginLat,  "n/a"))),
    lon360 = ifelse(lon < 0, lon + 360, lon),
    code   = suppressWarnings(as.numeric(SpeciesCode)),
    # Turtles (CC, CM, DC, EI) and monk seal (MS) use alphabetic codes.
    is_cetacean   = !is.na(code),
    is_identified = !CommonName %in% UNIDENTIFIED
  )

# ── the funnel ───────────────────────────────────────────────────────────────
# EffortType, per the InPort entity metadata: systematic (S), non-systematic
# (N), fine-scale (F), off (O). Only S is line-transect effort. "O" looks like
# "on" and is not.
keep <- with(d3, is_cetacean & is_identified &
               EEZ %in% "Hawaii" & EffortType %in% "S" & !is.na(best))

funnel <- tibble::tibble(
  step = c("all records", "cetaceans", "identified to species",
           "Hawaii EEZ", "systematic effort", "with group size"),
  n = c(
    nrow(d3),
    sum(d3$is_cetacean),
    sum(d3$is_cetacean & d3$is_identified),
    sum(d3$is_cetacean & d3$is_identified & d3$EEZ %in% "Hawaii"),
    sum(d3$is_cetacean & d3$is_identified & d3$EEZ %in% "Hawaii" &
          d3$EffortType %in% "S"),
    sum(keep)
  )
)
print(funnel)

usable <- filter(d3, keep)
stopifnot(nrow(usable) > 0, sum(is.na(usable$lon)) == 0)

# ── what the map will show, and what it will not ─────────────────────────────
# Months with no systematic sightings are not months with no survey activity:
# March and April carry non-systematic records, December carries none at all.
effort_by_month <- d3 |>
  filter(EEZ %in% "Hawaii") |>
  count(month, EffortType) |>
  tidyr::pivot_wider(names_from = EffortType, values_from = n, values_fill = 0) |>
  arrange(month)
print(effort_by_month)

print(count(usable, CommonName, sort = TRUE), n = 40)

# ── EEZ boundary ─────────────────────────────────────────────────────────────
# Marine Regions EEZ v12, CC BY 4.0. Not redistributed here — see README.
eez <- st_read(here("World_EEZ_v12_20231025", "eez_v12.shp"), quiet = TRUE)
hi_eez <- filter(eez, GEONAME == "United States Exclusive Economic Zone (Hawaii)")
stopifnot(nrow(hi_eez) == 1)

# ── PhyloPic silhouettes ─────────────────────────────────────────────────────
spp <- count(usable, Species, CommonName, sort = TRUE)

lookup <- purrr::map_dfr(seq_len(nrow(spp)), function(i) {
  u <- tryCatch(rphylopic::get_uuid(name = spp$Species[i], n = 1),
                error = function(e) NA_character_)
  tibble::tibble(CommonName = spp$CommonName[i],
                 Species    = spp$Species[i],
                 n          = spp$n[i],
                 uuid       = if (length(u)) u[1] else NA_character_)
})
print(lookup, n = 30)

# ── save ─────────────────────────────────────────────────────────────────────
saveRDS(usable, here("data", "derived", "usable.rds"))
saveRDS(hi_eez, here("data", "derived", "hi_eez.rds"))
saveRDS(lookup, here("data", "derived", "lookup.rds"))
write_csv(funnel,          here("data", "derived", "funnel.csv"))
write_csv(effort_by_month, here("data", "derived", "effort_by_month.csv"))