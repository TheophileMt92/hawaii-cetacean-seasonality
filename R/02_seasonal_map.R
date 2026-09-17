# ─────────────────────────────────────────────────────────────────────────────
# Seasonal distribution of cetacean sightings, Hawaiian Islands EEZ
# Monthly map of systematic shipboard survey records, 2006-2017
# Expects in the environment: usable, hi_eez, lookup
# ─────────────────────────────────────────────────────────────────────────────
library(sf); library(dplyr); library(ggplot2)
library(marmap); library(rphylopic); library(terra)

# ── palette ──────────────────────────────────────────────────────────────────
# Sequential single hue, deep -> shallow, luminance monotonic (checked).
RAMP   <- c("#08141F", "#0E2739", "#163E55", "#1F5875", "#2E7796", "#4E9BB5", "#7CBFD4")
BG     <- "#06101A"
PT     <- "#FFC24D"   # 11.6:1 on deep water; the dark stroke carries it on the shelf
STROKE <- "#08141F"
INK    <- "#E9F1F5"   # 16.7:1
LABEL  <- "#D2DEE6"   # 14.0:1  species names
SOFT   <- "#B6C6D2"   # 10.9:1  captions, key, counts
RULE   <- "#9FD9E8"

# ── extents ──────────────────────────────────────────────────────────────────
# Map occupies 176-210; 210-215 is a clean band for the species column, since
# the bathymetry simply stops at 210 and the rest renders as background.
MAP_L <- 176; MAP_R <- 210
XLIM  <- c(MAP_L, 215)
YLIM  <- c(15, 32.5)

# ── geometry ─────────────────────────────────────────────────────────────────
# The Hawaii EEZ crosses the antimeridian (8 sightings sit east of it), so
# everything works in 0-360. In plain -180/180 the bbox wraps the whole globe.
hi <- st_shift_longitude(hi_eez)

hi_df <- st_coordinates(hi) |>
  as.data.frame() |>
  transmute(x = X, y = Y, ring = paste(L1, L2, L3)) |>
  group_by(ring) |>
  # Vertices sitting exactly on 180 are the antimeridian split in the source
  # polygon, not a maritime boundary. Drop them and break the path there.
  mutate(seam = abs(x - 180) < 1e-4,
         brk  = cumsum(seam)) |>
  ungroup() |>
  filter(!seam) |>
  mutate(grp = paste(ring, brk)) |>
  select(x, y, grp)

# Dim outside the EEZ rather than clipping: the Hawaiian Ridge does not stop at
# the boundary, and a hard-edged cut would read as a mistake.
sf_use_s2(FALSE)
mask_df <- tryCatch({
  panel <- st_as_sfc(st_bbox(c(xmin = MAP_L, xmax = MAP_R,
                               ymin = YLIM[1], ymax = YLIM[2]), crs = 4326))
  st_difference(panel, st_union(st_geometry(hi))) |>
    st_coordinates() |> as.data.frame() |>
    transmute(x = X, y = Y, grp = paste(L1, L2, L3))
}, error = function(e) { message("EEZ mask skipped: ", conditionMessage(e)); NULL })
sf_use_s2(TRUE)

pos <- usable |>
  mutate(
    lon    = suppressWarnings(as.numeric(na_if(ShipBeginLong, "n/a"))),
    lat    = suppressWarnings(as.numeric(na_if(ShipBeginLat,  "n/a"))),
    lon360 = ifelse(lon < 0, lon + 360, lon),
    month  = lubridate::month(lubridate::mdy(Date))
  )

# ── bathymetry ───────────────────────────────────────────────────────────────
bathy <- getNOAA.bathy(lon1 = MAP_L, lon2 = -150, lat1 = YLIM[1], lat2 = YLIM[2],
                       resolution = 4, antimeridian = TRUE, keep = TRUE)

raw <- fortify.bathy(bathy) |>
  mutate(x = ifelse(x < 0, x + 360, x)) |>
  filter(between(x, MAP_L, MAP_R), between(y, YLIM[1], YLIM[2]))

# marmap stitches two halves either side of 180 with DIFFERENT cell sizes
# (0.06694 vs 0.06793, plus a 0.035 sliver at the seam). geom_raster assumes a
# uniform grid, so the mismatch renders as vertical striping. Resample onto one
# regular grid at 0.07 deg — coarser than both source spacings, so no holes.
# No crs: we never reproject, and asking for one trips over Anaconda's proj.db.
target <- rast(xmin = MAP_L, xmax = MAP_R, ymin = YLIM[1], ymax = YLIM[2],
               resolution = 0.07)
bath_r <- rasterize(as.matrix(raw[, c("x", "y")]), target,
                    values = raw$z, fun = mean)

bdf <- as.data.frame(bath_r, xy = TRUE)
names(bdf)[3] <- "z"          # layer name varies by terra version
sea  <- filter(bdf, z <= 0)
land <- filter(bdf, z >  0)

# If this ever fails, the grid is irregular again and the frames will stripe.
stopifnot(length(unique(round(diff(sort(unique(bdf$x))), 6))) == 1)

# ── silhouettes ──────────────────────────────────────────────────────────────
FALLBACK <- c("Stenella attenuata"    = "Stenella",
              "Balaenoptera edeni"    = "Balaenoptera",
              "Indopacetus pacificus" = "Ziphiidae")

sil <- lapply(seq_len(nrow(lookup)), function(i) {
  id <- lookup$uuid[i]
  if (is.na(id)) {
    nm <- FALLBACK[[lookup$Species[i]]]
    id <- tryCatch(get_uuid(name = nm, n = 1)[1], error = function(e) NA_character_)
  }
  if (is.na(id)) return(NULL)
  tryCatch(get_phylopic(id), error = function(e) NULL)
})
names(sil) <- lookup$CommonName

# Licence condition, not a courtesy. Put these in the repo README.
attrib_tbl <- purrr::map_dfr(na.omit(lookup$uuid), function(u) {
  a <- tryCatch(get_attribution(uuid = u), error = function(e) NULL)
  if (is.null(a)) return(NULL)
  tibble::tibble(uuid = u,
                 contributor = a$contributor %||% NA_character_,
                 license     = a$license     %||% NA_character_)
})
readr::write_csv(attrib_tbl, "phylopic_credits.csv")

# ── credits ──────────────────────────────────────────────────────────────────
# NMFS permit numbers are a binding condition of the data licence.
PERMITS <- paste(sort(unique(usable$Permit)), collapse = ", ")

CAPTION <- paste0(
  "Sightings: NOAA PIFSC Cetacean Research Program, WinCruz 2006\u20132017 (InPort 18141). NMFS permits: ", PERMITS, "\n",
  "Bathymetry: NOAA ETOPO via marmap \u00b7 EEZ: Marine Regions v12 \u00b7 ",
  "Silhouettes: PhyloPic (per-image credits in repo); some are genus-level stand-ins"
)

# ── layout ───────────────────────────────────────────────────────────────────
SIL_X   <- 212.4      # species column, centred in the right-hand band
SIL_TOP <- 31.3
SIL_GAP <- 2.60
SIL_H   <- 0.70
N_SIL   <- 5

key <- data.frame(x = seq(177, 181.5, length.out = 60), y = 16.4) |>
  mutate(z = seq(-6000, 0, length.out = 60))

MONTHS <- month.name

# ── one frame ────────────────────────────────────────────────────────────────
frame_plot <- function(m) {
  pts <- filter(pos, month == m)
  cnt <- count(pts, CommonName, sort = TRUE)
  top <- head(cnt, N_SIL)

  p <- ggplot() +
    geom_raster(data = sea,  aes(x, y, fill = z)) +
    geom_raster(data = land, aes(x, y), fill = "#20242A") +
    scale_fill_gradientn(colours = RAMP, limits = c(-6000, 0),
                         oob = scales::squish, guide = "none")

  if (!is.null(mask_df))
    p <- p + geom_polygon(data = mask_df, aes(x, y, group = grp),
                          fill = BG, alpha = .5)

  p <- p + geom_path(data = hi_df, aes(x, y, group = grp),
                     colour = RULE, linewidth = .3, alpha = .85)

  if (nrow(pts))
    p <- p + geom_point(data = pts, aes(lon360, lat), shape = 21,
                        fill = PT, colour = STROKE, stroke = .35,
                        size = 2.6, alpha = .95)

  if (nrow(top)) {
    for (i in seq_len(nrow(top))) {
      yy  <- SIL_TOP - (i - 1) * SIL_GAP
      img <- sil[[top$CommonName[i]]]
      if (!is.null(img))
        p <- p + add_phylopic(img = img, x = SIL_X, y = yy,
                              height = SIL_H, fill = INK, alpha = .92)
      p <- p + annotate("text", x = SIL_X, y = yy - 0.80,
                        label = stringr::str_wrap(top$CommonName[i], 14),
                        hjust = .5, vjust = 1, colour = LABEL,
                        size = 3.4, lineheight = .92)
    }
    if (nrow(cnt) > N_SIL)
      p <- p + annotate("text", x = SIL_X, y = SIL_TOP - N_SIL * SIL_GAP,
                        label = paste0("+ ", nrow(cnt) - N_SIL, " more species"),
                        hjust = .5, colour = SOFT, size = 3.2, fontface = "italic")
  }

  p +
    geom_tile(data = key, aes(x, y, fill = z), height = .32) +
    annotate("text", x = 177, y = 17.15, label = "Depth",
             hjust = 0, colour = SOFT, size = 3.6) +
    annotate("text", x = 177,   y = 15.8, label = "6,000 m", hjust = 0,
             colour = SOFT, size = 3.4) +
    annotate("text", x = 181.5, y = 15.8, label = "0", hjust = 1,
             colour = SOFT, size = 3.4) +
    annotate("text", x = MAP_R - .5, y = 16.8, label = toupper(MONTHS[m]),
             hjust = 1, colour = INK, size = 12, fontface = "bold") +
    annotate("text", x = MAP_R - .5, y = 15.9,
              label = if (nrow(pts)) paste(nrow(pts), "sightings")
                     else "no systematic sightings on record",
             hjust = 1, colour = if (nrow(pts)) SOFT else PT, size = 5) +
    coord_fixed(ratio = 1 / cos(mean(YLIM) * pi / 180),
                xlim = XLIM, ylim = YLIM, expand = FALSE) +
    labs(
      title    = "Seasonal distribution of cetacean sightings, Hawaiian Islands EEZ",
      subtitle = "Systematic shipboard survey records, 2006–2017, pooled by month.\nBlank months had no systematic line-transect effort, not necessarily no animals.",
      caption  = CAPTION, x = NULL, y = NULL
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.background = element_rect(fill = BG, colour = NA),
      plot.title      = element_text(colour = INK, face = "bold", size = 19,
                                     margin = margin(b = 4)),
      plot.subtitle   = element_text(colour = SOFT, size = 13, margin = margin(b = 10)),
      plot.caption    = element_text(colour = SOFT, size = 9, hjust = 0,
                                     lineheight = 1.35, margin = margin(t = 10)),
      plot.margin     = margin(18, 20, 14, 20)
    )
}

# ── check before rendering 240 frames ────────────────────────────────────────
ggsave("test_sep.png", frame_plot(9), width = 14.0, height = 8.2, dpi = 100)
ggsave("test_mar.png", frame_plot(3), width = 14.0, height = 8.2, dpi = 100)

# ── video ────────────────────────────────────────────────────────────────────
frames <- rep(1:12, each = 20)   # ~1.3 s per month at 15 fps, ~16 s total

invisible(av::av_capture_graphics(
  for (m in frames) print(frame_plot(m)),
  output = "hawaii_cetacean_seasonal.mp4",
  width = 1400, height = 820, framerate = 15
))