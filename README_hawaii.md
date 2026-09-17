# hawaii-cetacean-seasonality

Seasonal distribution of systematic cetacean sighting records in the Hawaiian
Islands EEZ, from NOAA PIFSC shipboard surveys, 2006–2017.

![Final frame](outputs/figures/still_september.png)

A twelve-frame animation, January to December, of every sighting in the public
WinCruz record that could support density estimation. Bathymetry underneath,
EEZ boundary outlined, the species seen that month shown down the right.

Full animation: [`outputs/hawaii_cetacean_seasonal.mp4`](outputs/hawaii_cetacean_seasonal.mp4)

## What the map shows, and what it does not

It shows **sighting records**, not animals. The public dataset contains
sightings only — there is no trackline or effort file — so the map cannot say
where cetaceans are, only where and when a sighting was written down.

Blank months are **months with no systematic line-transect effort**, which is
not the same as months with no survey activity, and certainly not months with
no animals. March and April both carry non-systematic records; December carries
none at all.

## From 1,492 records to 361 usable sightings

| step | n |
|---|---|
| all records | 1,492 |
| cetaceans (excludes turtles, monk seal) | 1,491 |
| identified to species | 1,042 |
| Hawaiian Islands EEZ | 882 |
| systematic effort | 370 |
| with a group-size estimate | **361** |

Three quarters of the file falls away before the data can carry a density
estimate. The largest remaining sample is **sperm whale at 56 sightings**,
across 24 species and twelve years. The conventional minimum for fitting a
species-specific detection function is 60–80 detections, so on this public
record not one species clears the bar alone — which is why detection functions
in practice are pooled across species by size class or guild.

## The trap in `EffortType`

From the InPort entity metadata:

> Did the sighting occur when the survey effort was systematic (**S**),
> non-systematic (**N**), fine-scale (**F**), or off (**O**)

`O` is **off**-effort. It reads like "on". Only `S` is line-transect effort and
only `S` belongs in a density analysis. Filtering on `O` yields 355 sightings
that look plausible and are the wrong ones.

## Seasonal coverage

Systematic sightings in the Hawaiian Islands EEZ, by month, 2006–2017 pooled:

| Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | Sep | Oct | Nov | Dec |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 46 | 0 | 0 | 57 | 6 | 35 | 34 | **114** | 46 | 22 | 0 |

September alone holds 31% of the record. July–November holds 70%. The winter
record is essentially one February — which is worth knowing, given that the
existing PIFSC humpback density model for this EEZ was built from a February
2009 survey, and that a new winter survey is the basis of current density
work.

`data/derived/effort_by_month.csv` gives the same breakdown across all four
effort types, before the species and group-size filters.

## Reproducing

```r
install.packages(c("dplyr", "readr", "sf", "here", "ggplot2", "terra",
                   "marmap", "rphylopic", "av", "purrr", "tidyr",
                   "stringr", "lubridate", "scales"))

source("R/01_prepare_data.R")   # downloads, filters, writes data/derived/
source("R/02_seasonal_map.R")   # renders outputs/hawaii_cetacean_seasonal.mp4
```

The sighting archive downloads automatically. The EEZ shapefile does not —
Marine Regions asks that you obtain it from them rather than take a
redistributed copy. Download **EEZ v12** from
[marineregions.org/downloads.php](https://www.marineregions.org/downloads.php)
and unpack it to `data/raw/World_EEZ_v12_20231025/`.

### Two things that will bite you

The Hawaiian Islands EEZ **crosses the antimeridian**, so its bounding box in
plain −180/180 spans the entire globe and eight sightings sit east of the line.
Everything here works in 0–360. The source polygon is also split at 180°, and
that seam survives `st_shift_longitude()` as an interior line that
`geom_path()` will happily draw across the map; `01` strips it.

`marmap::getNOAA.bathy(antimeridian = TRUE)` stitches two separately fetched
halves with **different cell sizes** (0.06694° and 0.06793°, plus a 0.035°
sliver at the join). `geom_raster()` assumes a regular grid, so the mismatch
renders as vertical striping across the whole panel. `02` resamples onto a
single regular 0.07° grid, and asserts regularity afterwards.

## Attribution

**Sightings** — NOAA Pacific Islands Fisheries Science Center, Cetacean
Research Program. *Shipboard Cetacean Surveys — Visual Surveys — WinCruz
Sighting Records*, [InPort 18141](https://www.fisheries.noaa.gov/inport/item/18141).
Collected under NMFS permits **PIFSC 15240, PIFSC 20311, SWFSC 14097,
SWFSC 774-1714**. Citing the permit numbers is a condition of use, not a
courtesy.

**Bathymetry** — NOAA ETOPO, retrieved with
[marmap](https://cran.r-project.org/package=marmap).

**EEZ boundary** — Flanders Marine Institute, *Maritime Boundaries Geodatabase*
v12 (2023), [marineregions.org](https://www.marineregions.org/), CC BY 4.0.

**Silhouettes** — [PhyloPic](https://www.phylopic.org/), retrieved with
[rphylopic](https://cran.r-project.org/package=rphylopic). Per-image
contributors and licences in [`outputs/phylopic_credits.csv`](outputs/phylopic_credits.csv).
Three species have no PhyloPic image and use a genus- or family-level stand-in:
pantropical spotted dolphin (*Stenella*), Bryde's whale (*Balaenoptera*),
Longman's beaked whale (Ziphiidae).
