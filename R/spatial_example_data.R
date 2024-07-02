# Script to generate the aggregation.rds data object

# Process Scotland intermediate zone shapefiles
# for the aggregation to areal units example

# This script does the following:
# 1) Downloads intermediate zone data + lookup tables from official
# government websites
# 2) Selects all intermediate zones in Glasgow City
# 3) Changes map units to km
# 4) Creates a Glasgow city boundary polygon
# 5) Creates an SPDE mesh and integration matrix for
# the data aggregation example
# 6) Saves the objects

# Updated to use sf instead of sp objects.

library(dplyr)
library(inlabru)
library(fmesher)
library(INLA)
library(ggplot2)
library(sf)

# Data download from:
# https://statistics.gov.scot/data/data-zone-lookup
# https://data.gov.uk/dataset/133d4983-c57d-4ded-bc59-390c962ea280/intermediate-zone-boundaries-2011

Lookup_URL <- "https://statistics.gov.scot/downloads/file?id=1ab6565e-10e0-4888-b91c-4ae6821b30d7%2FDatazone2011lookup+%288%29.csv"
Lookup_path <- here::here("Data", "raw", "DZ_lookup.csv")

IZ_URL <- "https://maps.gov.scot/ATOM/shapefiles/SG_IntermediateZoneBdry_2011.zip"
IZ_zip_path <- here::here("Data", "raw", "SG_IntermediateZone_Bdry_2011.zip")
IZ_path <- here::here("Data", "raw", "SG_IntermediateZone_Bdry_2011.shp")

if (!file.exists(IZ_zip_path)) {
  if (!dir.exists(dirname(IZ_zip_path))) {
    dir.create(dirname(IZ_zip_path), recursive = TRUE)
  }
  options(timeout = 600)
  download.file(
    url = IZ_URL,
    destfile = IZ_zip_path,
    method = "wget",
    extra = "--no-check-certificate"
  )
}
utils::unzip(IZ_zip_path,
  exdir = dirname(IZ_zip_path)
)
# file.remove(IZ_zip_path)

# lookup file to join local authority name
if (!file.exists(Lookup_path)) {
  if (!dir.exists(dirname(Lookup_path))) {
    dir.create(dirname(Lookup_path), recursive = TRUE)
  }
  download.file(
    url = Lookup_URL,
    destfile = Lookup_path
  )
}

shp <- st_read(IZ_path)
lookup <- read.csv(Lookup_path)
iz_la_lookup <- lookup %>%
  select(IZ2011_Code, LA_Code, LA_Name) %>%
  unique()

# Convert map units to km for smaller numbers
st_crs(shp)$units
st_crs(shp)$proj4string
new_crs <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"
shp2 <- st_transform(shp, st_crs(new_crs))
st_crs(shp2)$units
st_crs(shp2)$proj4string

# select Glasgow City intermediate zones
shp3 <- shp2 %>%
  mutate(ID = 1:nrow(shp)) %>%
  rename(IZ2011_Code = InterZone) %>%
  left_join(iz_la_lookup, by = "IZ2011_Code") %>%
  filter(LA_Name == "Glasgow City")

glasgow_iz <- shp3

# Glasgow boundary
bnd_sfc <- st_union(glasgow_iz)
bnd <- st_cast(st_sf(data.frame(geometry = bnd_sfc)), "POLYGON")
bnd_outer <- st_cast(
  st_sf(geometry = fm_nonconvex_hull(bnd, 5, dTolerance = 5 / 40 * 4)),
  "POLYGON"
)

boundary <- fm_as_segm_list(list(bnd, bnd_outer))

## Build the mesh:
mesh <- fm_mesh_2d_inla(
  boundary = boundary,
  max.edge = c(0.3, 4),
  min.angle = c(30, 21),
  max.n = c(48000, 16000), ## Safeguard against large meshes.
  max.n.strict = c(128000, 128000), ## Don't build a huge mesh!
  cutoff = 0.3, ## Filter away adjacent points.
  crs = fm_crs(bnd)
) ## Offset for extra boundaries, if needed.

glasgow_iz$ID <- 1:nrow(glasgow_iz)
ips <- fm_int(mesh, samplers = glasgow_iz)
ips$ID <- ips$.block

out <- list(
  glasgow_iz = glasgow_iz,
  glasgow_bnd = bnd,
  mesh = mesh,
  ips = ips
)

saveRDS(out, here::here("Data", "aggregation_data.rds"))

# delete intermediate files?
if (FALSE) {
  to_delete <- list.files(here::here("Data", "raw"),
    pattern = "SG_IntermediateZone_Bdry_2011",
    full.names = TRUE
  )
  to_delete <- append(
    to_delete,
    list(Lookup_path)
  )

  do.call(file.remove, to_delete)
}
