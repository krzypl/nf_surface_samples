library(tidyverse)
library(sf)
library(terra)
library(tmap)

#data source: https://www.gov.nl.ca/eccc/natural-areas/gis-data/?utm_source=chatgpt.com

eco <- st_read("data/ecoregions/GNL_Ecoregions_2022_Newfoundland.shp")

cores_coordinates <- read_delim("data/cores_coordinates.csv",
                                delim = ";",
                                escape_double = FALSE, trim_ws = TRUE)

cores_coordinates_sf <- vect(cores_coordinates, geom = c("longitude_core", "latitude_core"), crs = "EPSG:4326")

cores_coordinates_utm_prep <- project(cores_coordinates_sf, "EPSG:26720")

cores_coordinates_utm <- st_as_sf(cores_coordinates_utm_prep) 

eco_utm <- st_transform(eco, 26720)

#eco_utm <- st_make_valid(regional_nf_utm)

buffers <- st_buffer(cores_coordinates_utm, 1000)

eco_extract <- st_intersection(buffers, eco_utm)

tmap_mode("view")
tm_shape(eco_extract) +
  tm_polygons("ECO_NAME", palette = "Set3", title = "ecoregions") +
  tm_layout(frame = FALSE) +
  tm_shape(cores_coordinates_utm) +
  tm_dots(size = 0.3) +
  tm_text("lakeID")

eco_summary <- eco_extract %>% 
  group_by(lakeID) %>% 
  summarise(ecoregion_types = paste(unique(ECO_NAME), collapse = ", "))

ecoregions_matrix <- eco_summary %>% 
  st_drop_geometry()  %>% 
  separate_rows(ecoregion_types, sep = ",\\s*")  %>% 
  mutate(presence = 1)  %>% 
  pivot_wider(
    names_from = ecoregion_types,
    values_from = presence,
    values_fill = 0
  )

write_csv(ecoregions_matrix, "data/ecoregions_matrix.csv")
