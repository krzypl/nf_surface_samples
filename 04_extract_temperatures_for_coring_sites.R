library(tidyverse)
library(terra)
library(tmap)
library(sf)

nf_polygon <- st_read("data/nf_polygon.gpkg")

#extract summer temperature ----

temp_clip_large_summer_mean <- rast("data/temp_clip_large_expand.tif")

temp_clip_large_june <- rast("data/temp_clip_large_expand_june.tif")

temp_clip_large_july <- rast("data/temp_clip_large_expand_july.tif")

temp_clip_large_august <- rast("data/temp_clip_large_expand_august.tif")

cores_coordinates <- read_delim("data/cores_coordinates.csv",
                                delim = ";",
                                escape_double = FALSE, trim_ws = TRUE)

cores_coordinates_sf <- vect(cores_coordinates, geom = c("longitude_core", "latitude_core"), crs = "EPSG:4326")

temperatures_for_coring_sites <- tibble(
  lakeID = cores_coordinates$lakeID,
  x = cores_coordinates$longitude_core,
  y = cores_coordinates$latitude_core,
  summer_mean_temperature = terra::extract(temp_clip_large_summer_mean,
                                    cores_coordinates_sf, ID = FALSE)$mean,
  june_temperature = terra::extract(temp_clip_large_june,
                             cores_coordinates_sf,
                             ID = FALSE)$wc2.1_30s_tavg_06,
  july_temperature = terra::extract(temp_clip_large_july,
                             cores_coordinates_sf,
                             ID = FALSE)$wc2.1_30s_tavg_07,
  august_temperature = terra::extract(temp_clip_large_august,
                             cores_coordinates_sf,
                             ID = FALSE)$wc2.1_30s_tavg_08
)

write_csv(temperatures_for_coring_sites,
          "data/temperature_for_coring_sites.csv")

temp_hists <- temperatures_for_coring_sites %>% 
  pivot_longer(cols = summer_mean_temperature:august_temperature, names_to = "variable", values_to = "value") %>% 
  ggplot() +
  geom_histogram(aes(value), bins = 50) +
  facet_wrap(.~ variable, ncol = 2, scales = "free_x")

tmap_mode("view")

cores_map <- tm_shape(nf_polygon) +
  tm_polygons() +
  tm_shape(temp_clip_large_summer_mean) +
  tm_raster(col.scale = tm_scale_continuous(values = terrain.colors(9))) +
  tm_layout(main.title = "Mean Summer Temperature", 
            legend.outside = TRUE) +
  tm_shape(cores_coordinates_sf) +
  tm_dots(size = 1) +
  tm_text("lakeID", size = 0.7, auto.placement = TRUE)

#extract winter temperature ----

temp_clip_large_winter_mean <- rast("data/winter_clip_large_expand.tif")

temp_clip_large_january <- rast("data/temp_clip_large_expand_january.tif")

temp_clip_large_february <- rast("data/temp_clip_large_expand_february.tif")

temp_clip_large_decemeber <- rast("data/temp_clip_large_expand_december.tif")

cores_coordinates <- read_delim("data/cores_coordinates.csv",
                                delim = ";",
                                escape_double = FALSE, trim_ws = TRUE)

cores_coordinates_sf <- vect(cores_coordinates, geom = c("longitude_core", "latitude_core"), crs = "EPSG:4326")

winter_temperatures_for_coring_sites <- tibble(
  lakeID = cores_coordinates$lakeID,
  x = cores_coordinates$longitude_core,
  y = cores_coordinates$latitude_core,
  winter_mean_temperature = terra::extract(temp_clip_large_winter_mean,
                                           cores_coordinates_sf, ID = FALSE)$mean,
  january_temperature = terra::extract(temp_clip_large_january,
                                    cores_coordinates_sf,
                                    ID = FALSE)$wc2.1_30s_tavg_01,
  february_temperature = terra::extract(temp_clip_large_february,
                                    cores_coordinates_sf,
                                    ID = FALSE)$wc2.1_30s_tavg_02,
  december_temperature = terra::extract(temp_clip_large_december,
                                      cores_coordinates_sf,
                                      ID = FALSE)$wc2.1_30s_tavg_12
)

write_csv(winter_temperatures_for_coring_sites,
          "data/winter_temperature_for_coring_sites.csv")

winter_temp_hists <- winter_temperatures_for_coring_sites %>% 
  pivot_longer(cols = winter_mean_temperature:december_temperature, names_to = "variable", values_to = "value") %>% 
  ggplot() +
  geom_histogram(aes(value), bins = 50) +
  facet_wrap(.~ variable, ncol = 2, scales = "free_x")
