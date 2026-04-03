library(tidyverse)
library(terra)
library(tmap)
library(sf)
library(mapview)

#display a map ----

temp_clip <- rast("data/temp_clip.tif")

temp_clip_large <- rast("data/temp_clip_large_expand.tif")

nf_polygon <- st_read("data/nf_polygon.gpkg")

nf_polygon_large<- st_read("data/nf_polygon_large.gpkg")

temp_clip_x <- ifel(temp_clip_large >= 12.7 & temp_clip_large < 13.5,
                    temp_clip_large,      
                    NA)

summer_mean_x <- tm_shape(nf_polygon) +
  tm_polygons() +
  tm_shape(temp_clip_x) +
  tm_raster(palette = "-spectral", title = "10.3 - 10.5") +
  tm_layout(main.title = "Mean Summer Temperature", 
            legend.outside = TRUE)

mapview(temp_clip, 
        alpha = 0.4,
        na.color = NA)

mapview(temp_clip_x, 
        alpha = 0.4,
        na.color = NA)

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

#extract summer precipitation ----

prec_clip_large_summer_sum <- rast("data/prec_clip_large_expand_prec.tif")

prec_clip_large_june <- rast("data/prec_clip_large_expand_june.tif")

prec_clip_large_july <- rast("data/prec_clip_large_expand_july.tif")

prec_clip_large_august <- rast("data/prec_clip_large_expand_august.tif")

cores_coordinates <- read_delim("data/cores_coordinates.csv",
                                delim = ";",
                                escape_double = FALSE, trim_ws = TRUE)

cores_coordinates_sf <- vect(cores_coordinates, geom = c("longitude_core", "latitude_core"), crs = "EPSG:4326")

precipitation_for_coring_sites <- tibble(
  lakeID = cores_coordinates$lakeID,
  x = cores_coordinates$longitude_core,
  y = cores_coordinates$latitude_core,
  summer_sum_precipitation = terra::extract(prec_clip_large_summer_sum,
                                           cores_coordinates_sf, ID = FALSE)$sum,
  june_precipitation = terra::extract(prec_clip_large_june,
                                    cores_coordinates_sf,
                                    ID = FALSE)$wc2.1_2.5m_prec_06,
  july_precipitation = terra::extract(prec_clip_large_july,
                                    cores_coordinates_sf,
                                    ID = FALSE)$wc2.1_2.5m_prec_07,
  august_precipitation = terra::extract(prec_clip_large_august,
                                      cores_coordinates_sf,
                                      ID = FALSE)$wc2.1_2.5m_prec_08
)

write_csv(precipitation_for_coring_sites,
          "data/precipitation_for_coring_sites.csv")

##extract winter precipitation ----

prec_clip_large_winter_sum <- rast("data/winter_clip_large_expand_prec.tif")

prec_clip_large_january <- rast("data/prec_clip_large_expand_january.tif")

prec_clip_large_february <- rast("data/prec_clip_large_expand_february.tif")

prec_clip_large_december <- rast("data/prec_clip_large_expand_december.tif")

cores_coordinates <- read_delim("data/cores_coordinates.csv",
                                delim = ";",
                                escape_double = FALSE, trim_ws = TRUE)

cores_coordinates_sf <- vect(cores_coordinates, geom = c("longitude_core", "latitude_core"), crs = "EPSG:4326")

winter_precipitation_for_coring_sites <- tibble(
  lakeID = cores_coordinates$lakeID,
  x = cores_coordinates$longitude_core,
  y = cores_coordinates$latitude_core,
  winter_sum_precipitation = terra::extract(prec_clip_large_winter_sum,
                                            cores_coordinates_sf, ID = FALSE)$sum,
  january_precipitation = terra::extract(prec_clip_large_january,
                                      cores_coordinates_sf,
                                      ID = FALSE)$wc2.1_2.5m_prec_01,
  february_precipitation = terra::extract(prec_clip_large_february,
                                      cores_coordinates_sf,
                                      ID = FALSE)$wc2.1_2.5m_prec_02,
  december_precipitation = terra::extract(prec_clip_large_december,
                                        cores_coordinates_sf,
                                        ID = FALSE)$wc2.1_2.5m_prec_12
)

write_csv(winter_precipitation_for_coring_sites,
          "data/winter_precipitation_for_coring_sites.csv")


