library(terra)
library(sf)
library(tmap)
library(mapview)

#read data ----

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
