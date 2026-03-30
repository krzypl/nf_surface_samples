library(ClimDatDownloadR)
library(terra)
library(sf)
library(tmap)

# WorldClim.HistClim.download(
#   save.location = "./data",
#   parameter = c("temp"),
#   month.var = c(6:8),
#   resolution = "2.5m",
#   version.var = "2.1",
#   clipping = TRUE,
#   clip.extent = c(-60,-52, 46, 52),
#   convert.files.to.asc = FALSE,
#   stacking.data = FALSE,
#   keep.raw.zip = FALSE,
#   delete.raw.data = FALSE,
#   save.bib.file = TRUE
# )

# files <- list.files("data/temp/WorldClim_2.1_tavg_2.5m/clipped_2025-06-18_13-27-34", full.names = TRUE)

files <- list.files("data/map_raw", full.names = TRUE)

summer_temp <- vector("list", length = 12)

for (i in seq_along(summer_temp)) {
  summer_temp[[i]] <- rast(files[i])
}

summer_temp <- summer_temp[6:8]


#ponizej wszystko dziala, ale zajmuje duzo czasu, wiec zapisalem plik i go odczytuje
summer_stack <- rast(summer_temp)

summer_mean <- mean(summer_stack)

#plot(summer_mean, main = "Mean Temperature Map")

nf_polygon <- st_read("data/nf_polygon.gpkg")
nf_polygon_large <- st_read("data/nf_polygon_large.gpkg")

#temp_clip <- mask(crop(summer_mean, vect(nf_polygon)), vect(nf_polygon), touches = FALSE)
temp_clip_large <- mask(crop(summer_mean, vect(nf_polygon_large)), vect(nf_polygon_large), touches = TRUE)

local_mean <- focal(temp_clip_large, w = 3, fun = mean, na.rm = TRUE)

temp_clip_large_expand <- cover(temp_clip_large, local_mean) #this is to give values for cells in coastal areas that have NAs. The values are means from 3x3 matrix; the non-NA values remains unchanged

raster::writeRaster(temp_clip_large_expand, 
                    filename = "data/temp_clip_large_expand.tif",
                    overwrite = TRUE)

temp_clip_large_june <- mask(crop(summer_temp[[1]], vect(nf_polygon_large)), vect(nf_polygon_large), touches = TRUE)

local_mean_june <- focal(temp_clip_large_june, w = 3, fun = mean, na.rm = TRUE)

temp_clip_large_expand_june <- cover(temp_clip_large_june, local_mean_june) 

raster::writeRaster(temp_clip_large_expand_june, 
                    filename = "data/temp_clip_large_expand_june.tif",
                    overwrite = TRUE)

temp_clip_large_july <- mask(crop(summer_temp[[2]], vect(nf_polygon_large)), vect(nf_polygon_large), touches = TRUE)

local_mean_july <- focal(temp_clip_large_july, w = 3, fun = mean, na.rm = TRUE)

temp_clip_large_expand_july <- cover(temp_clip_large_july, local_mean_july) 

raster::writeRaster(temp_clip_large_expand_july, 
                    filename = "data/temp_clip_large_expand_july.tif",
                    overwrite = TRUE)

temp_clip_large_august <- mask(crop(summer_temp[[3]], vect(nf_polygon_large)), vect(nf_polygon_large), touches = TRUE)

local_mean_august <- focal(temp_clip_large_august, w = 3, fun = mean, na.rm = TRUE)

temp_clip_large_expand_august <- cover(temp_clip_large_august, local_mean_august) 

raster::writeRaster(temp_clip_large_expand_august, 
                    filename = "data/temp_clip_large_expand_august.tif",
                    overwrite = TRUE)

#raster::writeRaster(temp_clip, filename = "data/temp_clip.tif", overwrite = TRUE)

temp_clip <- rast("data/temp_clip.tif")

temp_clip_large <- rast("data/temp_clip_large.tif")

#temperatures in the same scale -----
#mean summer temp ----

summaer_mean_map <- tm_shape(temp_clip) +
  tm_raster(palette = "-spectral", title = "Min. temp = 10.3, max. temp = 15.4 (°C)",
            breaks = seq(6.4, 17.1, length.out = 35)) +
  tm_layout(main.title = "Mean Summer Temperature", 
            legend.outside = TRUE)

#june temp ----

temp_clip_june <- mask(crop(summer_temp[[1]], vect(nf_polygon)), vect(nf_polygon), touches = FALSE)

june_temp_map <- tm_shape(temp_clip_june) +
  tm_raster(palette = "-spectral", title = "Min. temp = 6.4, max. temp = 12.7 (°C)",
            breaks = seq(6.4, 17.1, length.out = 35)) +
  tm_layout(main.title = "June Temperature", 
            legend.outside = TRUE)


#July temperature ----

temp_clip_july <- mask(crop(summer_temp[[2]], vect(nf_polygon)), vect(nf_polygon), touches = FALSE)

july_temp_map <- tm_shape(temp_clip_july) +
  tm_raster(palette = "-spectral", title = "Min. temp = 11.3, max. temp = 17.1 (°C)",
            breaks = seq(6.4, 17.1, length.out = 35)) +
  tm_layout(main.title = "July Temperature", 
            legend.outside = TRUE)

#August temperature ----

temp_clip_august <- mask(crop(summer_temp[[3]], vect(nf_polygon)), vect(nf_polygon), touches = FALSE)

august_temp_map <- tm_shape(temp_clip_august) +
  tm_raster(palette = "-spectral", title = "Min. temp = 12.7, max. temp = 16.9 (°C)",
            breaks = seq(6.4, 17.1, length.out = 35)) +
  tm_layout(main.title = "August Temperature", 
            legend.outside = TRUE)

#combine maps ----

combined_map <- tmap_arrange(summaer_mean_map, june_temp_map, july_temp_map, august_temp_map, ncol = 2, nrow = 2)

tmap_save(combined_map,
          filename = "figures/combined_map.png", width = 10, height = 8)

#temperatures with season-addjusted temperature range ----

#mean summer temp ----

summaer_mean_map_a <- tm_shape(temp_clip) +
  tm_raster(palette = "-spectral", title = "Min. temp = 10.3, max. temp = 15.4 (°C)",
            breaks = seq(10.3, 15.4, length.out = 35)) +
  tm_layout(main.title = "Mean Summer Temperature", 
            legend.outside = TRUE)

#june temp ----

temp_clip_june <- mask(crop(summer_temp[[1]], vect(nf_polygon)), vect(nf_polygon), touches = FALSE)

june_temp_map_a <- tm_shape(temp_clip_june) +
  tm_raster(palette = "-spectral", title = "Min. temp = 6.4, max. temp = 12.7 (°C)",
            breaks = seq(6.4, 12.7, length.out = 35)) +
  tm_layout(main.title = "June Temperature", 
            legend.outside = TRUE)


#July temperature ----

temp_clip_july <- mask(crop(summer_temp[[2]], vect(nf_polygon)), vect(nf_polygon), touches = FALSE)

july_temp_map_a <- tm_shape(temp_clip_july) +
  tm_raster(palette = "-spectral", title = "Min. temp = 11.3, max. temp = 17.1 (°C)",
            breaks = seq(11.3, 17.1, length.out = 35)) +
  tm_layout(main.title = "July Temperature", 
            legend.outside = TRUE)

#August temperature ----

temp_clip_august <- mask(crop(summer_temp[[3]], vect(nf_polygon)), vect(nf_polygon), touches = FALSE)

august_temp_map_a <- tm_shape(temp_clip_august) +
  tm_raster(palette = "-spectral", title = "Min. temp = 12.7, max. temp = 16.9 (°C)",
            breaks = seq(12.7, 16.9, length.out = 35)) +
  tm_layout(main.title = "August Temperature", 
            legend.outside = TRUE)

#combine maps ----

combined_map_a <- tmap_arrange(summaer_mean_map_a, june_temp_map_a, july_temp_map_a, august_temp_map_a, ncol = 2, nrow = 2)

tmap_save(combined_map_a,
          filename = "figures/combined_map_season_a.png", width = 10, height = 8)
