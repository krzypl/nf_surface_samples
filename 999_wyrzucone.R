#do plotu w quatro pokazujacego dane klimatyczne:
```{r}
summer_temperature <- read_csv("data/temperature_for_coring_sites.csv")

winter_temperature <- read_csv("data/winter_temperature_for_coring_sites.csv")

summer_precipitation <- read_csv("data/precipitation_for_coring_sites.csv")

winter_precipitation <- 
  read_csv("data/winter_precipitation_for_coring_sites.csv")

list_climate <- list(summer_temperature,
                     winter_temperature,
                     summer_precipitation,
                     winter_precipitation)

climate_data <- reduce(list_climate, left_join,
                       by = "lakeID") %>%
  select(lakeID, summer_mean_temperature, winter_mean_temperature,
         summer_sum_precipitation, winter_sum_precipitation) %>% 
  pivot_longer(cols = summer_mean_temperature:winter_sum_precipitation,
               names_to = "variable", values_to = "value")

distribution_of_climate_vars <- ggplot(climate_data) +
  geom_histogram(aes(x = value)) +
  facet_wrap(.~ variable, scales = "free")

distribution_of_climate_vars

```


#stare site selection ----
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

canada <- ne_states(country = "Canada", returnclass = "sf")

newfoundland <- canada[canada$name == "Newfoundland and Labrador", ]

nfl_parts <- st_cast(newfoundland, "POLYGON")

nfl_parts$area <- st_area(nfl_parts)

nfl_parts <- nfl_parts %>% arrange(desc(area))

nf_polygon <- nfl_parts[2,] %>% select(geometry)

st_write(nf_polygon, "data/nf_polygon.gpkg")

nf_proj <- st_transform(nf_polygon, crs = 26922) 

area <- st_area(nf_proj)

cell_area <- as.numeric(area) / 80
cell_size <- sqrt(cell_area)  # This gives you cell size in meters

grid <- st_make_grid(nf_proj, cellsize = cell_size, square = TRUE)

grid <- st_as_sf(grid)

grid_clipped <- st_intersection(grid, nf_proj)

grid_clipped <- st_transform(grid_clipped, crs = 4326)

n_grids <- nrow(grid_clipped)

centroids <- st_centroid(grid_clipped)

centroids <- st_transform(centroids, crs = 4326)

centroid_coords <- st_coordinates(centroids)

randys_sites <- read_csv("data/randys_sites_coords.csv")

randys_sites_sf <- st_as_sf(randys_sites, coords = c("x", "y"), crs = 4326)

plot(st_geometry(grid_clipped), border = 'blue')
plot(centroids, add = TRUE, col = 'black', pch = 16)
plot(randys_sites_sf, add = TRUE, col = "red", pch = 18)

#check for historical climate data: https://github.com/paleolimbot/rclimateca, or obtain gridded data from https://worldclim.org/ and then plot a map


library(leaflet)
library(leafem)


leaflet() %>%
  addProviderTiles("OpenStreetMap") %>%
  addRasterImage(temp_clip, 
                 opacity = 0.6)

pal <- colorNumeric("Spectral", values(temp_clip), na.color = NA)

r_proj <- project(temp_clip, "EPSG:3857")

leaflet() %>%
  addProviderTiles(providers$OpenTopoMap, group = "Topo") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addRasterImage(r_proj, opacity = 0.9, layerId = "raster") %>%
  addLayersControl(
    baseGroups = c("Topo", "Satellite"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  addMouseCoordinates() %>%
  addImageQuery(r_proj, project = TRUE)
