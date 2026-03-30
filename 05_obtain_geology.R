#remotes::install_github("yonghah/esri2sf")
library(esri2sf)
library(sf)
library(tidyverse)
library(terra)
library(tmap)

#local geology ----

#to byloby dobre, ale jest niewielkie pokrycie, dla wielu obszarow nie ma danych, wiec wykorzystuje nastepny podklad, czyli regional geology.

# layer_url <- "https://dnrmaps.gov.nl.ca/arcgis/rest/services/GeoAtlas/Surficial_Geology_All/MapServer/11"
# 
 nf_bbox_wgs <- c(xmin = -59.0, ymin = 46.5, xmax = -52.0, ymax = 51.0)
# 
 bbox_coords_matrix <- matrix(
   c(
     nf_bbox_wgs["xmin"], nf_bbox_wgs["ymin"],
     nf_bbox_wgs["xmin"], nf_bbox_wgs["ymax"],
     nf_bbox_wgs["xmax"], nf_bbox_wgs["ymax"],
     nf_bbox_wgs["xmax"], nf_bbox_wgs["ymin"],
     nf_bbox_wgs["xmin"], nf_bbox_wgs["ymin"]  
   ), 
   ncol = 2, byrow = TRUE
 )
# 
 bbox_poly <- st_sfc(st_polygon(list(bbox_coords_matrix)), crs = 4326)
# 
 bbox_proj <- st_transform(bbox_poly, 26720)
 bbox_obj <- st_bbox(bbox_proj)

# layer_url <- "https://dnrmaps.gov.nl.ca/arcgis/rest/services/GeoAtlas/Surficial_Geology_All/MapServer/11"
# 
 # surficial_nf <- esri2sf::esri2sf(
 #   url = layer_url,
 #   bbox = bbox_obj,
 #   crs = 26720  # CRS serwisu
 # )
# surficial_nf_wgs <- st_transform(surficial_nf, 4326)

#st_write(surficial_nf_wgs, "data/geology/surficial_geology.gpkg", delete_dsn = TRUE)

surficial_nf_wgs <- st_read("data/geology/surficial_geology.gpkg")

#regional geology ----

#wykorzystuje to zamiast surficial, bo ma lepsze pokrycie, chociaz tez nie jest idealnie i nie wszystkie jeziora maja znana geolgie. postarac sie im przypisac cos na podstawie obserwacji terenowej plus na podstawie innych danych (np geologia w skali 1M)

# layer_url_regional <- "https://dnrmaps.gov.nl.ca/arcgis/rest/services/GeoAtlas/Surficial_Geology_All/MapServer/12"
# 
# 
# regional_nf <- esri2sf::esri2sf(
#   url = layer_url_regional,
#   bbox = bbox_obj,
#   crs = 26720  # CRS serwisu
# )
# 
# regional_nf_wgs <- st_transform(regional_nf, 4326)
# 
# st_write(regional_nf_wgs, "data//geology/regional_geology.gpkg", delete_dsn = TRUE)

regional_nf_wgs <- st_read("data/geology/regional_geology.gpkg")

tmap_mode("view")
tm_shape(regional_nf_wgs) +
  tm_polygons("GENETIC250", palette = "Set3", title = "Surficial Geology") +
  tm_layout(frame = FALSE)

#bedrock geology ----

#to wykorzystuje do warstwy bedrock. tutaj mam dokladniejszcze oznaczenie co to za skaly, wiec to bedzie wsparcie dla tamtej geologii powierzchniowej. Nie jest pełne pokrycie

# layer_url_bedrock_geology <- "https://dnrmaps.gov.nl.ca/arcgis/rest/services/GeoAtlas/Bedrock_Geology_All/MapServer/19"
# 
# detailed_geology_nf <- esri2sf::esri2sf(
#   url = layer_url_bedrock_geology,
#   bbox = bbox_obj,
#   crs = 26720  # CRS serwisu
# )
# 
# bedrock_geology_nf_wgs <- st_transform(detailed_geology_nf, 4326)
# 
# st_write(bedrock_geology_nf_wgs, "data/geology/bedrock_geology.gpkg", delete_dsn = TRUE)

bedrock_geology_nf_wgs <- st_read("data/geology/bedrock_geology.gpkg")

#extract regional geology for coring sites ----

cores_coordinates <- read_delim("data/cores_coordinates.csv",
                                delim = ";",
                                escape_double = FALSE, trim_ws = TRUE)

cores_coordinates_sf <- vect(cores_coordinates, geom = c("longitude_core", "latitude_core"), crs = "EPSG:4326")

cores_coordinates_utm_prep <- project(cores_coordinates_sf, "EPSG:26720")

cores_coordinates_utm <- st_as_sf(cores_coordinates_utm_prep) 

regional_nf_utm <- st_transform(regional_nf_wgs, 26720)

regional_nf_utm <- st_make_valid(regional_nf_utm)

buffers <- st_buffer(cores_coordinates_utm, 1000)

regional_extract <- st_intersection(buffers, regional_nf_utm)

tmap_mode("view")
tm_shape(regional_extract) +
  tm_polygons("GENETIC250", palette = "Set3", title = "regional Geology") +
  tm_layout(frame = FALSE) +
  tm_shape(cores_coordinates_utm) +
  tm_dots(size = 0.3)

regional_geology_summary <- regional_extract %>% 
  group_by(lakeID) %>% 
  summarise(geology_type = paste(unique(GENETIC250), collapse = ", "))

lakes_missing_regional <- cores_coordinates %>% 
  filter(!lakeID %in% unique(regional_geology_summary$lakeID)) %>% 
  pull(lakeID) #data for these three lakes is missing, but they are all located in the area that can be easily classified as exposed bedrock

lakes_add <- tibble(
  lakeID = lakes_missing_regional,
  geology_type = "exposed bedrock"
)

regional_geology_summary_final <- st_drop_geometry(regional_geology_summary) %>%
  add_row(lakes_add)

#surface geology final ----
#classification scheme:
# 1 Exposed bedrock
# Indicates bedrock is visible or can be easily exposed.
# 
# Keywords:
#   
#   exposed bedrock
# till veneer (thin discontinuous cover)
# 
# 2 Sediment cover
# Bedrock is covered by sediment.
# 
# Keywords:
#   
#   concealed bedrock
# till blanket
# ridged till
# hummocky terrain
# glaciofluvial gravel and sand
# marine clay
# sand, gravel and diamicton
# alluvium
# colluvium
# 
# 3 Bog
# Peatland/wetland surfaces.
# 
# Keyword:
#   
#   bog

classify_surface <- function(x){
  
  x <- tolower(x)
  
  exposed <- str_detect(x, "exposed bedrock|till veneer")
  
  sediment <- str_detect(x,
                         "concealed bedrock|till blanket|ridged till|hummocky terrain|glaciofluvial gravel and sand|marine clay|sand, gravel and diamicton|alluvium|colluvium")
  
  bog <- str_detect(x, "bog")
  
  tibble(
    exposed_bedrock = as.integer(exposed),
    covered_bedrock = as.integer(sediment),
    bog = as.integer(bog)
  )
  
}

surface_geology_matrix <- regional_geology_summary_final %>%
  bind_cols(classify_surface(.$geology_type))

surface_long <- surface_geology_matrix %>% select(!geology_type) %>% 
       pivot_longer(cols = exposed_bedrock:bog, names_to = "surface", values_to = "value")

surface_counts <- surface_long %>%
  group_by(surface) %>%
  summarize(count = sum(value))

ggplot(surface_counts, aes(x = surface, y = count)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(x = "Surface type", y = "Number of occurances", title = "Surface geology counts")

write_csv(surface_geology_matrix, "data/surface_geology_matrix.csv")
  

#extract bedrock geology for coring sites ----

bedrock_geology_nf_utm <- st_transform(bedrock_geology_nf_wgs, 26720)

bedrock_geology_nf_utm <- st_make_valid(bedrock_geology_nf_utm)

bedrock_geology_extract <- st_intersection(buffers, bedrock_geology_nf_utm) %>% 
  mutate(geology = ifelse(D_ROCKTYPE == " ", G_ROCKTYPE, D_ROCKTYPE))
   
tmap_mode("view")
tm_shape(bedrock_geology_extract) +
  tm_polygons("geology", palette = "Set3", title = "bedrock geology") +
  tm_layout(frame = FALSE) +
  tm_shape(cores_coordinates_utm) +
  tm_dots(size = 0.3)

bedrock_geology_summary <- bedrock_geology_extract %>% 
  group_by(lakeID) %>% 
  summarise(geology_type = paste(unique(geology), collapse = ", "))

#add missing data from lower resolution maps ----

bedrock_missing_nf_wgs <- st_read("data/geology/bedrock_1M.gpkg") #mialem problem, zeby sciagnac przez R, wiec pobralem przez QGISa; problem byl z ukladem wspolrzednych. Po dlugich staraniach udalo sie to ogarnac.

bedrock_missing_nf_utm <- st_transform(bedrock_missing_nf_wgs, 26720)

bedrock_missing_nf_utm <- st_cast(bedrock_missing_nf_utm, "MULTIPOLYGON")

bedrock_missing_nf_utm <- st_make_valid(bedrock_missing_nf_utm)

bedrock_missing_extract <- st_intersection(buffers, bedrock_missing_nf_utm)

tmap_mode("view")
tm_shape(bedrock_missing_extract) +
  tm_polygons("LITHOLOGY", palette = "brewer.set3", title = "bedrock geology") +
  tm_layout(frame = FALSE) +
  tm_shape(cores_coordinates_utm) +
  tm_dots(size = 0.3)

lakes2stay_bedrock_missing <- cores_coordinates %>% 
  filter(!lakeID %in% unique(bedrock_geology_summary$lakeID)) %>% 
  pull(lakeID)

bedrock_missing_geology_summary <- bedrock_missing_extract %>% 
  filter(lakeID %in% lakes2stay_bedrock_missing) %>% 
  group_by(lakeID) %>% 
  summarise(geology_type = paste(unique(LITHOLOGY), collapse = ", "))

bedrock_merged_summary <- bedrock_geology_summary %>%
  add_row(bedrock_missing_geology_summary)

#bedrock final ----
# Classification Rules in the Code
# Carbonate → "carbonate", "limestone", "dolostone", "marble", "breccia"
# Igneous Felsic → "plutonic felsic", "volcanic felsic"
# Igneous Mafic → "plutonic mafic", "hypabyssal mafic", "volcanic mafic", "ophiolite"
# Igneous Undifferentiated → general "plutonic", "volcanic", "submarine mafic to felsic volcanics"
# Siliciclastic → "siliciclastic" (sandstone, siltstone, shale, psammite, phyllite, conglomerate, etc.)
# Others → everything else (sedimentary, schist, gneiss, migmatite, amphibolite, melange, chert)


classify_rocks <- function(lithology){
  
  lithology <- tolower(lithology)
  
  carbonate <- str_detect(lithology,
                          "carbonate|limestone|dolostone|marble|carbonate breccia")
  
  felsic <- str_detect(lithology,
                       "plutonic felsic|volcanic felsic")
  
  mafic <- str_detect(lithology,
                      "plutonic mafic|volcanic mafic|hypabyssal mafic|ophiolite|ultramafic")
  
  igneous_general <- str_detect(lithology,
                                "\\bplutonic\\b|\\bvolcanic\\b|submarine mafic to felsic volcanics")
  
  igneous_undiff <- igneous_general & !felsic & !mafic
  
  siliciclastic <- str_detect(lithology,
                              "siliciclastic|sandstone|siltstone|shale|conglomerate|turbidite|psammite|greywacke|graywacke|phyllite|mudstone")
  
  others <- str_detect(lithology,
                       "sedimentary|schist|gneiss|migmatite|amphibolite|melange|chert")
  
  tibble(
    carbonate = as.integer(carbonate),
    igneous_felsic = as.integer(felsic),
    igneous_mafic = as.integer(mafic),
    igneous_undifferentiated = as.integer(igneous_undiff),
    siliciclastic = as.integer(siliciclastic),
    others = as.integer(others)
  )
}

bedrock_geology_matrix <- bedrock_merged_summary %>%
  bind_cols(classify_rocks(.$geology_type))

bedrock_long <- bedrock_geology_matrix %>% select(!geology_type) %>% 
  pivot_longer(cols = carbonate:others, names_to = "bedrock", values_to = "value")

bedrock_counts <- bedrock_long %>%
  group_by(bedrock) %>%
  summarize(count = sum(value))

ggplot(bedrock_counts, aes(x = bedrock, y = count)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(x = "Bedrock type", y = "Number of occurances", title = "Bedrock geology counts")

write_csv(st_drop_geometry(bedrock_geology_matrix), "data/bedrock_geology_matrix.csv")
