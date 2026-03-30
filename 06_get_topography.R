library(tidyverse)
library(sf)
library(terra)
library(elevatr)
library(readxl)

cores_coordinates <- read_delim("data/cores_coordinates.csv",
                                delim = ";",
                                escape_double = FALSE, trim_ws = TRUE)

#relief ----

#relief is from field observations
#four classes are distinguished for relief: very low (flat), low (rolling), low to moderate (rolling to hilly), moderete (hilly)


relief <- tribble(
  ~lakeID, ~relief,
  "NF_001", "moderate",
  "NF_002", "moderate",
  "NF_003", "low to moderate",
  "NF_004", "moderate",
  "NF_005", "low",
  "NF_006", "moderate",
  "NF_007", "very low",
  "NF_008", "very low",
  "NF_009", "moderate",
  "NF_010", "low to moderate",
  "NF_011", "moderate",
  "NF_012", "very low",
  "NF_013", "very low",
  "NF_014", "low",
  "NF_015", "low",
  "NF_016", "very low",
  "NF_017", "very low",
  "NF_018", "moderate",
  "NF_019", "moderate",
  "NF_020", "moderate",
  "NF_021", "low to moderate",
  "NF_022", "low",
  "NF_023", "low",
  "NF_024", "moderate",
  "NF_025", "very low",
  "NF_026", "very low",
  "NF_027", "low",
  "NF_028", "moderate",
  "NF_029", "low",
  "NF_030", "very low",
  "NF_031", "very low",
  "NF_032", "very low",
  "NF_033", "very low",
  "NF_034", "very low",
  "NF_035", "moderate",
  "NF_036", "low to moderate",
  "NF_037", "low to moderate",
  "NF_038", "moderate",
  "NF_039", "moderate",
  "NF_040", "moderate",
  "NF_041", "very low",
  "NF_042", "very low",
  "NF_043", "very low",
  "NF_044", "very low",
  "NF_045", "very low",
  "NF_046", "low to moderate",
  "NF_047", "low",
  "NF_048", "low",
  "NF_049", "very low",
  "NF_050", "low to moderate",
  "NF_051", "low",
  "NF_052", "very low",
  "NF_053", "low",
  "NF_054", "very low",
  "NF_055", "low",
  "NF_056", "low",
  "NF_057", "low",
  "NF_058", "low to moderate",
  "NF_059", "very low",
  "NF_060", "very low",
  "NF_061", "low to moderate",
  "NF_062", "moderate",
  "NF_063", "moderate",
  "NF_064", "moderate",
  "NF_065", "low to moderate",
  "NF_066", "moderate",
  "NF_067", "moderate",
  "NF_068", "low",
  "NF_069", "moderate",
  "NF_070", "moderate"
  )

ggplot(relief, aes(x = relief)) +
     geom_bar()

write_csv(relief, "data/relief.csv")

#elevetion ----

cores_coordinates_sf <- st_as_sf(cores_coordinates, coords = c("longitude_core", "latitude_core"), crs = 4326) 

elev <- get_elev_point(cores_coordinates_sf, src = "aws")

ggplot(elev, aes(x = elevation)) +
  geom_histogram()

elevetion_gps <- read_excel("data/elevetion_gps.xlsx") #tutaj sa tez dane dla powierzchni jezior

ggplot(elevetion_gps, aes(x = elevetion)) +
  geom_histogram()

xy <- tibble(gps = elevetion_gps$elevetion,
             srtm = elev$elevation
               )

ggplot(xy, aes(x = gps, y = srtm)) +
  geom_point() +
  geom_smooth(method = "lm")

plot(elev$elevation, elevetion_gps$elevetion)


#jest troche roznica pomiedzy tym, co pokazuja dane z gpsa i tym co pokazuja dane z srtm. ale to nie sa jakies diametralne rozjady, ogolny obraz sie zgadza. Dla SRTM sa wartosci ujemne i dlatego mysle, ze powininem wziac z gps.



