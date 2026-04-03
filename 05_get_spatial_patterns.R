library(adespatial)

cores_coordinates <- read_delim("data/cores_coordinates.csv",
                                delim = ";",
                                escape_double = FALSE, trim_ws = TRUE) %>% 
  select(-lakeID)

lakeID <- read_delim("data/cores_coordinates.csv",
                                          delim = ";",
                                          escape_double = FALSE, trim_ws = TRUE) %>% 
  select(lakeID) %>% 
  pull(lakeID)

cores_coordinates_sf <- vect(cores_coordinates, geom = c("longitude_core", "latitude_core"), crs = "EPSG:4326")

cores_coordinates_utm_prep <- project(cores_coordinates_sf, "EPSG:26720")

cores_coordinates_utm <- st_as_sf(cores_coordinates_utm_prep)

cores_coordinates_metric  <- cores_coordinates_utm %>%
  st_coordinates() %>%     
  as_tibble() %>%          
  rename(x = X, y = Y)     

cores_dbmem_prep <- dbmem(cores_coordinates_metric)

cores_dbmem <- as_tibble(cores_dbmem_prep)

(thr <- give.thresh(dist(cores_coordinates_metric)))

water_chemistry <- read_rds("data/water_chemistry_for_ordinations.rds") %>% 
  select(-LakeID) %>% 
  mutate(across(everything(), ~ as.numeric(scale(.))))

anova(rda(water_chemistry, cores_coordinates_metric))
water_chemistry_det <- as_tibble(resid(lm(as.matrix(water_chemistry) ~ .,
                                data = cores_coordinates_metric)))

cores_dbmem_rda <- rda(water_chemistry_det ~., cores_dbmem)
anova(cores_dbmem_rda)

(cores_R2a <- RsquareAdj(cores_dbmem_rda)$adj.r.squared)

(cores_dbmem_fwd <- forward.sel(water_chemistry_det, cores_dbmem,
                                adjR2thresh = cores_R2a)) 

(nb_sig_dbmem <- nrow(cores_dbmem_fwd)) # Number of signif. dbMEM  # Identity of the significant dbMEM in increasing order 
(dbmem_sign <- sort(cores_dbmem_fwd[ ,2]))  # Write the significant dbMEM to a new object dbmem.red <- mite.dbmem[ ,c(dbmem.sign)]

spatial_factors <- tibble(MEM_broad = cores_dbmem$MEM5, #tutaj i nizej sprawdzic czy to dobre memy sa, tj. te, ktore wyszly jako istotne statystycznie
                 MEM_fine = cores_dbmem$MEM9,
                 xcoordinate = cores_coordinates_metric$x,
                 ycoordinate = cores_coordinates_metric$y,
                 lakeID = lakeID)

write_csv(spatial_factors, "data/spatial_factors.csv")

xyz <- tibble(x = cores_coordinates_metric$x, y = cores_coordinates_metric$y,
              mem5 = cores_dbmem$MEM5, mem9 = cores_dbmem$MEM9)


moran.randtest(mem[,9], listw)