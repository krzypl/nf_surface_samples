library(vegan)
library(tidyverse)
library(ggvegan)
library(readxl)
library(forecast)

#plan na analizy: po utworzeniu wszystkich zbiorow danych (chemia wody, geologia, topografia, macierz dla roslin ladowych, macierz dla roslin wodnych) zrobic nastepujace wykresy:
# PCA/CA dla samych danych dla chemii wody
# RDA/CCA gdzie zmienne wyjasniane to bedzie chemia wody, a zmienne wyjasniajace to bedzie wszystko poza roslinami wodnymi (zastanowic sie jeszcze jak to zrobic z roslinami ladowymi, w zaleznosci od tego ile ich bedzie, czy wlaczyc je w ogolny  zbior wyjasniajacy, czy, kiedy ich bedzie za duzo, cos innego wymyslic). I tutaj dwa wykresy:
## 1) wykres z zaznaczonymi jeziorami jako punktami i wektorami jako zmiennymi wyjaśnianymi
## 2) wykres z zaznaczonymi jezioramo jako punktami i wektorami jako zmiennymi wyjaśniającymi
#Osobno zrobić RDA z roslinami wodnymi jako zmiennymi wyjasnianymi


water_chemistry <- read_excel("data/water_chemistry_data.xlsx")

water_chemistry[-1] <-  lapply(water_chemistry[-1], as.numeric)

water_chemistry_missing_data <- water_chemistry %>% 
  pivot_longer(cols = "Sodium":"Langelier Index (5°C)", names_to = "variable", values_to = "value") %>% 
  group_by(variable) %>% 
  summarise(nna = sum(is.na(value)))

water_chemistry_variables_2_remove <- water_chemistry_missing_data %>% 
  filter(nna > 4) %>% 
  pull(variable)

water_chemistry$Iron[is.na(water_chemistry$Iron)] <- 0.02 #dodaje wartosc ponizej poziomu detekcji dla brakujacych danych Fe

water_chemistry$Zinc[is.na(water_chemistry$Zinc)] <- 0.0009 #dodaje wartosc ponizej poziomu detekcji dla brakujących danych Zn

water_chemistry_red <- water_chemistry %>% 
  select(!all_of(water_chemistry_variables_2_remove))

field_data <- read_excel("data/nf_field_measurements.xlsx")

plot(water_chemistry_red$pH, field_data$pH) 
reg <- lm(water_chemistry_red$pH ~ field_data$pH)
summary(reg) # nie jest idealnie. ciezko powiedziec czmu jest lepiej bardziej wierzyc. labo czy teren? dla ulatwienia biore na teraz lab

plot(water_chemistry_red$Conductivity, field_data$conductivity) 
reg2 <- lm(water_chemistry_red$Conductivity ~ field_data$conductivity) 
summary(reg2) # tutaj jest bardzo duza zgodnosc. dla ulatwienia biore na teraz lab.



water_chemistry_red_distrib_plot <- water_chemistry_red %>% 
  pivot_longer(cols = "Sodium":"Ion Sum", names_to = "variable", values_to = "value") %>%
  ggplot() +
  geom_histogram(aes(x = value)) +
  facet_wrap(.~variable, scales = "free") #wiekszosc rozkladow mocno skosne

water_chemistry_red2 <- water_chemistry_red %>% 
  dplyr::select(!c("Percent Difference", "Theoretical Conductivity")) %>% 
  rename(anion_sum = "Anion Sum",
         bicarbonate = "Bicarbonate (as CaCO3)",
         DOC = "Carbon - Dissolved Organic",
         TOC = "Carbon - Total Organic",
         carbonate = "Carbonate (as CaCO3)",
         hardness = "Hardness (as CaCO3)",
         Ca = Calcium,
         cation_sum = "Cation Sum",
         Cl = "Chloride",
         Fe = "Iron",
         ion_sum = "Ion Sum",
         Mg = "Magnesium",
         Mn = "Manganese",
         K = "Potassium",
         Na = "Sodium",
         Zn = "Zinc"
         
  )

lakeID <- water_chemistry_red$LakeID

field_data <- read_excel("data/nf_field_measurements.xlsx")

field_data_2add <- field_data %>% 
  select(lakeID, DO_percent) #nie uwzgledniam temperatury wody, bo jesli to robie, to mam oczywisty wplyw temperatury powietrza i sztucznie moj model wyjasnia duzo zmiennosci

sediment_data <- read_excel("data/lake_sediment_data.xlsx")



#zrezygnowalem z BoxCoxa, zeby indywidualnie dobrac transformacje
water_chemistry_boxcox_prep <- water_chemistry_red2 %>%
  left_join(field_data_2add, by = c("LakeID" = "lakeID")) %>% 
  left_join(sediment_data, by = c("LakeID" = "lakeID"))

write_csv(water_chemistry_boxcox_prep, "data/water_chemistry_for_distributions.csv")


water_chemistry_boxcox <- water_chemistry_boxcox_prep %>%
#  mutate(across(-1, ~ BoxCox(.x, lambda = BoxCox.lambda(.x))))
  mutate(Na = log(Na),
         K = log(K),
         Ca = log(Ca),
         Mg = log(Mg),
         Fe = log(Fe),
         Mn = log(Mn),
         Zn = log(Zn),
         Cl = log(Cl),
         Turbidity = log(Turbidity),
         Conductivity = log(Conductivity),
         LOI950 = log(LOI950),
         sand_grain = log(sand_grain+1))

water_chemistry_red_distrib_plot_transformed <- water_chemistry_boxcox %>% 
  pivot_longer(cols = Na:sand_grain, names_to = "variable", values_to = "value") %>%
  ggplot() +
  geom_histogram(aes(x = value)) +
  facet_wrap(.~variable, scales = "free") #duzo lepiej wygladaja te rozklady teraz dla zmiennych, ktore wykorzystuje w analizach. nie jest idealnie, ale do PCA to sie nadaje. pozosta


water_chemistry_zscore_prep <- water_chemistry_boxcox %>%
  select(-c(hardness, `Hydroxide (as CaCO3)`, bicarbonate, anion_sum, cation_sum, ion_sum, carbonate, Cl, TOC)) #usuwam zmienne kolinearne; CLma bardzo silna korelacje z Na

write_rds(water_chemistry_zscore_prep, "data/water_chemistry_for_ordinations.rds")

water_chemistry_zscore <- water_chemistry_zscore_prep %>% 
  dplyr::select(!LakeID) %>% 
  mutate(across(everything(), ~ as.numeric(scale(.))))

correlation_matrix_water_chemistry <- cor(water_chemistry_zscore)
  
water_chemistry_pca <- rda(water_chemistry_zscore) 

screeplot(water_chemistry_pca, bstick = TRUE)

water_chemistry_fort <- fortify(water_chemistry_pca, axes = c(1,2), scaling = "species") #scaling jest na species, wiec wiernie beda odzwierciedlone przede wszystkim relacje pomiedzy zeminnymi

water_chemistry_sites <- 
  water_chemistry_fort[water_chemistry_fort$score %in% "sites",] %>% 
  mutate(lakeID = lakeID)

water_chemistry_sp <- water_chemistry_fort[water_chemistry_fort$score %in% "species",]

water_chemistry_ve_prep <- water_chemistry_pca$CA$eig / water_chemistry_pca$tot.chi * 100
(water_chemistry_PC1_ve <- round(((water_chemistry_ve_prep / sum(water_chemistry_ve_prep))[c(1)]) * 100, digits = 1))
(water_chemistry_PC2_ve <- round(((water_chemistry_ve_prep / sum(water_chemistry_ve_prep))[c(2)]) * 100, digits = 1))

#PCA biplot ----

water_chemistry_pca_plot <- ggplot() +
  labs(y = paste("PC2 (", water_chemistry_PC2_ve, "%)", sep = ""), x = paste("PC1 (", water_chemistry_PC1_ve, "%)", sep = ""), title = "PCA biplot for water chemistry data") +
  geom_segment(data = water_chemistry_sp,
               color = "black", linewidth = 0.7,
               aes(x = 0, y = 0, xend = pc1, yend = pc2),
               arrow = 
                 grid::arrow(length = grid::unit(0.25, "cm"))) +
  geom_point(data = water_chemistry_sites, aes(x = pc1, y = pc2)) +
  ggrepel::geom_text_repel(data = water_chemistry_sites, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = pc1, y = pc2, 
                               label = lakeID)) +
  ggrepel::geom_text_repel(data = water_chemistry_sp, color = "red",
                           size = 4, segment.alpha = 0,
                           aes(x = pc1, y = pc2, 
                               label = label)) +
  geom_vline(xintercept = 0, color = 'black', linewidth = 0.6,linetype=2) + 
  geom_hline(yintercept = 0, color = 'black', linewidth = 0.6,linetype=2) +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(color = "grey80", linewidth = 0.2)) +
  theme(legend.position = "right") +
  coord_fixed()

ggsave("figures/water_chemistry_pca.svg",
       plot = water_chemistry_pca_plot, width = 8, height = 6)

#RDA ----

elevation_gps <- read_excel("data/elevation_gps.xlsx") %>% 
  rename(dtto = distance_to_the_ocean)
bedrock_geology <- read_csv("data/bedrock_geology_matrix.csv") %>% 
  select(!geology_type)
surface_geology <- read_csv("data/surface_geology_matrix.csv") %>% 
  select(-c(geology_type, covered_bedrock))
ecoregions <- read_csv("data/ecoregions_matrix.csv")
spatial_factors <- read_csv("data/spatial_factors.csv")
  
relief <- read_csv("data/relief.csv")
relief$relief <- as.factor(relief$relief)
temperatures_for_coring_sites <- read_csv("data/temperature_for_coring_sites.csv")
winter_temperatures_for_coring_sites <- read_csv("data/winter_temperature_for_coring_sites.csv")
precipitation_for_coring_sites <- read_csv("data/precipitation_for_coring_sites.csv")
winter_precipitation_for_coring_sites <- read_csv("data/winter_precipitation_for_coring_sites.csv")

x <- temperatures_for_coring_sites %>% select(summer_mean_temperature:july_temperature)
cor(x) #nie jest niespodzianka, ze temperatury sa ze soba bardzo silnie skorelowane, dlatego do modelu wykorzystuje tylko srednia temperature lata

env_data_prep <- elevation_gps %>% 
  left_join(temperatures_for_coring_sites, by = c("lakeID" = "lakeID")) %>% 
  left_join(winter_temperatures_for_coring_sites, by = c("lakeID" = "lakeID")) %>% 
  left_join(precipitation_for_coring_sites, by = c("lakeID" = "lakeID")) %>% 
  left_join(winter_precipitation_for_coring_sites, by = c("lakeID" = "lakeID")) %>% 
  select(lakeID, elevation, dtto, area_ha, summer_mean_temperature, winter_mean_temperature, summer_sum_precipitation, winter_sum_precipitation) %>% 
  mutate(depth = field_data$maximum_measured_depth_m)

env_data_raw_plot <- env_data_prep %>% 
  pivot_longer(cols = elevation:depth,
               names_to = "variable", values_to = "value") %>% 
  ggplot() +
  geom_histogram(aes(x = value)) +
facet_wrap(.~variable, scales = "free")

ggsave("figures/env_data_raw_plot.svg",
       plot = env_data_raw_plot, width = 8, height = 6)

env_data <- env_data_prep %>%
#  mutate(across(-1, ~ BoxCox(.x, lambda = BoxCox.lambda(.x)))) %>% 
  mutate(area_ha = log(area_ha),
         depth = log(depth),
         dtto = log(dtto),
         elevation = sqrt(elevation),
         temp_amplitude =temperatures_for_coring_sites$july_temperature -
           winter_temperatures_for_coring_sites$january_temperature,
         across(-1, ~ as.numeric(scale(.)))) %>% 
  left_join(bedrock_geology, by = c("lakeID" = "lakeID")) %>% 
  left_join(surface_geology, by = c("lakeID" = "lakeID")) %>% 
  left_join(ecoregions, by = c("lakeID" = "lakeID")) %>% 
  left_join(spatial_factors, by = c("lakeID" = "lakeID")) %>% 
  left_join(relief, by = c("lakeID" = "lakeID")) %>%
  select(-lakeID)

env_data_transformed_plot <- env_data %>% 
  select(area_ha, depth, dtto, elevation) %>% 
  mutate(lakeID = lakeID) %>%
  pivot_longer(!lakeID, names_to = "variable", values_to = "value") %>% 
  ggplot() +
  geom_histogram(aes(x = value)) +
  facet_wrap(.~variable, scales = "free")

ggsave("figures/env_data_transformed_plot.svg",
       plot = env_data_transformed_plot, width = 8, height = 6)


rda_all <- rda(water_chemistry_zscore ~ ., data = env_data)
summary(rda_all)
(R2a_all <- RsquareAdj(rda_all)$adj.r.squared)
mod0 <- rda(water_chemistry_zscore ~ 1, data = env_data)

set.seed(12)
step_forward <- ordistep(mod0,
                         scope = formula(rda_all),
                         direction = "forward",
                         permutations = how(nperm = 999) ) 
(RsquareAdj(step_forward)$adj.r.squared)  


env_variables_forward_selected <- env_data %>% 
  select(c(carbonate, ycoordinate, depth, xcoordinate, MEM_broad, MEM_fine,
           `Northern Peninsula Forest`, `Maritime Barrens`, `Eastern Hyper-Oceanic Barrens`))

rda_fs <- rda(water_chemistry_zscore ~ ., data = env_variables_forward_selected)
(R2a_all_fs <- RsquareAdj(rda_fs)$adj.r.squared)

set.seed(12)
all_anova_fs <- anova(rda_fs, permutations = 999)

set.seed(12)
all_anova_cca <- anova.cca(rda_fs, by = "axis")

set.seed(12)
rda_carbonate <- rda(water_chemistry_zscore ~ carbonate + Condition(ycoordinate + depth + xcoordinate + MEM_broad + MEM_fine + `Northern Peninsula Forest` + `Maritime Barrens` + `Eastern Hyper-Oceanic Barrens`), data = env_variables_forward_selected)

(R2a_carbonate_fs <- RsquareAdj(rda_carbonate)$adj.r.squared)

set.seed(12)
carbonate_anova_fs <- anova(rda_carbonate, permutations = 999)

set.seed(12)
rda_xy_coordinate <- rda(water_chemistry_zscore ~ ycoordinate + xcoordinate + Condition(carbonate + depth + MEM_broad + MEM_fine + `Northern Peninsula Forest` + `Maritime Barrens` + `Eastern Hyper-Oceanic Barrens`), data = env_variables_forward_selected)

(R2a_xy_coordinate_fs <- RsquareAdj(rda_xy_coordinate)$adj.r.squared)

set.seed(12)
xy_coordinate_anova_fs <- anova(rda_xy_coordinate, permutations = 999)

set.seed(12)
rda_depth <- rda(water_chemistry_zscore ~ depth + Condition(ycoordinate + xcoordinate + carbonate + MEM_broad + MEM_fine + `Northern Peninsula Forest` + `Maritime Barrens` + `Eastern Hyper-Oceanic Barrens`), data = env_variables_forward_selected)

(R2a_depth_fs <- RsquareAdj(rda_depth)$adj.r.squared)

set.seed(12)
area_anova_fs <- anova(rda_depth, permutations = 999)

set.seed(12)
rda_MEM_broad <- rda(water_chemistry_zscore ~ MEM_broad + Condition(ycoordinate + xcoordinate + carbonate + depth + MEM_fine + `Northern Peninsula Forest` + `Maritime Barrens` + `Eastern Hyper-Oceanic Barrens`), data = env_variables_forward_selected)

(R2a_MEM_broad_fs <- RsquareAdj(rda_MEM_broad)$adj.r.squared)

set.seed(12)
MEM_broad_anova_fs <- anova(rda_MEM_broad, permutations = 999)

set.seed(12)
rda_MEM_fine <- rda(water_chemistry_zscore ~ MEM_fine + Condition(ycoordinate + xcoordinate + carbonate + depth + MEM_broad + `Northern Peninsula Forest` + `Maritime Barrens` + `Eastern Hyper-Oceanic Barrens`), data = env_variables_forward_selected)

(R2a_MEM_fine_fs <- RsquareAdj(rda_MEM_fine)$adj.r.squared)

set.seed(12)
MEM_fine_anova_fs <- anova(rda_MEM_fine, permutations = 999)

set.seed(12)
rda_np_forest <- rda(water_chemistry_zscore ~ `Northern Peninsula Forest` + Condition(ycoordinate + xcoordinate + carbonate + depth + MEM_broad + MEM_fine + `Maritime Barrens` + `Eastern Hyper-Oceanic Barrens`), data = env_variables_forward_selected)

(R2a_np_forest_fs <- RsquareAdj(rda_np_forest)$adj.r.squared)

set.seed(12)
np_forest_anova_fs <- anova(rda_np_forest, permutations = 999)

set.seed(12)
rda_m_barrens <- rda(water_chemistry_zscore ~ `Maritime Barrens` + Condition(ycoordinate + xcoordinate + carbonate + depth + MEM_broad + MEM_fine + `Northern Peninsula Forest` + `Eastern Hyper-Oceanic Barrens`), data = env_variables_forward_selected)

(R2a_m_barrens_fs <- RsquareAdj(rda_m_barrens)$adj.r.squared)

set.seed(12)
m_barrens_anova_fs <- anova(rda_m_barrens, permutations = 999)

set.seed(12)
rda_ho_barrens <- rda(water_chemistry_zscore ~ `Eastern Hyper-Oceanic Barrens` + Condition(ycoordinate + xcoordinate + carbonate + depth + MEM_broad + MEM_fine + `Northern Peninsula Forest` + `Maritime Barrens`), data = env_variables_forward_selected)

(R2a_ho_barrens_fs <- RsquareAdj(rda_ho_barrens)$adj.r.squared)

set.seed(12)
ho_barrens_anova_fs <- anova(rda_ho_barrens, permutations = 999)

#partial RDA column plot -----
rda_tprep <- tibble(
  "Carbonate" = c(R2a_carbonate_fs, carbonate_anova_fs$`Pr(>F)`[[1]]),
  "XY_trend" = c(R2a_xy_coordinate_fs, xy_coordinate_anova_fs$`Pr(>F)`[[1]]),
  "MEM_broad " = c(R2a_MEM_broad_fs, MEM_broad_anova_fs$`Pr(>F)`[[1]]),
  "MEM_fine" = c(R2a_MEM_fine_fs, MEM_fine_anova_fs$`Pr(>F)`[[1]]),
  "NP_forest" = c(R2a_np_forest_fs, np_forest_anova_fs$`Pr(>F)`[[1]]),
  "M_Barrens" = c(R2a_m_barrens_fs, m_barrens_anova_fs$`Pr(>F)`[[1]]),
  "HO_Barrens" = c(R2a_ho_barrens_fs, ho_barrens_anova_fs$`Pr(>F)`[[1]]),
  "All" = c(R2a_all_fs, all_anova_fs$`Pr(>F)`[[1]]),
  "Unexplained" = c(1 - R2a_all_fs, NA),
  "Metric" = c("prop_of_var_explained", "pval")
) %>% 
  pivot_longer(cols = Carbonate:Unexplained, names_to = "variable", values_to = "value")

prda_plot <- ggplot(filter(rda_tprep, Metric == "prop_of_var_explained")) +
  geom_col(aes(x = variable, y = value))

ggsave("figures/prda_plot.svg",
       plot = prda_plot, width = 8, height = 6)

#variation partitioning ----
spatial_patterns <- env_variables_forward_selected %>% 
  select(ycoordinate, xcoordinate, MEM_broad, MEM_fine)

local_environment <- env_variables_forward_selected %>% 
  select(carbonate, `Northern Peninsula Forest`, `Maritime Barrens`, `Eastern Hyper-Oceanic Barrens`)

morphometric_paramters <- env_variables_forward_selected %>% 
  select(depth)

# climate <- env_data %>% 
#   select(summer_mean_temperature, winter_mean_temperature)

varpart_model <- varpart(
  water_chemistry_zscore,
  spatial_patterns,
  local_environment,
  morphometric_paramters
#  climate
)


svg("figures/venn_plot_rda.svg", width = 8, height = 6)
plot(
  varpart_model,
  bg = c("skyblue","orange","forestgreen"),
  Xnames = c("Spatial patterns","Local environmant","Morphometric parameters"),
  id.size = 1.2
)
dev.off()

#RDA biplot ----
rda_fort <- fortify(rda_fs, axes(1,2), scaling = "species") # First and second axes

rda_water_chemistry_sp <- rda_fort %>% 
  filter(score == "species")

rda_expl_var <- rda_fort %>% 
  filter(score == "biplot")

rda_sites <- rda_fort %>% 
  filter(score == "sites")


rda_ve_prep_sites <- rda_fs$CCA$eig / rda_fs$tot.chi*100
(rda1_ve <- round(((rda_ve_prep_sites / sum(rda_ve_prep_sites))[c(1)] * R2a_all_fs) * 100, digits = 1))
(rda2_ve <- round(((rda_ve_prep_sites / sum(rda_ve_prep_sites))[c(2)] * R2a_all_fs) * 100, digits = 1))


rda_plot_vectors_only <- ggplot() +
  labs(y = paste("RDA2 (", rda2_ve, "%)", sep = ""), 
       x = paste("RDA1 (", rda1_ve, "%)", sep = "")) +
  
  geom_segment(data = rda_water_chemistry_sp ,
               color = "red", linewidth = 0.4,
               aes(x = 0, y = 0, xend = rda1, yend = rda2),
               arrow = grid::arrow(length = grid::unit(0.20, "cm"))) +
  
  ggrepel::geom_text_repel(data = rda_water_chemistry_sp, color = "red",
                           size = 4, segment.alpha = 0,
                           aes(x = rda1, y = rda2, label = label, fontface = "italic"),
                           max.overlaps = 30) +
  
  geom_segment(data = rda_expl_var,
               color = "darkblue", linewidth = 1,
               aes(x = 0, y = 0, xend = rda1, yend = rda2),
               arrow = grid::arrow(length = grid::unit(0.20, "cm"))) +
  
  ggrepel::geom_text_repel(data = rda_expl_var, color = "darkblue",
                           size = 5, segment.alpha = 0,
                           aes(x = rda1, y = rda2, label = label)) +
  
  geom_vline(xintercept = 0, color = 'grey65', linewidth = 0.6, linetype = 2) + 
  geom_hline(yintercept = 0, color = 'grey65', linewidth = 0.6, linetype = 2) +
  
  theme_bw() +
  scale_size(range = 2) +
  coord_equal()


ggsave("figures/water_chemistry_rda_vectors_only.svg",
       plot = rda_plot_vectors_only, width = 8, height = 6)


rda_plot_vectors_and_sites <- ggplot() +
  labs(y = paste("RDA2 (", rda2_ve, "%)", sep = ""), 
       x = paste("RDA1 (", rda1_ve, "%)", sep = "")) +
  
  geom_point(data = rda_sites,
             color = "black",
             aes(x = rda1, y = rda2)) +
  
  ggrepel::geom_text_repel(data = rda_sites, color = "black",
                           size = 4, segment.alpha = 0,
                           aes(x = rda1, y = rda2, label = lakeID),
                           max.overlaps = 30) +
  
  geom_segment(data = rda_expl_var,
               color = "darkblue", linewidth = 1,
               aes(x = 0, y = 0, xend = rda1, yend = rda2),
               arrow = grid::arrow(length = grid::unit(0.20, "cm"))) +
  
  ggrepel::geom_text_repel(data = rda_expl_var, color = "darkblue",
                           size = 5, segment.alpha = 0,
                           aes(x = rda1, y = rda2, label = label)) +
  
  geom_vline(xintercept = 0, color = 'grey65', linewidth = 0.6, linetype = 2) + 
  geom_hline(yintercept = 0, color = 'grey65', linewidth = 0.6, linetype = 2) +
  
  theme_bw() +
  scale_size(range = 2) +
  coord_equal()


ggsave("figures/water_chemistry_rda_vectors_and_sites.svg",
       plot = rda_plot_vectors_and_sites, width = 8, height = 6)
