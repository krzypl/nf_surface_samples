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
water_chemistry_boxcox <- water_chemistry_red2 %>%
  left_join(field_data_2add, by = c("LakeID" = "lakeID")) %>% 
  left_join(sediment_data, by = c("LakeID" = "lakeID")) %>% 
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


water_chemistry_zscore <- water_chemistry_boxcox %>% 
  dplyr::select(!LakeID) %>% 
  mutate(across(everything(), ~ as.numeric(scale(.)))) %>% 
  select(-c(hardness, `Hydroxide (as CaCO3)`, bicarbonate, anion_sum, cation_sum, ion_sum, carbonate, Cl, TOC)) #usuwam zmienne kolinearne; CLma bardzo silna korelacje z Na

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
  theme(legend.position = "right")

#RDA ----

elevation_gps <- read_excel("data/elevation_gps.xlsx")
bedrock_geology <- read_csv("data/bedrock_geology_matrix.csv") %>% 
  select(!geology_type)
surface_geology <- read_csv("data/surface_geology_matrix.csv") %>% 
  select(-c(geology_type, covered_bedrock))
relief <- read_csv("data/relief.csv")
relief$relief <- as.factor(relief$relief)
temperatures_for_coring_sites <- read_csv("data/temperature_for_coring_sites.csv")

x <- temperatures_for_coring_sites %>% select(summer_mean_temperature:july_temperature)
cor(x) #nie jest niespodzianka, ze temperatury sa ze soba bardzo silnie skorelowane, dlatego do modelu wykorzystuje tylko srednia temperature lata

env_data_prep <- elevation_gps %>% 
  left_join(temperatures_for_coring_sites, by = c("lakeID" = "lakeID")) %>% 
  select(lakeID, elevation, area_ha, summer_mean_temperature) %>% 
  mutate(depth = field_data$maximum_measured_depth_m)

env_data_raw_plot <- env_data_prep %>% 
  pivot_longer(cols = elevation:depth,
               names_to = "variable", values_to = "value") %>% 
  ggplot() +
  geom_histogram(aes(x = value)) +
facet_wrap(.~variable, scales = "free")

env_data <- env_data_prep %>%
#  mutate(across(-1, ~ BoxCox(.x, lambda = BoxCox.lambda(.x)))) %>% 
  mutate(area_ha = log(area_ha),
         depth = log(depth),
         elevation = sqrt(elevation),
         across(-1, ~ as.numeric(scale(.)))) %>% 
  left_join(bedrock_geology, by = c("lakeID" = "lakeID")) %>% 
  left_join(surface_geology, by = c("lakeID" = "lakeID")) %>% 
  left_join(relief, by = c("lakeID" = "lakeID")) %>%
  select(-lakeID)

rda_all <- rda(water_chemistry_zscore ~ ., data = env_data)
summary(rda_all)
(R2a_all <- RsquareAdj(rda_all)$adj.r.squared)
mod0 <- rda(water_chemistry_zscore ~ 1, data = env_data)

set.seed(12)
step_forward <- ordistep(mod0,
                         scope = formula(rda_all),
                         direction = "forward",
                         permutations = how(nperm = 999) ) #forward selection retains carbonate, june_temperature, depth
(RsquareAdj(step_forward)$adj.r.squared)  


env_variables_forward_selected <- env_data %>% 
  select(c(carbonate, june_temperature, depth))

rda_fs <- rda(water_chemistry_zscore ~ ., data = env_variables_forward_selected)
(R2a_all_fs <- RsquareAdj(rda_fs)$adj.r.squared)

set.seed(12)
all_anova_fs <- anova(rda_fs, permutations = 999)

set.seed(12)
all_anova_cca <- anova.cca(rda_fs, by = "axis") #two axis significant

set.seed(12)
rda_carbonate <- rda(water_chemistry_zscore ~ carbonate + Condition(june_temperature, depth), data = env_variables_forward_selected)

(R2a_carbonate_fs <- RsquareAdj(rda_carbonate)$adj.r.squared)

set.seed(12)
rda_jun_temp <- rda(water_chemistry_zscore ~ june_temperature + Condition(carbonate, depth), data = env_variables_forward_selected)

(R2a_jun_temp_fs <- RsquareAdj(rda_jun_temp)$adj.r.squared)

set.seed(12)
rda_depth <- rda(water_chemistry_zscore ~ depth + Condition(june_temperature, carbonate), data = env_variables_forward_selected)

(R2a_depth_fs <- RsquareAdj(rda_depth)$adj.r.squared)

