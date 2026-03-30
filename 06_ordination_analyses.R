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
  select(lakeID, DO_percent, water_temperature)

sediment_data <- read_excel("data/lake_sediment_data.xlsx")


water_chemistry_boxcox <- water_chemistry_red2 %>%
  left_join(field_data_2add, by = c("LakeID" = "lakeID")) %>% 
  left_join(sediment_data, by = c("LakeID" = "lakeID")) %>% 
  mutate(across(-1, ~ BoxCox(.x, lambda = BoxCox.lambda(.x))))

lambdas <- water_chemistry_red2 %>%
  summarise(across(-LakeID, BoxCox.lambda)) # to sa wartosci lambda, ktore zostaly wkorzystane do transformacji boxa-coxa

water_chemistry_red_distrib_plot_transformed <- water_chemistry_boxcox %>% 
  pivot_longer(cols = Na:water_temperature, names_to = "variable", values_to = "value") %>%
  ggplot() +
  geom_histogram(aes(x = value)) +
  facet_wrap(.~variable, scales = "free") #duzo lepiej wygladaja te rozklady teraz. nie jest idealnie, ale do PCA to sie nadaje


water_chemistry_zscore <- water_chemistry_boxcox %>% 
  dplyr::select(!LakeID) %>% 
  mutate(across(everything(), ~ as.numeric(scale(.)))) %>% 
  select(-c(hardness, `Hydroxide (as CaCO3)`, bicarbonate, anion_sum, cation_sum, ion_sum, carbonate)) #usuwam zmienne kolinearne
  
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


