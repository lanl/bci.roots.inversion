
#-----------------------------------------------------
# Title: Preparing isotopic data for inverse model evaluation
# Author : Rutuja Chitra-Tarak
# Original date: December 18, 2019
#-----------------------------------------------------

rm(list=ls())
gc()
pacman::p_load(tidyverse, scales, ggpmisc)
# graphics info
theme_set(theme_bw())
theme_update(text = element_text(size = 14),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)
iso <- read.csv("data-raw/traits/isotopes/Meinzer1999_Table1_Xylem_Sap_deltaD_Fig4.csv", na.strings = c("NA", ""), header = T, row.names = NULL, check.names = F)
tlp <- read.csv("data-raw/traits/HydraulicTraits_Kunert/tlp_sp_mean.csv", na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)

head(iso)
# adding species codes
load(file.path("data-raw/CTFScensuses/bci.spptable.rdata"))
write.csv(bci.spptable, file.path("data-raw/CTFScensuses/bci.spptable.csv"), row.names = FALSE)

iso$genus.sp <- paste(iso$Genus, iso$Species, sep = " ")
bci.spptable$genus.sp <- paste(bci.spptable$Genus, bci.spptable$Species, sep = " ")
matchrows = match(iso$genus.sp, bci.spptable$genus.sp)
iso$genus.sp

iso$sp <- bci.spptable$sp[matchrows]

iso$sp
## some codes are not found, because of change of latin names, so searching in synonymous name-field should return the match
match("Lonchocarpus", bci.spptable$Genus) #"Lonchocarpus latifolius"

grep("Lonchocarpus", bci.spptable$syn)
iso$sp[which(iso$Genus == "Lonchocarpus")] <- bci.spptable$sp[grep("Lonchocarpus", bci.spptable$syn)]
head(iso)

#"Guapira" ; iso$Family[11]
grep("Guapira", bci.spptable$syn)
iso$sp[which(iso$Genus == "Guapira")] <- bci.spptable$sp[grep("Guapira", bci.spptable$syn)]# "Guapira" "campochagres"
iso$sp

## the following species code need to be matched with that in ds.bestfit
# isotopic uptake depth is available for Quararibea asterolepis subsp. stenophylla, but
# modelled uptake depth is for Quararibea asterolepis, so making them synonymous for comparison
iso$sp[which(iso$sp == "quara1")] <- "quaras"
## adding phonology score

iso <- iso %>% mutate(se = SE,
                      source = "Meinzer et al.1999 Fig. 4",
                      location = "BCI",
              PhenoRank = recode(Phenology, `Evergreen` = 1, `Facultatively deciduous` = 2,
                        Brevideciduous = 3,  `Partially deciduous` = 4,  `Deciduous` = 5),
              Phenology = factor(Phenology, levels = c("Evergreen", "Facultatively deciduous",
                                                  "Brevideciduous", "Partially deciduous", "Deciduous")))

iso.sp <- iso$sp #[!is.na(iso$Xylem_sap_deltaD_permil)]
save(iso.sp, file = "results/sp_with_isotopic_record.Rdata")
load("results/sp_with_isotopic_record.Rdata")
write.csv(iso, file.path(paste0("data-raw/traits/isotopes/Meinzer1999_Table1_Xylem_Sap_deltaD_Fig4_sp_code.csv")), row.names = FALSE)

### More data from March1997 from Fig 5B (includes data in Fig 5A)
iso.2.raw <- read.csv("data-raw/traits/isotopes/Meinzer1999_Xylem_Sap_deltaD_March97_DBH_Fig5B.csv", na.strings = c("NA", ""), header = T, row.names = NULL, check.names = F)
iso.change <- read.csv("data-raw/traits/isotopes/Meinzer1999_Xylem_Sap_deltaD_by_delta_SapFlow_&_monthly_LeafFall_Fig7AB.csv", na.strings = c("NA", ""), header = T, row.names = NULL, check.names = F)

iso.2.raw <- iso.2.raw %>% mutate(source = "Meinzer et al.1999 Fig. 5A", location = "BCI")
iso.change <- iso.change %>% mutate(source = "Meinzer et al.1999 Fig. 7AB", location = "BCI")
head(iso.2.raw)
iso.2 <- iso.2.raw %>%
  group_by(sp, source) %>%
  summarise(Xylem_sap_deltaD_permil = mean(Xylem_sap_deltaD_permil, na.rm = TRUE),
            n = n(),
            SD = sd(Xylem_sap_deltaD_permil, na.rm = TRUE),
            SE = SD/sqrt(n),
            DBH = mean(DBH, na.rm = TRUE))

iso.2 <- iso.2 %>% left_join(iso %>% select(sp, Phenology), by = "sp") %>%
  left_join(iso.change, by = "sp") %>%
  arrange(Phenology)
View(iso.2)

#### Plotting soil isotopic variation from Meinzer et al and Jackson et al -----------
iso.soil.1 <- read.csv("data-raw/traits/isotopes/Meinzer1999_Fig2A_soil_deltaD_BCI_data_Mar_Apr_1997.csv", na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)
iso.soil.2 <- read.csv("data-raw/traits/isotopes/Oecologia 1995 Jackson_fig2_soil_deltaD_Gigante_data_dry_season_1992.csv", na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)

depth.m1 <- lm(depth ~ soil.deltaD, data = iso.soil.1)
summ.depth.m1 <- summary(depth.m1)
soil.label <- expression('Soil '*delta~""^2*"H (\u2030)"*'')
## expression(paste(delta^{2}, "H (\u2030)")))
depth.m1.label = paste0("Depth = ", round(depth.m1$coefficients[1], 2), " + ",
                        round(depth.m1$coefficients[2], 2),
                        " * Soil delta 2H\nR-squared = ", round(summ.depth.m1$r.squared, 2),
                        ", p-val = ", round(summ.depth.m1$coefficients[2, 4], 4))

g1 <- ggplot(iso.soil.1, aes(y = depth, x = soil.deltaD)) +
  geom_point(size = 3) +
  xlab(soil.label) + ylab("Depth (cm)") +
  scale_y_reverse()
g1 +   geom_smooth(method = "lm", se = FALSE) +
  geom_errorbarh(aes(xmax = soil.deltaD + se, xmin = soil.deltaD - se), size = 0.5) +
  geom_text(aes(x = -40, y = 0, label = depth.m1.label), color = "blue")
ggsave(file.path(paste0("figures/UDI_confidence/Meinzer_etal_1999_Depth_vs_soil_deltaD.jpeg")), height = 5, width = 6, units='in')

g1 +   geom_smooth(method = "loess", se = FALSE, span = 0.7) +
  geom_errorbarh(aes(xmax = soil.deltaD + se, xmin = soil.deltaD - se), size = 0.5)

depth.m2 <- lm(depth ~ soil.deltaD, data = iso.soil.2)
summ.depth.m2 <- summary(depth.m2)
depth.m2.label = paste0("Depth = ", round(depth.m2$coefficients[1], 2), " + ",
                        round(depth.m2$coefficients[2], 2),
                        " * Soil delta 2H\nR-squared = ", round(summ.depth.m2$r.squared, 2),
                        ", p-val = ", round(summ.depth.m2$coefficients[2, 4], 4))
g1 %+% iso.soil.2 + geom_smooth(method = "glm", se = FALSE) +
  geom_text(aes(x = -35, y = 5, label = depth.m2.label), color = "blue")
ggsave(file.path(paste0("figures/UDI_confidence/Jackson_etal_1995_Depth_vs_soil_deltaD.jpeg")), height = 5, width = 6, units='in')


## adding Xylem sap's soil uptake depth ----
# how to incorporate se in this?
iso <- iso %>% mutate(depth = predict.lm(depth.m1, newdata = data.frame(soil.deltaD = Xylem_sap_deltaD_permil))/100,
                      depth.se = depth - predict.lm(depth.m1,
                                                    newdata = data.frame(soil.deltaD = Xylem_sap_deltaD_permil + SE))/100) # from cm to m
head(iso)


xylem.label <- expression('Xylem Sap '*delta~""^2*"H (\u2030)"*'')
change.xylem.label <- expression('Change in Xylem Sap '*delta~""^2*"H (\u2030)"*'')

iso.2.raw <- iso.2.raw %>% left_join(iso %>% select(sp, Phenology), by = "sp") %>%
  left_join(iso.change, by = "sp") %>% left_join(tlp %>% mutate(sp = as.character(sp)), by = "sp")


ggplot(iso.2.raw, aes(y =  Xylem_sap_deltaD_permil, x = DBH, color = tlp)) +
  geom_text(aes(y =  Xylem_sap_deltaD_permil + 1.5, x = DBH, label = sp),
            size = 2) +
  geom_point(aes(shape = Phenology), size = 3) +
  scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
  scale_y_continuous(limits = c(-70, -18), breaks = seq(-70, -20, by = 10)) +
  xlab("DBH (cm)") + ylab(xylem.label) +
  ggtitle("Species March 1997 Xylem Sap deltaD Vs DBH")
ggsave(file.path(paste0("data-raw/traits/isotopes/Meinzer1999_Xylem_Sap_deltaD_March97_DBH_Fig5B.jpeg")), height = 5, width = 8, units ='in')

ggplot(iso.2.raw, aes(y =  delta_sapflow_percent_per_day, x = delta_xylem_sap_deltaD_permil_per_day, color = tlp)) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5, color = "black") +
  geom_text(aes(y =  delta_sapflow_percent_per_day + 0.05, label = sp),
            size = 4) +
  geom_point(aes(shape = Phenology), size = 3) +
  scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
  scale_x_continuous(limits = c(-0.6, 0.1)) +
  xlab(change.xylem.label) + ylab("Change in daily sap flow (% per day)") +
  ggtitle("Sap flow Vs. Xylem Sap deltaD\nWithin species variation in 1997 dry season")
ggsave(file.path(paste0("data-raw/traits/isotopes/Meinzer1999_Xylem_Sap_deltaD_by_delta_SapFlow_&_monthly_LeafFall_Fig7A.jpeg")), height = 5, width = 7, units ='in')

ggplot(iso.2.raw, aes(y = SE_mean_monthly_leaf_fall, x = delta_xylem_sap_deltaD_permil_per_day, color = tlp)) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5, color = "black") +
  geom_text(aes(y = SE_mean_monthly_leaf_fall + 0.002, label = sp),
            size = 4) +
  geom_point(aes(shape = Phenology), size = 3) +
  scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
  scale_x_continuous(limits = c(-0.6, 0.1)) +
  xlab(change.xylem.label) + ylab("SE of mean monthly leaf fall") +
  ggtitle("Leaf fall Vs. Xylem Sap deltaD\nWithin species variation in 1997 dry season")
ggsave(file.path(paste0("data-raw/traits/isotopes/Meinzer1999_Xylem_Sap_deltaD_by_delta_SapFlow_&_monthly_LeafFall_Fig7B.jpeg")), height = 5, width = 7, units ='in')


iso.3 <- read.csv("data-raw/traits/isotopes/Oecologia 1995 Jackson _Fig3_Fig4.csv", na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)
iso.3 <- iso.3 %>% arrange(Xylem_sap_deltaD_permil) %>%
  left_join(tlp %>% mutate(sp = as.character(sp)), by = "sp") %>%
  mutate(source = "Jackson et al. 1995")

iso.udi <- left_join(iso.3, udi, by = "sp")
iso.sp <- unique(iso.udi$sp[!is.na(iso.udi$udi.med.rsq) & !is.na(iso.udi$Xylem_sap_deltaD_permil)])
iso.udi <- iso.udi[order(iso.udiudi.med.rsq), ]

## Combine all isotopic records:
View(iso)

iso.all.list <- list(iso, iso.2.raw, iso.3, iso.change)
names(iso.all.list) <- c(iso$source[1], iso.2.raw$source[1], iso.3$source[1], iso.change$source[1])

iso.all <- iso %>% select(sp, location, source, Xylem_sap_deltaD_permil, se) %>%
  rbind(iso.3 %>% select(sp, location, source, Xylem_sap_deltaD_permil, se)) %>%
  rbind(iso.2.raw %>% mutate(se = NA) %>% select(sp, location, source, Xylem_sap_deltaD_permil, se)) %>%
  left_join(iso.change %>% select(-source, -location), by = "sp")
View(iso.all)
save(iso.all, file = "results/all_isotopic_record.Rdata")

all.iso.sp <- unique(as.character(iso.all$sp))
all.iso.sp <- all.iso.sp[!is.na(all.iso.sp)]

save(all.iso.sp, file = "results/all_sp_with_isotopic_record.Rdata")

