
#-----------------------------------------------------
# Title: Preparing isotopic data for inverse model evaluation
# Author : Rutuja Chitra-Tarak
# Original date: December 18, 2019
#-----------------------------------------------------

#************************
# ------ Libraries -----
#************************

rm(list=ls())
gc()
#*******************************************
####   Load Libraries, Prep for graphics, folders  ####
#*******************************************
#### Written with R version 4 ###
#*******************************************
if (!require("groundhog")) install.packages("groundhog"); library(groundhog)
groundhog.folder <- paste0("groundhog.library")
if(!dir.exists(file.path(groundhog.folder))) {dir.create(file.path(groundhog.folder))}
set.groundhog.folder(groundhog.folder)
groundhog.day = "2021-01-01"
pkgs=c('tidyverse', 'scales', 'ggpmisc')
groundhog.library(pkgs, groundhog.day)


# graphics info
theme_set(theme_bw())
theme_update(text = element_text(size = 14),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)

#************************
## ----- Lood Data -----
#************************
iso <- read.csv("data-raw/traits/isotopes/Meinzer1999_Table1_Xylem_Sap_deltaD_Fig4.csv", na.strings = c("NA", ""), header = T, row.names = NULL, check.names = F)
tlp <- read.csv("data-raw/traits/HydraulicTraits_Kunert/tlp_sp_mean.csv", na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)
deci <- read_excel(file.path("data-raw/traits/nomenclature_R_20190524_Rready_Osvaldo Calderon & JoeWright_expert_opinion.xlsx"))
deci.level_key <- c("Evg" = "1", "DF" = "2", "DB" = "3", "DO" = "4") #c(a = "apple", b = "banana", c = "carrot")

deci <- deci %>% mutate(sp = tolower(sp6)) %>%
  dplyr::select(sp4, sp, deciduous) %>%
  subset(deciduous %in% c("E", "DF", "DO", "DB")) %>%
  mutate(deciduousness = recode_factor(as.factor(deciduous), `E` = "Evergreen", `DO` = "Obligate Deciduous", `DF` = "Facultative Deciduous",
                                       `DB` = "Brevideciduous"), ordered = TRUE) %>%
  transform(deciduousness = factor(deciduousness,
                                   levels = c("Evergreen", "Brevideciduous",
                                              "Facultative Deciduous", "Obligate Deciduous"), ordered = TRUE)) %>%
  dplyr::select(sp4, sp, deciduous, deciduousness)

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
                                                  "Brevideciduous", "Partially deciduous", "Deciduous"), ordered = TRUE)) %>%
  left_join(deci, dplyr::select(-sp4), by = "sp")

iso.sp <- iso$sp #[!is.na(iso$Xylem_sap_deltaD_permil)]
save(iso.sp, file = "results/sp_with_isotopic_record.Rdata")
load("results/sp_with_isotopic_record.Rdata")
write.csv(iso, file.path(paste0("data-raw/traits/isotopes/Meinzer1999_Table1_Xylem_Sap_deltaD_Fig4_sp_code.csv")), row.names = FALSE)

iso <- iso %>%
  unite("deci_sp", deciduous, sp, remove = FALSE) %>%
  mutate(deci_sp.plot = factor(deci_sp, levels=unique(deci_sp[order(Xylem_sap_deltaD_permil)]), ordered=TRUE))

xylem.label <- expression(delta^2*H[xylem]~"( \u2030)")
g1 <- ggplot(iso %>% subset(!is.na(Xylem_sap_deltaD_permil)),
             aes(x = deci_sp.plot, y = Xylem_sap_deltaD_permil)) +
  geom_col(aes(fill = deciduousness)) +
  guides(fill = guide_legend(title = "")) +
  theme(legend.position = "top") +
  ylab(xylem.label) + xlab("Species") +
  geom_errorbar(aes(ymax = Xylem_sap_deltaD_permil + se, ymin = Xylem_sap_deltaD_permil - se), width = 0.2, size = 0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(breaks = seq(from = -60, to = 0, by = 5)) +
  theme(panel.grid.major = element_line(colour = "gray", size = 0.1))
ggsave(file.path(paste0("data-raw/traits/isotopes/Meinzer1999_Table1_Xylem_Sap_deltaD_Fig4.jpeg")),
       plot = g1, height = 5, width = 6, units='in')
g2 <- g1 + coord_flip() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))
ggsave(file.path(paste0("data-raw/traits/isotopes/Meinzer1999_Table1_Xylem_Sap_deltaD_Fig4_flipped.jpeg")),
       plot = g2, height = 5, width = 6, units='in')

### More data from March1997 from Fig 5B (includes data in Fig 5A)
iso.2.raw <- read.csv("data-raw/traits/isotopes/Meinzer1999_Xylem_Sap_deltaD_March97_DBH_Fig5B.csv", na.strings = c("NA", ""), header = T, row.names = NULL, check.names = F)
iso.change <- read.csv("data-raw/traits/isotopes/Meinzer1999_Xylem_Sap_deltaD_by_delta_SapFlow_&_monthly_LeafFall_Fig7AB.csv", na.strings = c("NA", ""), header = T, row.names = NULL, check.names = F)

iso.2.raw <- iso.2.raw %>% mutate(source = "Meinzer et al.1999 Fig. 5B", location = "BCI")
iso.change <- iso.change %>% mutate(source = "Meinzer et al.1999 Fig. 7AB", location = "BCI")
head(iso.2.raw)
iso.2 <- iso.2.raw %>%
  group_by(sp, source) %>%
  summarise(se = sd(Xylem_sap_deltaD_permil, na.rm = TRUE)/sqrt(n()),
            Xylem_sap_deltaD_permil = mean(Xylem_sap_deltaD_permil, na.rm = TRUE),
            n = n(),
            DBH = mean(DBH, na.rm = TRUE))

iso.2 <- iso.2 %>% left_join(iso %>% dplyr::select(sp, Phenology), by = "sp") %>%
  left_join(iso.change %>% rename(change.source = source), by = "sp") %>%
  arrange(Phenology)
# View(iso.2)

#### Plotting soil isotopic variation from Meinzer et al and Jackson et al -----------
iso.soil.1 <- read.csv("data-raw/traits/isotopes/Meinzer1999_Fig2A_soil_deltaD_BCI_data_Mar_Apr_1997.csv", na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)
iso.soil.2 <- read.csv("data-raw/traits/isotopes/Oecologia 1995 Jackson_fig2_soil_deltaD_Gigante_data_dry_season_1992.csv", na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)

depth.m1 <- lm(depth ~ soil.deltaD, data = iso.soil.1)
summ.depth.m1 <- summary(depth.m1)
soil.label <- delta^2*H[soil]~"( \u2030)"
## expression(paste(delta^{2}, "H (\u2030)")))
# depth.m1.label1 = bquote(Depth*''==''*.(round(depth.m1$coefficients[1], 2)) ~''+''~
#                            .(round(depth.m1$coefficients[2], 2)) ~ delta^{2}*H)

depth.m1.label1 = paste0("Depth = ", round(depth.m1$coefficients[2], 2),
                             "*Soil delta 2H ", round(depth.m1$coefficients[1], 2))
depth.m1.label2 = paste0("R2 = ", round(summ.depth.m1$r.squared, 2),
                         "\np = ", round(summ.depth.m1$coefficients[2, 4], 4))
# label.df <- data.frame(x = -40, y= 0, label = deparse(depth.m1.label1))
formula <- 'y ~ x'
g1 <- ggplot(iso.soil.1, aes(y = depth, x = soil.deltaD)) +
  xlab(soil.label) + ylab("Depth (cm)") +
  scale_y_reverse(breaks = seq(from = 0, to = 100, by = 10)) +
  scale_x_continuous(breaks = seq(from = -60, to = -10, by = 5)) +
  coord_cartesian(xlim = c(-60, -10)) +
  theme(panel.grid.major = element_line(colour = "gray", size = 0.1)) +
  theme(text = element_text(size = 10)) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5, formula = formula) +
  geom_errorbarh(aes(xmax = soil.deltaD + se, xmin = soil.deltaD - se), size = 0.3, width = 0.2) +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 1, size = 2) +
  stat_poly_eq(aes(label = stat(eq.label)),
             npcx = 0.05, npcy = 0.95, rr.digits = 2,
             formula = formula, parse = TRUE, size = 3) +
  stat_poly_eq(aes(label = stat(adj.rr.label)),
               npcx = 0.05, npcy = 0.87, rr.digits = 2,
               formula = formula, parse = TRUE, size = 3) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = ifelse(p.value < 0.001, sprintf('italic(p)~"< 0.001"'),
                                     sprintf('italic(p)~"="~%.2f',stat(p.value)))),
                  parse = TRUE, npcx = 0.05, npcy = 0.75, size = 3)
ggsave(file.path(paste0("figures/UDI_confidence/Meinzer_etal_1999_Depth_vs_soil_deltaD.jpeg")),
       height = 3, width = 3.2, units='in')
ggsave(file.path(paste0("figures/UDI_confidence/Meinzer_etal_1999_Depth_vs_soil_deltaD.tiff")),
       height = 3, width = 3.2, units='in')

g1 +   geom_smooth(method = "loess", se = FALSE, span = 0.7, formula = formula) #+
  # geom_errorbarh(aes(xmax = soil.deltaD + se, xmin = soil.deltaD - se), size = 0.5)

depth.m2 <- lm(depth ~ soil.deltaD, data = iso.soil.2)
summ.depth.m2 <- summary(depth.m2)

depth.m2.label1 = paste0("Depth = ", round(depth.m2$coefficients[2], 2),
                         "*Soil delta 2H ", round(depth.m2$coefficients[1], 2))
depth.m2.label2 = paste0("R2 = ", round(summ.depth.m2$r.squared, 2),
                         "\np = ", round(summ.depth.m2$coefficients[2, 4], 4))
# but iso.soil.2 does not have se
g1 %+% iso.soil.2 + geom_smooth(method = "glm", se = FALSE, formula = formula) +
  geom_text(x = -35, y= 0, label = depth.m2.label1, color = "black", size = 3) +
  geom_text(x = -10, y= -85, label = depth.m2.label2, color = "black", size = 4) +
  xlim(-60,0) +  ylim(100, 0) + theme(text = element_text(size = 20))
ggsave(file.path(paste0("figures/UDI_confidence/Jackson_etal_1995_Depth_vs_soil_deltaD.jpeg")),
       height = 3, width = 3.5, units='in')


## adding Xylem sap's soil uptake depth ----
# how to incorporate se in this?
iso <- iso %>% mutate(depth = predict.lm(depth.m1, newdata = data.frame(soil.deltaD = Xylem_sap_deltaD_permil))/100,
                      depth.se = depth - predict.lm(depth.m1,
                                                    newdata = data.frame(soil.deltaD = Xylem_sap_deltaD_permil + SE))/100) # from cm to m
head(iso)

change.xylem.label <- expression('Change in '*delta^2*H[xylem]~"( \u2030)"*day^-1)

iso.2.raw <- iso.2.raw %>%
  left_join(iso %>% dplyr::select(sp, Phenology, deciduousness), by = "sp") %>%
  left_join(iso.change %>% rename(change.source = source) %>%
              dplyr::select(-location), by = "sp") %>%
  left_join(tlp %>% mutate(sp = as.character(sp)), by = "sp")

ggplot(iso.2.raw, aes(y =  Xylem_sap_deltaD_permil, x = DBH, color = tlp)) +
  geom_text(aes(y =  Xylem_sap_deltaD_permil + 1.5, x = DBH, label = sp),
            size = 2) +
  geom_point(aes(shape = deciduousness), size = 3) +
  scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
  scale_y_continuous(limits = c(-70, -18), breaks = seq(-70, -20, by = 10)) +
  xlab("DBH (cm)") + ylab(xylem.label) +
  guides(shape = guide_legend(title = "Deciduousness")) +
  ggtitle("Species March 1997 Xylem Sap deltaD Vs DBH") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))
ggsave(file.path(paste0("data-raw/traits/isotopes/Meinzer1999_Xylem_Sap_deltaD_March97_DBH_Fig5B.jpeg")), height = 5, width = 7, units ='in')

ggplot(iso.2.raw, aes(y =  delta_sapflow_percent_per_day, x = delta_xylem_sap_deltaD_permil_per_day, color = tlp)) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5, color = "black", formula = formula) +
  geom_text(aes(y =  delta_sapflow_percent_per_day + 0.05, label = sp),
            size = 4) +
  geom_point(aes(shape = deciduousness), size = 3) +
  guides(shape = guide_legend(title = "Deciduousness")) +
  scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
  scale_x_continuous(limits = c(-0.6, 0.1)) +
  xlab(change.xylem.label) + ylab("Change in daily sap flow (% per day)") +
  ggtitle("Sap flow Vs. Xylem Sap deltaD\nWithin species variation through 1997 dry season") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))
ggsave(file.path(paste0("data-raw/traits/isotopes/Meinzer1999_Xylem_Sap_deltaD_by_delta_SapFlow_&_monthly_LeafFall_Fig7A.jpeg")),
       height = 5, width = 7, units ='in')

ggplot(iso.2.raw, aes(y = SE_mean_monthly_leaf_fall, x = delta_xylem_sap_deltaD_permil_per_day, color = tlp)) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5, color = "black", formula = formula) +
  geom_text(aes(y = SE_mean_monthly_leaf_fall + 0.002, label = sp),
            size = 4) +
  geom_point(aes(shape = deciduousness), size = 3) +
  guides(shape = guide_legend(title = "Deciduousness")) +
  scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
  scale_x_continuous(limits = c(-0.6, 0.1)) +
  xlab(change.xylem.label) + ylab("SE of mean monthly leaf fall") +
  ggtitle("Leaf fall Vs. Xylem Sap deltaD\nWithin species variation through 1997 dry season") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))
ggsave(file.path(paste0("data-raw/traits/isotopes/Meinzer1999_Xylem_Sap_deltaD_by_delta_SapFlow_&_monthly_LeafFall_Fig7B.jpeg")), height = 5, width = 7, units ='in')


iso.3 <- read.csv("data-raw/traits/isotopes/Oecologia 1995 Jackson _Fig3_Fig4.csv", na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)
iso.3 <- iso.3 %>% arrange(Xylem_sap_deltaD_permil) %>%
  left_join(tlp %>% mutate(sp = as.character(sp)), by = "sp") %>%
  mutate(source = "Jackson et al. 1995") %>%
  left_join(deci, dplyr::select(-sp4), by = "sp")  %>%
  unite("deci_sp", deciduous, sp, remove = FALSE) %>%
  mutate(deci_sp.plot = factor(deci_sp, levels=unique(deci_sp[order(Xylem_sap_deltaD_permil)]), ordered=TRUE))


# iso.udi <- left_join(iso.3, udi, by = "sp")
# iso.sp <- unique(iso.udi$sp[!is.na(iso.udi$udi.med.rsq) & !is.na(iso.udi$Xylem_sap_deltaD_permil)])
# iso.udi <- iso.udi[order(iso.udiudi.med.rsq), ]

g3 <- ggplot(iso.3 %>% subset(!is.na(Xylem_sap_deltaD_permil)),
             aes(x = deci_sp.plot, y = Xylem_sap_deltaD_permil)) +
  geom_col(aes(fill = deciduousness), color = "gray") +
  guides(fill = guide_legend(title = "")) +
  theme(legend.position = "top") +
  ylab(xylem.label) + xlab("Species") +
  geom_errorbar(aes(ymax = Xylem_sap_deltaD_permil + se, ymin = Xylem_sap_deltaD_permil - se), width = 0.2, size = 0.5) +
  theme(axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(file.path(paste0("data-raw/traits/isotopes/Oecologia 1995 Jackson _Fig3_Fig4.jpeg")),
       plot = g3, height = 5, width = 6, units='in')
g4 <- g3 + coord_flip() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5), axis.text.y = element_text(size = 20))
ggsave(file.path(paste0("data-raw/traits/isotopes/Oecologia 1995 Jackson _Fig3_Fig4_flipped.jpeg")),
       plot = g4, height = 5, width = 6, units='in')

iso.1.3.join <- iso %>% subset(!is.na(Xylem_sap_deltaD_permil)) %>%
  dplyr::select(sp, Xylem_sap_deltaD_permil, deciduous, deciduousness, se, source) %>%
  bind_rows(iso.3 %>% dplyr::select(sp, Xylem_sap_deltaD_permil, deciduous, deciduousness, se, source)) %>%
  unite("deci_sp", deciduous, sp, remove = FALSE) %>%
  mutate(deci_sp.plot = factor(deci_sp, levels=unique(deci_sp[order(Xylem_sap_deltaD_permil)]), ordered=TRUE))

g5 <- g3 %+% subset(iso.1.3.join, !is.na(deciduousness)) +
  facet_wrap(. ~ source) +
  geom_col(aes(fill = deciduousness)) +
  coord_flip() +
  theme(text = element_text(size = 20)) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, size = 20),
        axis.text.y = element_text(size = 12))
ggsave(file.path(paste0("data-raw/traits/isotopes/Oecologia 1995 Jackson_Fig3_Fig4_& Meinzer 1999_Fig4_flipped.jpeg")),
       plot = g5, height = 5, width = 8, units='in')
save(iso.1.3.join, file = "data-raw/traits/isotopes/Oecologia 1995 Jackson_Fig3_Fig4_& Meinzer 1999_Fig4.Rdata")
load(file = "data-raw/traits/isotopes/Oecologia 1995 Jackson_Fig3_Fig4_& Meinzer 1999_Fig4.Rdata")
## Combine all isotopic records:
# View(iso)

iso.all.list <- list(iso, iso.2.raw, iso.3, iso.change)
names(iso.all.list) <- c(iso$source[1], iso.2.raw$source[1], iso.3$source[1], iso.change$source[1])

iso.all <- iso %>% dplyr::select(sp, location, source, Xylem_sap_deltaD_permil, se) %>%
  rbind(iso.3 %>% dplyr::select(sp, location, source, Xylem_sap_deltaD_permil, se)) %>%
  rbind(iso.2.raw %>% mutate(se = NA) %>% dplyr::select(sp, location, source, Xylem_sap_deltaD_permil, se)) %>%
  left_join(iso.change %>% dplyr::select(-source, -location), by = "sp")
View(iso.all)
save(iso.all, file = "results/all_isotopic_record.Rdata")

all.iso.sp <- unique(as.character(iso.all$sp))
all.iso.sp <- all.iso.sp[!is.na(all.iso.sp)]

save(all.iso.sp, file = "results/all_sp_with_isotopic_record.Rdata")

