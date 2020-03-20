##-------------------------
## Edited: Jan 23 2020
## Author : Rutuja
## Title : Spatial hydrological niche
##-------------------------
rm(list = ls())

if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, devtools, data.table, reshape2, fgeo, remotes)
# remotes::install_github("forestgeo/bciex")
# had to load XQuartz for Mac first
# install.packages("installr")
# available.packages("utils")
# packageUrl<- "https://cran.r-project.org/src/contrib/Archive/R.utils/R.utils_2.6.0.tar.gz"
# install.packages(packageUrl, repos=NULL, type='source')
# graphics info
theme_set(theme_bw())
theme_update(text = element_text(size = 14),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)
require(scales)
rev_sqrt_trans <- function() {
  scales::trans_new(
    name = "rev_sqrt",
    transform = function(x) -sqrt(abs(x)),
    inverse = function(x) x^2);
}
file.path.spatial <- file.path("figures", "spatial")
if(!dir.exists(file.path.spatial)) {dir.create(file.path.spatial)}

###------loading trees UDI for those with isotopic data--------

# load interval and working.iter
load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
load(file.path("results/4.1GLUEsetup_part2_BCI.RData")) # has working.iter and growth and si matrix

intervals <- info$intervals
n.ensembles <- growth_by_si.info$n.ensembles
growth.type <- growth_by_si.info$growth.type
growth.selection <- growth_by_si.info$growth.selection
dbh.residuals <- "on"#growth_by_si.info$dbh.residuals
si.type <- growth_by_si.info$si.type
goodness.fit <- 0.3
soil.depths <- unique(info$root.param.long$depth)
dryseason <- "on"
root.selection <- "on"
rm(info); rm(growth_by_si.info)
##
# load species specific UDImean
iso.udi <- read.csv(file.path(paste0("results/iso.udi_cor", goodness.fit, "_", si.type, "_",
                                        n.ensembles, "_", growth.type, "_", growth.selection, "_",
                                        dbh.residuals, "_", intervals, "_id_dryseason_", dryseason,
                                        "_root.selection_", root.selection, ".csv", sep = "")))

iso.udi.sub <- iso.udi %>% subset(tlplevel == "comm" & !is.na(Xylem_sap_deltaD_permil)) %>%
  select(sp, udi.best, udi, Xylem_sap_deltaD_permil, Phenology)
## loading UDI for splevel
load(file = paste("results/splevel/ds.bestfit_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_id_dryseason_", dryseason, "_root.selection_", root.selection, ".Rdata", sep = ""))
ds <- ds.bestfit

### Dristribution on Panama Isthumus
## Panama rainfall gradient preference


moist <- read.csv("data-raw/Condit_et_al_2013/TreeCommunityDrySeasonSpeciesResponse.csv")
moist <- moist %>% mutate(sp = paste0(tolower(str_sub(species, 1, 4)), str_sub(genus, 1, 2))) %>%
  rename(moist = Moist, moist.2 = Moist.2)

occur <- read.csv("data-raw/Condit_et_al_2013/TreeCommunityDrySeasonSpeciesOccurrence.csv")
deficit <- read.csv("data-raw/Condit_et_al_2013/TreeCommunityDrySeasonSurveySite.csv")

occur <- occur %>% mutate(sp = paste0(tolower(str_sub(species, 1, 4)), str_sub(genus, 1, 2))) %>%
  left_join(deficit %>% select(site, dry.season.moisture), by = "site")
head(occur)

###
#------------------
# load tree locations data
load("data-raw/CTFScensuses/BCI.tree8.Rdata")

# bci_elevation gives elevation for a 5 x 5 m grid (corners of a 5 x 5 m quadrat)
str(bci.tree8)
trees <- select(bci.tree8, c(treeID, tag, sp, quadrat, gx, gy, status, dbh, ba))
## just defining size class for the alive trees for now
## adding tree size class:
cut.breaks <- c(10, 50, 100, 300, max(trees$dbh, na.rm = T)) # about 45 trees per interval are > 1000 mm in dbh
cut.labels <- c("tiny", "small", "medium", "large")
trees <- trees %>% mutate(size.class = forcats::fct_explicit_na(cut(dbh, include.lowest = TRUE, breaks = cut.breaks,
                                                 labels = cut.labels, right = TRUE)),
                          sp_size = paste0(sp, size.class))
trees$x5 <- as.numeric(cut(trees$gx, breaks = seq(0.0, 1000.0, by = 5))) ## somehow labels does not accept equal length as breaks, but complains that it needs equal numbers, so 0 removed, but now relabelling
## (0, 5] is 1; to make (0, 5] as 0
trees$x5 <- (trees$x5 -1)*5
summary(trees$x5)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#     0.0   235.0   485.0   491.5   745.0   995.0     101
trees$y5 <- as.numeric(cut(trees$gy, breaks = seq(0, 500, by = 5))) ## somehow labels does not accept equal length as breaks, but complains that it needs equal numbers, so 0 removed, but now relabelling
trees$y5 <- (trees$y5 - 1)*5
# so tentatively quadrat is given by (trees$x, trees$y)
## so most South-West quadrat is (0,0)
## assigning elevation of its south-west corner (0, 0)
## (ideally a mean of four corners should be taken)
trees$xy5 <- paste(trees$x5, trees$y5, sep = ",")
## SWP map in mid-dry season for a regular and a drought year
swp.map.reg <- read.table("data-raw/Kupers_et_al/Output/BCI_SWP_map_mid_dry_season_regular.txt", header = TRUE)
swp.map.dry <- read.table("data-raw/Kupers_et_al/Output/BCI_SWP_map_mid_dry_season_drought.txt", header = TRUE)
swp.map <- swp.map.reg %>% mutate(swp.reg = swp) %>% full_join(swp.map.dry %>% mutate(swp.dry = swp) , by = c("x", "y")) %>%
  mutate(xy5 = paste(x - 2.5, y - 2.5, sep = ","))
head(swp.map)
trees <- left_join(trees, select(swp.map, xy5, swp.reg, swp.dry), by = "xy5")

hab.swp <- trees %>% #subset(size.class %in% c("medium", "large")) %>%
  group_by(sp) %>%
  summarise(med.swp.reg = median(swp.reg, na.rm = TRUE),
            med.swp.dry = median(swp.dry, na.rm = TRUE),
            sd.swp.reg = sd(swp.reg, na.rm = TRUE),
            sd.swp.dry = sd(swp.dry, na.rm = TRUE))
head(hab.swp)
write.csv(hab.swp, file.path("data-raw/sp.plot.hab.swp.csv"), row.names = FALSE)

hab.swp.n <- trees %>% #subset(size.class %in% c("medium", "large")) %>%
  group_by(sp, xy5) %>%
  summarise(n = n(), swp.reg = median(swp.reg, na.rm = TRUE), swp.dry = mean(swp.dry, na.rm = TRUE)) %>%
  ungroup(sp, xy5) %>% subset(!is.na(swp.dry))
hab.swp.tau <- hab.swp.n %>%
  subset(n > 0) %>%
  group_by(sp) %>%
  summarise(tau.swp.reg = as.numeric(cor.test(swp.reg, n, method = c("kendall"), exact = FALSE)$estimate),
            tau.swp.dry = as.numeric(cor.test(swp.dry, n, method = c("kendall"), exact = FALSE)$estimate))
#####---------
## adding quadrat specific elevation & habitat
load("data-raw/CTFScensuses/CTFSElev_bci.rdata")
str(CTFSElev_bci) ## at 5 x 5 m quadrat
bci_elevation <- CTFSElev_bci[[1]]
bci_elevation$xy5 <- paste(bci_elevation$x, bci_elevation$y, sep = ",")
trees <- left_join(trees, select(bci_elevation, xy5, elev), by = "xy5")
### Similarly for habitat
load("data-raw/CTFScensuses/bci_habitat.Rda")
head(bci_habitat) ## at 20 x 20 m quadrat
trees$x20 <- as.numeric(cut(trees$gx, breaks = seq(0, 1000, by = 20))) ## somehow labels does not accept equal length as breaks, but complains that it needs equal numbers, so 0 removed, but now relabelling
trees$x20 <- (trees$x20 -1)*20
trees$y20 <- as.numeric(cut(trees$gy, breaks = seq(0, 2000, by = 20))) ## somehow labels does not accept equal length as breaks, but complains that it needs equal numbers, so 0 removed, but now relabelling
trees$y20 <- (trees$y20 - 1)*20
bci_habitat$xy20 <- paste(bci_habitat$x, bci_habitat$y, sep = ",")
trees$xy20 <- paste(trees$x20, trees$y20, sep = ",")
trees <- left_join(trees, select(bci_habitat, xy20, habitat), by = "xy20")
head(trees)

alsebl <- trees %>% subset(sp == "alsebl")
View(alsebl)
## highest rank stable across size class, so not using group_by(sp, size.class, habitat)
level_key <- c(swamp = "10", stream = "8", slope = "6", young = "1", mixed = "1", low_plateau = "1", hi_plateau = "1")

hab.pc <- trees %>% #subset(size.class %in% c("medium", "large")) %>%
  group_by(habitat) %>%
  summarise(n.trees = n()) %>%
  mutate(plot.tree.ratio = round(n.trees/sum(n.trees), 2)) %>%
  arrange(desc(plot.tree.ratio)) %>%
  mutate(wetness = as.numeric(recode(habitat , !!!level_key))) %>%
  subset(!is.na(habitat))
hab.n <- trees %>% #subset(size.class %in% c("medium", "large")) %>%
  group_by(sp, habitat) %>%
  summarise(n.trees = n()) %>%
  mutate(tree.ratio = round(n.trees/sum(n.trees), 2)) %>%
  left_join(hab.pc %>% select(-n.trees), by = "habitat") %>%
  mutate(dif.tree.ratio = tree.ratio - plot.tree.ratio,
         max.hab = ifelse(dif.tree.ratio == max(dif.tree.ratio, na.rm = TRUE), habitat, NA)) %>%
  # left_join(tlp, by = "sp") %>% ##
  # ungroup(sp, habitat) %>%
  # mutate(sp = fct_reorder(sp, tlp, .desc = TRUE)) %>%
  arrange(sp, desc(dif.tree.ratio))
hab.n.iso <- hab.n %>%
  subset(sp %in% unique(iso.udi.sub$sp))
View(hab.n.iso)

ggplot(hab.n.iso %>% subset(!is.na(habitat)), aes(fill = dif.tree.ratio, y = sp, x = habitat)) +
  geom_tile() + scale_fill_viridis_c("Habitat Pref", trans = "reverse", option = "viridis")
ggsave(file.path(file.path.spatial, paste0("habitat_pref_for_iso_sp.jpeg")), height = 5, width = 8, units ='in')

hab.wet <- hab.n %>%
  group_by(sp) %>%
  summarise(hab.wet = sum(dif.tree.ratio*wetness, na.rm = TRUE)) %>%
  arrange(desc(hab.wet))
hab.wet.iso <- hab.wet %>%
  subset(sp %in% unique(iso.udi.sub$sp))
View(hab.wet.iso)

udi.wet.iso <- hab.wet.iso %>% left_join(iso.udi.sub, by = "sp")
ggplot(udi.wet.iso, aes(Xylem_sap_deltaD_permil, hab.wet)) +
  geom_point(aes(colour = Phenology)) +
  geom_text(aes(x =  Xylem_sap_deltaD_permil,
                y = hab.wet + diff(range(udi.wet.iso$hab.wet, na.rm = TRUE))/20, label = sp), size = 4)
ggplot(udi.wet.iso, aes(hab.wet, udi.best)) +
  geom_point(aes(colour = Phenology)) +
  geom_text(aes(y =  udi.best, x = hab.wet + diff(range(udi.wet.iso$hab.wet, na.rm = TRUE))/20,
                label = sp), size = 4) + scale_y_reverse() + scale_x_reverse()
# udi.wet <- ds %>% left_join(hab.wet, by = "sp")
# ggplot(udi.wet, aes(udi.best, hab.wet)) +
#   geom_point() +
#   geom_smooth(method = "loess")
## checkign against species tlp
traits.indi <- read.csv("data-raw/traits/hydraulic_traits_panama_kunert.csv") # Nobby's data
tlp <- traits.indi %>% group_by(sp) %>% select(-idividual, -ind_ID) %>%
  dplyr::summarise(tlp = mean(mean_TLP_Mpa, na.rm = TRUE))

udi.wet.iso <-  udi.wet.iso %>%
  left_join(tlp, by = "sp")
ggplot(udi.wet.iso, aes(tlp, hab.wet)) +
  geom_point(aes(colour = Phenology)) +
  geom_text(aes(x = tlp, y = hab.wet + diff(range(udi.wet.iso$hab.wet, na.rm = TRUE))/20, label = sp), size = 4) +
  scale_x_reverse()
ggsave(file.path(file.path.spatial, paste0("habitat.wetness_vs_tlp_for_iso_sp.jpeg")), height = 4, width = 6, units ='in')


udi.wet <- hab.wet %>% left_join(ds, by = "sp") %>%
  left_join(tlp, by = "sp")
ggplot(udi.wet %>% subset(!is.na(tlp)), aes(tlp, hab.wet)) +
  geom_point() +
  geom_smooth(method = "loess") +
  #geom_text(aes(x = tlp, y = hab.wet + diff(range(udi.wet$hab.wet, na.rm = TRUE))/20, label = sp), size = 4) +
  scale_x_reverse()

#
# m2 <- lm(udi ~ Xylem_sap_deltaD_permil, data = iso.udi.sub)
# rank.corr <- cor.test(x=iso.udi.sub$udi, y = -iso.udi.sub$Xylem_sap_deltaD_permil,
#                       method = 'spearman', exact = FALSE)
# summ.m2 <- summary(m2)
# m2.label = paste0("R-squared = ", round(summ.m2$r.squared, 3),
#                   ", p-val = ", round(summ.m2$coefficients[2, 4], 2),
#                   "\nSpearman's Rho = ", round(as.numeric(rank.corr$estimate), 2), ", p-val = ", round(as.numeric(rank.corr$p.value), 2))

##habitat_vs_tlp
trees <- trees %>%
  left_join(tlp, by = "sp")

ggplot(trees %>% subset(!is.na(habitat)), aes(y = tlp, x = habitat)) +
  geom_boxplot()
ggsave(file.path(file.path.spatial, paste0("habitat_vs_tlp.jpeg")), height = 4, width = 6, units ='in')

##habitat_vs_wsg
## adding wood specific gravity
load(file.path("data-raw/CTFScensuses/bci.spptable.rdata"))
trees <- trees %>%
  left_join(bci.spptable %>% select(sp, wsg), by = "sp")

ggplot(trees %>% subset(!is.na(habitat)), aes(y = wsg, x = habitat)) +
  geom_boxplot()
ggsave(file.path(file.path.spatial, paste0("habitat_vs_wsg.jpeg")), height = 4, width = 6, units ='in')
