##-------------------------
## Date : 22 Jul 2018
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
#------------------
## map trees at BCI color coded by UDImean
# load tree locations data
load("data-raw/CTFScensuses/BCI.tree8.Rdata")
# load interval and working.iter
load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
load(file.path("results/GLUEsetup_part2_BCI.RData")) # has working.iter and growth and si matrix

intervals <- info$intervals
n.ensembles <- growth_by_si.info$n.ensembles
growth.type <- growth_by_si.info$growth.type
growth.selection <- growth_by_si.info$growth.selection
dbh.residuals <- growth_by_si.info$dbh.residuals
si.type <- growth_by_si.info$si.type
goodness.fit <- 0.3
soil.depths <- unique(info$root.param.long$depth)
dryseason <- "on"
root.selection <- "on"
rm(info); rm(growth_by_si.info)
##
# load species specific UDImean
sp.info.ud <- read.csv(file.path(paste0("results/iso.udi_cor", goodness.fit, "_", si.type, "_",
                                                  n.ensembles, "_", growth.type, "_", growth.selection, "_",
                                                  dbh.residuals, "_", intervals, "_id_dryseason_", dryseason,
                                                  "_root.selection_", root.selection, ".csv", sep = "")))

head(sp.info.ud)
sp.ud <- sp.info.ud %>% subset(tlplevel == "sp") %>% select(sp, udi.best)

# load(file = paste("results/Sp_UDI_UD_cm_>30cm.Rdata", sep = ""))

head(bci.tree8)
load("data-raw/CTFScensuses/bci_habitat.Rda")
head(bci_habitat)
bci.tree8 <- left_join(bci.tree8, sp.ud, by = "sp")
length(unique(sp.ud$sp))
# 57
length(unique(bci.tree8$sp))
# 328
alltrees <- length(unique(bci.tree8$treeID[bci.tree8$status != "D"])) # alive trees
# 221527 # including dead trees the number is 423617
UDtrees <- length(unique(bci.tree8$treeID[bci.tree8$sp %in% sp.ud$sp]))
# 114393
# proportion of alive trees that have an update depth estimate
UDtrees/alltrees*100
# 51%
trees.sub <- subset(bci.tree8, sp %in% sp.ud$sp)
ggplot(trees.sub, aes(x = gx, y = gy, colour = udi.best)) +
  geom_point(size = 1) +
  guides(colour = guide_legend(title = expression(UDI[mean]~(cm)))) +
  scale_color_gradient(name = "", high = "blue", low = "yellow") +
  my.theme + my.bg + my.adjust + theme(panel.grid.major.x = element_blank()) + theme(legend.text = element_text(size = 14, face = "plain"))
  theme(legend.position = "top")
ggsave(file.path(paste("figures/spatial/tree_UDI_map_#i74VKZ.jpeg")), height = 5, width = 10, units ='in')
string(1, 6)

## get elevation data for BCI
load("data-raw/CTFScensuses/CTFSElev_bci.rdata")
str(CTFSElev_bci) ## at 5 x 5 m quadrat
bci_elevation <- CTFSElev_bci[[1]]
# str(bci_habitat)
length(unique(bci_elevation$x))
length(unique(bci_elevation$y))
ele.wide <- dcast(bci_elevation, x ~ y)
par(pin = c(6, 3))
contour(as.matrix(ele.wide[,-1]))
ggsave(file.path(paste("figures/spatial/BCI_contour_map_#DmCpIo.jpeg")), height = 5, width = 10, units ='in')
string(1, 6)

#### testing if elevation (a proxy for distance from ground water) predicts UDI---------
## A hypothesis based on Fan et al 1017 PNAS
## adding quadrat specific elevation
# bci_elevation gives elevation for a 5 x 5 m grid (corners of a 5 x 5 m quadrat)
str(bci.tree8)
trees <- select(bci.tree8, c(treeID, tag, sp, quadrat, gx, gy, status, dbh, ba, udi.best))
trees$x5 <- as.numeric(cut(trees$gx, breaks = seq(0, 1000, by = 5))) ## somehow labels does not accept equal length as breaks, but complains that it needs equal numbers, so 0 removed, but now relabelling
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
bci_elevation$xy5 <- paste(bci_elevation$x, bci_elevation$y, sep = ",")
trees <- left_join(trees, select(bci_elevation, xy5, elev), by = "xy5")
head(trees)

ggplot(trees, aes(x = elev, y = udi.best)) +
  geom_point(size = 1) +
  my.theme + my.bg + my.adjust + theme(panel.grid.major.x = element_blank()) + theme(legend.text = element_text(size = 14, face = "plain")) +
  xlab("Quadrat Elevation (m)") + ylab(expression(UDI[mean]~(cm)))
ggsave(file.path(paste("figures/spatial/UDI_by_elevation_#cWPevW.jpeg")), height = 5, width = 10, units ='in')
string(1, 6)

# expect that trres on higher elevation should have higher udi.best
lm.ele <- lm(udi.best ~ elev, data = trees)
summary(lm.ele)
## So small effect size and in the opposite direction to what is expected
# Residuals:
#   Min      1Q  Median      3Q     Max
# -39.222 -14.437   5.081  10.163  41.003
#
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept) 73.754807   1.055887   69.85   <2e-16 ***
#   elev        -0.161176   0.007308  -22.06   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 19.99 on 114367 degrees of freedom
# (309248 observations deleted due to missingness)
# Multiple R-squared:  0.004235,	Adjusted R-squared:  0.004227
# F-statistic: 486.4 on 1 and 114367 DF,  p-value: < 2.2e-16
###-----------

# #### testing if habitat (a proxy for distance from ground water) is associated with UDI---------
## adding quadrat specific habitat type
# bci_habitat gives elevation for a 20 x 20 m grid (corners of a quadrat)
trees$x20 <- as.numeric(cut(trees$gx, breaks = seq(0, 1000, by = 20))) ## somehow labels does not accept equal length as breaks, but complains that it needs equal numbers, so 0 removed, but now relabelling
## (0, 20] is 1; to make (0, 20] as 0
trees$x20 <- (trees$x20 - 1)*20
summary(trees$x20)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#      20     240     500     504     760    1000     101
trees$y20 <- as.numeric(cut(trees$gy, breaks = seq(0, 500, by = 20))) ## somehow labels does not accept equal length as breaks, but complains that it needs equal numbers, so 0 removed, but now relabelling
trees$y20 <- (trees$y20 - 1)*20
trees$xy20 <- paste(trees$x20, trees$y20, sep = ",")
bci_habitat$xy20 <- paste(bci_habitat$x, bci_habitat$y, sep = ",")
trees <- left_join(trees, select(bci_habitat, xy20, habitat), by = "xy20")
head(trees)

ggplot(trees, aes(x = gx, y = gy, colour = habitat)) +
  geom_point(size = 1) +
  my.theme + my.bg + my.adjust + theme(panel.grid.major.x = element_blank()) + theme(legend.text = element_text(size = 14, face = "plain"))
ggsave(file.path(paste("figures/spatial/tree_on_habitat_#mzOHBf.jpeg")), height = 5, width = 10, units ='in')

lm.hab <- lm(udi.best ~ habitat, data = trees)
# Residuals:
#   Min      1Q  Median      3Q     Max
# -39.998 -14.677   4.741  11.435  43.537
#
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)         49.6526     0.1639 303.027  < 2e-16 ***
#   habitatlow_plateau   1.4115     0.1848   7.637 2.24e-14 ***
#   habitatmixed         1.9884     0.3201   6.212 5.23e-10 ***
#   habitatslope         0.2598     0.2005   1.296    0.195
#   habitatstream        3.5380     0.4103   8.622  < 2e-16 ***
#   habitatswamp         5.3651     0.4674  11.480  < 2e-16 ***
#   habitatyoung        -4.0614     0.3331 -12.192  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 19.98 on 114362 degrees of freedom
# (309248 observations deleted due to missingness)
# Multiple R-squared:  0.004852,	Adjusted R-squared:  0.004799
# F-statistic: 92.92 on 6 and 114362 DF,  p-value: < 2.2e-16

## so swamp and stream trees have (~5 cm)deeper roots than those on high-pleateau
## young patch has 4 cm shallower roots

bwplot(udi.best ~ habitat, data = trees, ylim = c(100,0))
ggsave(file.path(paste("figures/spatial/udi.best_by_habitat_bwplot_#JuWZDV.jpeg")), height = 5, width = 10, units ='in')

ggplot(subset(trees, habitat != "NA"),
       aes(x = habitat, y = udi.best, color = habitat)) +
  geom_violin() +
  guides(colour = "none") +
  scale_y_reverse( limits = c(100, 0)) +
  my.theme + my.bg + my.adjust + theme(panel.grid.major.x = element_blank()) +
  theme(axis.text = element_text(size = 14, face = "plain"))
ggsave(file.path(paste("figures/spatial/udi.best_by_habitat_boxplot_#JuWZDV.jpeg")), height = 5, width = 10, units ='in')

save(trees, file = file.path("data/BCI.tree8_ele_hab.Rdata"))
###-------------

## what are habitat associations of species?-------
## are species strongly associated with a habitat show specific; stream/swamp associated shallow udi.best, pleateu asociated deep?
## hab associatations by correlations
alive <- subset(trees, status != "A")
contingency.tab <- table(alive$sp, alive$habitat)
# ## if doesnt work; update packages, R and RStudio
# library(fgeo)
# install_github("forestgeo/fgeo.habitat")

