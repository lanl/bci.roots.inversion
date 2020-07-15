
#*******************************************
####   Load Libraries, Prep for graphics, folders  ##
#*******************************************

if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, readxl, forcats, scales, data.table, ggpmisc)

# graphics info
theme_set(theme_bw())
theme_update(text = element_text(size = 14),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)
#******************************************************
## folders
#******************************************************
figures.folder <- paste0("figures/PhenoDemoTraitsPsi")
if(!dir.exists(file.path(figures.folder))) {dir.create(file.path(figures.folder))}
results.folder <- paste0("results/PhenoDemoTraitsPsi")
if(!dir.exists(file.path(results.folder))) {dir.create(file.path(results.folder))}
#******************************************************
## Load Deciduousness-----
#******************************************************
load(file = file.path(results.folder, "deci.prepped.Rdata"))

#******************************************************
### Load LWP -----
#******************************************************
load(file = file.path(results.folder, "lwp.min.Rdata"))
load(file = file.path(results.folder, "lwp.diff.Rdata"))

#******************************************************
## Load BCI traits---
#******************************************************
bci.traits <- read.csv("data-raw/traits/BCITRAITS_20101220.csv") %>%
  rename(form1 = GRWFRM1., sp = SP., SG100C_AVG = SG100C_AVG) %>% mutate(sp = tolower(sp))

#******************************************************
### Load hydraulic traits -----
#******************************************************
load(file = file.path(results.folder, "hyd.traits.all.RData"))
load(file = file.path(results.folder, "kunert.traits.all.RData"))
load(file = file.path(results.folder, "hyd.traits.key.long.RData"))
load(file = file.path(results.folder, "kunert.traits.key.long.RData"))

#******************************************************
## Load Isotopic data-----
#******************************************************
load(file = "data-raw/traits/isotopes/Oecologia 1995 Jackson_Fig3_Fig4_& Meinzer 1999_Fig4.Rdata")
leafless_mar.apr <- read.csv("data-raw/traits/isotopes/Meinzer_1999_isotope_sp_leafless_in_mar_april.csv")
# load(file = "results/all_isotopic_record.Rdata")
iso.2.raw <- read.csv("data-raw/traits/isotopes/Meinzer1999_Xylem_Sap_deltaD_March97_DBH_Fig5B.csv", na.strings = c("NA", ""), header = T, row.names = NULL, check.names = F)

#******************************************************
## Load Demographic data----
#******************************************************
# load(file.path("results/GLUEsetup_part2_BCI.RData"))
## growth rates when dbh.residuals = "on" are residuals from a dbh mixed effects model (for spp) of
## growth. A median residual for each sp_size is caluclated only when at least data from
# 3 trees are present across all census intervals.
# Medians within sp_size are then centered and scaled. {residual - E(residual)/sd(residual)}
growth.type <- "med"
dbh.residuals <- "on"
solar.residuals <- "off"
growth.selection <- "size_class_predefined_cc_scaled"
intervals <- 5
growth.name <- load(file =  paste0("results/sp_size.", growth.type, "_growth_dbh.residuals_", dbh.residuals, "_ci_", intervals, "_", growth.selection, ".Rdata"))
# growth.name <- load(file =  paste0("results/sp_size.", growth.type, "_growth_dbh.residuals_", dbh.residuals, "_", intervals, "_", growth.selection, ".Rdata"))

growth <- get(growth.name); rm(growth.name)

load(file = file.path(results.folder, "grate.long_by_species-size_deciduousness.Rdata"))
load(file = file.path(results.folder, "mrate.long_by_species-size_deciduousness.Rdata"))
load(file = file.path(results.folder, "adult.mrate.long_by_species-size_deciduousness.Rdata"))

load(file = file.path("results/demo.sp_size.Rdata"))

#******************************************************
### Load Psi from ELM-FATES-------
#******************************************************
load(file = file.path(results.folder, "psi.prepped.Rdata"))
load(file = file.path(results.folder, "psi.interpolated.depths_layers_combined.Rdata"))
#******************************************************
### Load climate data-------
#******************************************************
load(file = file.path(results.folder, "clim.daily_with_pet.PM.Rdata"))
load(file = file.path(results.folder, "gpp.models.Rdata"))

#******************************************************
### Load K by Psi models -------
#******************************************************
load(file = file.path(results.folder, "k_by_psi.models.Rdata"))
load(file = file.path(results.folder, "data.model.AB.Rdata"))
load(file = file.path(results.folder, "sp.soft.filled.Rdata"))

#******************************************************
### Load Leaf Cohort tracking data from the crane sites------
#******************************************************
load(file = file.path(results.folder, "sp.leaf_cover.Rdata"))
load(file = file.path(results.folder, "sp.leaf_cover.for.model.Rdata"))

load(file = file.path(results.folder, "cohort.Rdata"))
load(file = file.path(results.folder, "coh.sp.summ.Rdata"))
