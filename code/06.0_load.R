#-----------------------------------------------------
# Title: Loading data generated from
# script 05.0_Prepare_data_for_correlative analyses.R
# This data is used by 0.7.0_PhenoDemoTraitsPsi.R and 08.0_manuscript.rmd
# Author : Rutuja Chitra-Tarak
# Original date: Summer, 2020
#-----------------------------------------------------

rm(list=ls())

#*******************************************
####   Load Libraries, Prep for graphics, folders  ##
#*******************************************

if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, readxl, forcats, scales, data.table, ggpmisc, GGally)

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
## Load ELM-FATES best-fit parameters -----
#******************************************************
params.obj.top.few <- read.csv(file = file.path("results/2019-10-14_5000/params.obj.top.few_100.csv"), header = TRUE)
p.best <- params.obj.top.few[,1:13]
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
  dplyr::rename(form1 = GRWFRM1., sp = SP., SG100C_AVG = SG100C_AVG) %>% mutate(sp = tolower(sp))

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
# growth.type <- "stats"
# dbh.residuals <- "on"
# solar.residuals <- "off"
# growth.selection <- "size_class_predefined_cc_scaled"
# intervals <- 5
# growth.name <- load(file =  paste0("results/sp_size.", growth.type, "_growth_dbh.residuals_", dbh.residuals, "_", intervals, "_", growth.selection, ".Rdata"))
# # growth.name <- load(file =  paste0("results/sp_size.", growth.type, "_growth_dbh.residuals_", dbh.residuals, "_", intervals, "_", growth.selection, ".Rdata"))

# load(file.path("results/GLUEsetup_part2_BCI.RData"))
## growth rates when dbh.residuals = "on" are residuals from a dbh mixed effects model (for spp) of
## growth. A median residual for each sp_size is caluclated only when at least data from
# 3 trees are present across all census intervals.
# Medians within sp_size are then centered and scaled. {residual - E(residual)/sd(residual)}
growth.type <- "stats"
dbh.residuals <- "on"
solar.residuals <- "off"
growth.selection <- "size_class_predefined_cc_scaled"
intervals <- 5
growth.name <- load(file =  paste0("results/sp_size.", growth.type, "_growth_dbh.residuals_", dbh.residuals, "_", intervals, "_", growth.selection, ".Rdata"))
# growth.name <- load(file =  paste0("results/sp_size.", growth.type, "_growth_dbh.residuals_", dbh.residuals, "_", intervals, "_", growth.selection, ".Rdata"))

growth <- get(growth.name); rm(growth.name)

## No. of trees by sp for grpwth data
g.n <- lapply(growth[grep("large", names(growth))], as.data.frame) %>%
  bind_rows(.id = "sp_size") %>%
  separate(sp_size, c("sp", "size", sep = "_"), remove = FALSE, extra = "drop", fill = "right") %>%
  dplyr::select(-sp_size, -"_") %>%
  group_by(sp) %>%
  summarise(n = mean(trees), .groups = "drop_last")

load(file = file.path(results.folder, "grate.long_by_species-size_deciduousness.Rdata"))
load(file = file.path(results.folder, "mrate.long_by_species-size_deciduousness.Rdata"))
load(file = file.path(results.folder, "adult.mrate.long_by_species-size_deciduousness.Rdata"))

load(file = file.path("results/demo.sp_size.Rdata"))
load(file = file.path("results/demo.sp.Rdata"))
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
load(file = file.path(results.folder, "rain.man.stats.Rdata"))
#******************************************************
### Load K by Psi models -------
#******************************************************
load(file = file.path(results.folder, "k_by_psi.models.Rdata"))
load(file = file.path(results.folder, "gap.models.Rdata"))
load(file = file.path(results.folder, "bci.AB.Rdata"))
load(file = file.path(results.folder, "data.model.AB.Rdata"))
load(file = file.path(results.folder, "sp.soft.filled.Rdata"))
load(file = file.path(results.folder, "sp.exp.param.Rdata"))
#******************************************************
### Load Leaf Life-span and Cover estiamtes------
#******************************************************
load(file = file.path(results.folder, "sp.leaf_cover.Rdata"))
# load(file = file.path(results.folder, "sp.leaf_cover.mean.Rdata")) # without leaf gain
load(file = file.path(results.folder, "sp.leaf_cover.for.model.Rdata"))
load(file = file.path(results.folder, "gap.models.ll.Rdata"))

soil.depths <- unique(psi$depth)

#******************************************************
### Supporting data referenced in the manuscript
#******************************************************
# Number of species ------
richness.database <- read.csv("data-raw/global_tree_search_trees_1_4.csv", header = TRUE)
richness <- nrow(richness.database)
tot.rich.k <- signif(richness/1000, 2)
tropics.low <- signif(round(40000/richness*100, 0), 1)
tropics.high <- signif(round(53000/richness*100, 0), 1)

## Climate

clim <- clim.daily %>%
  subset(Year != 1984) %>% #since only one day available
  group_by(Year) %>%
  summarise(Precip = sum(Precip, na.rm = TRUE),
            pet.PM = sum(pet.PM, na.rm = TRUE), .groups = "drop_last")
range(clim$Year)

## LMA~lamina~ and LMA~
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

disc.cf <- as.numeric(round(gap.models$LMA.LAM.DISC$coefficients, 2))
disc.r2 <- round(summary(gap.models$LMA.LAM.DISC)$adj.r.squared, 2)
disc.n <- length(gap.models$LMA.LAM.DISC$residuals)
disc.p <- ifelse(broom::glance(gap.models$LMA.LAM.DISC)$p.value < 0.001,
                 paste0("< 0.001"), paste0("= ", signif(broom::glance(gap.models$LMA.LAM.DISC)$p.value, 2)))

lma.cf <- as.numeric(round(gap.models$LMA.LAM.LMA$coefficients, 2))
lma.r2 <- round(summary(gap.models$LMA.LAM.LMA)$adj.r.squared, 2)
lma.n <- length(gap.models$LMA.LAM.LMA$residuals)
lma.p <- ifelse(broom::glance(gap.models$LMA.LAM.LMA)$p.value < 0.001,
                 paste0("< 0.001"), paste0("= ", signif(broom::glance(gap.models$LMA.LAM.LMA)$p.value, 2)))

c.lma.gaps <- table(bci.AB$sp.LMA.sub)
c.wsg.gaps <- table(bci.AB$sp.WSG.sub)

leaf.cf <- as.numeric(round(gap.models$LMA.LAM.LEAF$coefficients, 2))
leaf.r2 <- round(summary(gap.models$LMA.LAM.LEAF)$adj.r.squared, 2)
leaf.n <- length(gap.models$LMA.LAM.LEAF$residuals)
leaf.p <- ifelse(broom::glance(gap.models$LMA.LAM.LEAF)$p.value < 0.001,
                paste0("< 0.001"), paste0("= ", signif(broom::glance(gap.models$LMA.LAM.LEAF)$p.value, 2)))


chave.cf <- as.numeric(round(gap.models$WSG.100.CHAVE$coefficients, 2))
chave.r2 <- round(summary(gap.models$WSG.100.CHAVE)$adj.r.squared, 2)
chave.n <- length(gap.models$WSG.100.CHAVE$residuals)
chave.p <- ifelse(broom::glance(gap.models$WSG.100.CHAVE)$p.value < 0.001,
                paste0("< 0.001"), paste0("= ", signif(broom::glance(gap.models$WSG.100.CHAVE)$p.value, 2)))

## A & B data models based on soft traits
# https://stats.stackexchange.com/questions/95939/how-to-interpret-coefficients-from-a-polynomial-model-fit
# models had raw = TRUE
acf <- vector() ; bcf <- vector()
for (i in 1: (length(k_by_psi.models$A.B.LMA.LAM$coefficients) - 1)) {
    acf[i] <- as.numeric(round(k_by_psi.models$A.B.LMA.LAM$coefficients[i], 2))
}
acf.6 <- as.numeric(round(k_by_psi.models$A.B.LMA.LAM$coefficients[6], 3))

for (i in 1: length(k_by_psi.models$B.WSG100.LMA$coefficients)) {
  if (i == 6){
    bcf[i] <- formatC(k_by_psi.models$B.WSG100.LMA$coefficients[i], format = "f")
  } else if (i == 4) {
    bcf[i] <- as.numeric(round(k_by_psi.models$B.WSG100.LMA$coefficients[i], 3))
  } else {
    bcf[i] <- as.numeric(round(k_by_psi.models$B.WSG100.LMA$coefficients[i], 2))
  }
}

afit.r2 <- signif(broom::glance(k_by_psi.models$A.B.LMA.LAM)$adj.r.squared, 2)
afit.p <- ifelse(broom::glance(k_by_psi.models$A.B.LMA.LAM)$p.value < 0.001,
                 paste0("< 0.001"), paste0("= ", signif(broom::glance(k_by_psi.models$A.B.LMA.LAM)$p.value, 2)))

bfit.r2 <- signif(broom::glance(k_by_psi.models$B.WSG100.LMA)$adj.r.squared, 2)
bfit.p <- ifelse(broom::glance(k_by_psi.models$B.WSG100.LMA)$p.value < 0.001,
       paste0("< 0.001"), paste0("= ", signif(broom::glance(k_by_psi.models$B.WSG100.LMA)$p.value, 2)))

a.lma.r2 <- signif(broom::glance(k_by_psi.models$A.LMA.LAM)$adj.r.squared, 2)
a.lma.p <- ifelse(broom::glance(k_by_psi.models$A.LMA.LAM)$p.value < 0.001,
                 paste0("< 0.001"), paste0("= ", signif(broom::glance(k_by_psi.models$A.LMA.LAM)$p.value, 2)))
a.wsg.r2 <- signif(broom::glance(k_by_psi.models$A.WSG100)$adj.r.squared, 2)
a.wsg.p <- ifelse(broom::glance(k_by_psi.models$A.WSG100)$p.value < 0.001,
                  paste0("< 0.001"), paste0("= ", signif(broom::glance(k_by_psi.models$A.WSG100)$p.value, 2)))
a.b.r2 <- signif(broom::glance(k_by_psi.models$A.B)$adj.r.squared, 2)
a.b.p <- ifelse(broom::glance(k_by_psi.models$A.B)$p.value < 0.001,
                  paste0("< 0.001"), paste0("= ", signif(broom::glance(k_by_psi.models$A.B)$p.value, 2)))

b.lma.r2 <- signif(broom::glance(k_by_psi.models$B.LMA.LAM)$adj.r.squared, 2)
b.lma.p <- ifelse(broom::glance(k_by_psi.models$B.LMA.LAM)$p.value < 0.001,
                  paste0("< 0.001"), paste0("= ", signif(broom::glance(k_by_psi.models$B.LMA.LAM)$p.value, 2)))
b.wsg.r2 <- signif(broom::glance(k_by_psi.models$B.WSG100)$adj.r.squared, 2)
b.wsg.p <- ifelse(broom::glance(k_by_psi.models$B.WSG100)$p.value < 0.001,
                  paste0("< 0.001"), paste0("= ", signif(broom::glance(k_by_psi.models$B.WSG100)$p.value, 2)))

b.wsg.cf <- as.numeric(round(k_by_psi.models$B.WSG100$coefficients, 2))
a.wsg.cf <- as.numeric(round(k_by_psi.models$A.WSG100$coefficients, 2))

vpd.cf <- as.numeric(round(gpp.models$eq.gpp.vpd$coefficients, 2))

#************************
### Leaf lifespan from LMA model
#************************

ll.lma.cf <- as.numeric(round(gap.models.ll$LMA.lifespan$coefficients, 2))
ll.lma.r2 <- round(summary(gap.models.ll$LMA.lifespan)$r.squared, 2)
ll.lma.n <- length(gap.models.ll$LMA.lifespan$residuals)
ll.lma.p <- ifelse(broom::glance(gap.models.ll$LMA.lifespan)$p.value < 0.001,
                  paste0("< 0.001"), paste0("= ", signif(broom::glance(gap.models.ll$LMA.lifespan)$p.value, 2)))

#******************************************************
### Tables-----
#******************************************************
hypo.table <-
  data.frame(
    Variable = c(
      "*K*~max,stem~",
      "$\\mathrm{\\Psi}$~88,stem~",
      "$\\mathrm{\\Psi}$~tlp~",
      "$\\mathrm{\\Psi}$~min~ - $\\mathrm{\\Psi}$~88,stem~"
    ),
    Deeper.ERD = c("Higher", "Less negative", "Less negative", "More negative"),
    Shallower.ERD = c("Lower", "More negative", "More negative", "Positive, or less negative")
  )
symbols.table <-
  data.frame(
    Symbol = c(
      "$\\mathrm{\\Psi}_{\\textrm{soil},z}$",
      "$\\mathrm{\\Psi}$~leaf~, $\\mathrm{\\Psi}$~stem~",
      "$\\mathrm{\\Psi}$~tlp~",
      "$\\mathrm{\\Psi}$~crit~ or $\\mathrm{\\Psi}$~50,leaf~",
      "$\\mathrm{\\Psi}$~88,stem~",
      "$\\mathrm{\\Psi}$~min~",
      "$\\mathrm{\\Psi}$~min~ - $\\mathrm{\\Psi}$~88,stem~",
      "*K*~leaf~",
      "*K*~max,leaf~",
      "*K*~max,stem~",
      "WSG",
      "LMA",
      "$\\mathrm{\\delta}$^2^H~xylem~"
    ),
    Definition = c(
      "Soil water potential at depth $z$",
      "Water potential of leaf, or stem, respectively",
      "Bulk leaf turgor loss point, the $\\mathrm{\\Psi}$~leaf~ where turgor potential = 0",
      "$\\mathrm{\\Psi}$~leaf~ at 50% loss of leaf conductance",
      "$\\mathrm{\\Psi}$~stem~ at 88% loss of stem conductivity",
      "Seasonal minimum water potential, the most negative $\\mathrm{\\Psi}$~leaf~ measured at midday in the dry season",
      "Above-ground hydraulic safety margin",
      "Leaf-area specific hydraulic conductance of leaf",
      "Maximum leaf area-specific hydraulic conductance of leaf",
      "Maximum stem area-specific hydraulic conductivity of stem",
      "Wood specific gravity",
      "Leaf mass per unit area",
      "$\\mathrm{\\delta}$^2^H of tree xylem sap"
    ),
    Units = c(
      "MPa",
      "MPa",
      "MPa",
      "MPa",
      "MPa",
      "MPa",
      "MPa",
      "mmol m^-2^ s^-1^ MPa^-1^",
      "mmol m^-2^ s^-1^ MPa^-1^",
      "kg m^-1^ s^-1^ MPa^-1^",
      "g cm^-3^",
      "g m^-2^",
      "‰"
    )

  )
param.table <-
  data.frame(
    Parameter = c(
      "fates_leaf_BB_slope",
      "fates_leaf_slatop",
      "fates_leaf_vcmax25top",
      "fates_roota_par",
      "fates_rootb_par",
      "fates_smpsc",
      "aveDTB",
      "fmax",
      "HKSAT_ADJ",
      "HKSAT",
      "fpi_max"
    ),
    Description = c(
      "stomatal slope parameter, as per Ball-Berry",
      "Specific Leaf Area (SLA) at the top of canopy, projected area basis",
      "Maximum carboxylation rate of Rubisco at $25^\\circ$C, canopy top",
      "ELM rooting distribution parameter",
      "ELM rooting distribution parameter",
      "Soil water potential at full stomatal closure",
      "Distance to bedrock",
      "The maximum fractional saturated area",
      "Adjusting factor for soil hydraulic conductivity",
      "Soil hydraulic conductivity profile",
      "Maximum interception fraction of precipitation"
    ),
    Global.minimum =
      c(
        "4",
        "0.0042",
        "22.7",
        "5.94",
        "0.217",
        "113000",
        "3",
        "0.01",
        "1",
        "0.007 for =< 12.5 cm; 5.56e-05 for >= 60 cm",
        "0.05"
      ),
    Global.maximum =
      c(
        "16",
        "0.0400",
        "92.5",
        "7.4",
        "1.956",
        "242000",
        "18",
        "0.80",
        "8",
        "0.014 for =< 12.5 cm; 8e-04 for >= 60 cm",
        "0.44"
      ),
    Best_fit.minimum = c(
      round(min(p.best$fates_leaf_BB_slope), 1),
      round(min(p.best$fates_leaf_slatop), 4),
      round(min(p.best$fates_leaf_vcmax25top), 1),
      round(min(p.best$fates_roota_par), 2),
      round(min(p.best$fates_rootb_par), 3),
      round(min(p.best$fates_smpsc), 0),
      round(min(p.best$aveDTB), 1),
      round(min(p.best$FMAX), 2),
      round(min(p.best$HKSAT_ADJ), 1),
      paste0(round(min(p.best$HKSAT_12.5), 3), " for =< 12.5 cm; ", round(min(p.best$HKSAT_60), 4), " for >= 60 cm"),
      round(min(p.best$fpi_max), 2)
    ),
    Best_fit.maximum = c(
      round(max(p.best$fates_leaf_BB_slope), 1),
      round(max(p.best$fates_leaf_slatop), 4),
      round(max(p.best$fates_leaf_vcmax25top), 1),
      round(max(p.best$fates_roota_par), 2),
      round(max(p.best$fates_rootb_par), 3),
      round(max(p.best$fates_smpsc), 0),
      round(max(p.best$aveDTB), 1),
      round(max(p.best$FMAX), 2),
      round(max(p.best$HKSAT_ADJ), 1),
      paste0(round(max(p.best$HKSAT_12.5), 3), " for =< 12.5 cm; ", round(max(p.best$HKSAT_60), 4), " for >= 60 cm"),
      round(max(p.best$fpi_max), 2)
    ),
    Units =
      c(
        "unitless",
        "m^2^ gC^-1^",
        "umol CO^2^ m^2^ s^-1^",
        "m^-1^",
        "m^-1^",
        "mm",
        "m",
        "unitless",
        "unitless",
        "mm s^-1^",
        "unitless"
      ),
    Rationale.and.references =
      c(
        "The range of Ball Berry parameter fitted across several studies as compiled in Table 2 of [@Medlyn:2012]",
        "Based on the observed range of Leaf Mass Area (LMA)--24.97 to 235.8 g m^-2^-- for individual tree variation across 51 tree species in Barro Colorado Island (this study), and assumption of 50% carbon content in biomass.",
        "Observed range of $V_{cmax}$ at  $25^\\circ$C for tropical tree species under relative radiation of 50% or greater, thereby excluding highly shaded leaves. Values derived with a conversion specific for CLM model are used. [@Ali:2015]",
        "Parameter a in Eq. 2 in [@Zeng:2001] that regulates the shape of the rooting profile. Range corresponds to this parameter specified for BATS (or IGBP) land cover classified Deciduous Broadleaf Trees (5.9 m^-1^ and Evergreen Broadleaf Trees  (7.4 m^-1^) as given in Table 1 of [@Zeng:2001].",
        "Parameter b in Eq. 2 in [@Zeng:2001] that regulates the depth of the rooting profile. Chosen range of *b* is derived using this equation so as to fit the observed range of rooting depth (d~r~) of 2 - 18 m for Tropical Deciduous Forest (mean ± Standard Error (SE); 3.7 ± 0.5, *n* = 5 trees; min = 2, max = 4.7)  and Tropical Evergreen Forest (mean $\\pm$ S.E.; 7.3 $\\pm$ 0.5, n = 3 trees and 3 communities; min = 2, max = 18 m) combined [@Canadell:1996]. Besides the direct observation of roots at 18 m included by [@Nepstad:1994] in Paragominas, eastern Amazonia that is included in the above study; in Tapajos, eastern Amazonia water extraction by roots was also inferred up to 18 m. [@Davidson:2011]",
        "Based on observed range of -1.13 to -2.42 MPa for leaf turgor loss point for individual tree variation across 49 BCI tree species. This study.",
        "Ben Turner; pers. comm.",
        "Empirical",
        "To account for high macroporosity and direct flow paths in tropical soils that is not accounted for by small soil core samples [@Broedel:2017; @Kinner:2004; @Tomasella:1998]",
        "For Conrad catchment observed median Ksat (95% CI) for 12.5 cm depth is 38.3 mm hr^-1^ (25.4 - 51.2, n = 75), while for 60 cm depth 0.7 mm hr^-1^ (0.2 - 1.2, n = 40) [@Kinner:2004]. For reference, a storm with 12.5 mm hr^-1^ rainfall intensity has a 0.2 probability of occurring in any given rainfall event.",
        "Based on throughfall data from [@Zimmermann:2010kzl] and precipitation data for BCI from STRI Physical Monitoring program, defined for precipitation events greater than 10 mm."
      )
  )

eddy.table <- data.frame(Variable = c("Rain",
                                              "Temperature flux",
                                              "H~2~O flux",
                                              "CO~2~ flux",
                                              "Turbulent intensity",
                                              "Fiction velocity",
                                              "H~2~O 4th moment",
                                              "CO~2~ 4th moment"),
                         Unit = c("mm", "C m^-2^ s^-1^", "mmol m^-2^ s^-1^", "mumol m^-2^ s^-1^",
                                  "-", "m s^-1^", "mmol m^-2^", "mumol m^-2^"),
# symbol = c(),
                          Criterion = c("During and 30 minutes after",
                                        "> -0.13 & < 0.51",
                                        "> -5 & < 30",
                                        "> -60 & < 20",
                                        "< 5",
                                        "> 0.15",
                                        "< 55",
                                        "< 55"))
