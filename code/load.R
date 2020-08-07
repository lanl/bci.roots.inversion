
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

#******************************************************
### Tables
#******************************************************

symbols.table <-
  data.frame(
    Symbol = c(
      "$\\Psi_{soil, z}$",
      "$\\Psi_{leaf}$, $\\Psi_{stem}$",
      "$\\Psi_{tlp}$",
      "$K_{leaf}$, $K_{stem}$",
      "$K_{max, leaf}$, $K_{leaf, max}$",
      "$\\Psi_{88, leaf}$",
      "$\\Psi_{50, stem}$",
      "$\\Psi_{88, stem}$",
      "$\\Psi_{min}$"
    ),
    Definition = c(
      "Soil water potential at depth z",
      "Water potential of leaf or stem",
      "Bulk leaf turgor loss point, the $\\Psi_{leaf}$ where turgor potential = 0",
      "Hydraulic conductivity of leaf or stem",
      "Maximum area-specific hydraulic conductivity of leaf or stem",
      "$\\Psi_{leaf}$ at 88% loss of leaf conductivity",
      "$\\Psi_{stem}$ at 50% loss of stem conductivity",
      "$\\Psi_{stem}$ at 88% loss of stem conductivity",
      "Seasonal minimum water potential, the most negative $\\Psi_{leaf}$ measured at midday"
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
    Minimum =
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
    Maximum =
      c("16",
        "0.0400",
        "92.5",
        "7.4",
        "1.956",
        "242000",
        "18",
        "0.80",
        "8",
        "0.014 for =< 12.5 cm; 0.0003 for >= 60 cm",
        "0.44"
      ),
    Units =
      c(
        "unitless",
        "m^2^gC^-1^",
        "umol CO^2^m^2^s^-1^",
        "m^-1^",
        "m^-1^",
        "mm",
        "m",
        "unitless",
        "unitless",
        "mms^-1^",
        "unitless"
      ),
    `Method of determination with reference` =
      c(
        "The range of Ball Berry parameter fitted across several studies as compiled in Table 2 of [@Medlyn:2012]",
        "Based on the observed range of Leaf Mass Area (LMA)--24.97 to 235.8 gm^-2^-- for individual tree variation across 51 tree species in Barro Colorado Island (this study), and assumption of 50% carbon content in biomass.",
        "Observed range of $V_{cmax}$ at  $25^\\circ$C for tropical tree species under relative radiation of 50% or greater, thereby excluding highly shaded leaves. Values derived with a conversion specific for CLM model are used. [@Ali:2015]",
        "Parameter a in Eq. 2 in [@Zeng:2001] that regulates the shape of the rooting profile. Range corresponds to this parameter specified for BATS (or IGBP) land cover classified Deciduous Broadleaf Trees (5.9 m^-1^ and Evergreen Broadleaf Trees  (7.4 m^-1^) as given in Table 1 of [@Zeng:2001].",
        "Parameter b in Eq. 2 in [@Zeng:2001] that regulates the depth of the rooting profile. Chosen range of *b* is derived using this equation so as to fit the observed range of rooting depth (d~r~) of 2 - 18 m for Tropical Deciduous Forest (mean ± Standard Error (SE); 3.7 ± 0.5, *n* = 5 trees; min = 2, max = 4.7)  and Tropical Evergreen Forest (mean $\\pm$ S.E.; 7.3 $\\pm$ 0.5, n = 3 trees and 3 communities; min = 2, max = 18 m) combined [@Canadell:1996]. Besides the direct observation of roots at 18 m included by [@Nepstad:1994] in Paragominas, eastern Amazonia that is included in the above study; in Tapajos, eastern Amazonia water extraction by roots was also inferred up to 18 m. [@Davidson:2011]",
        "Based on observed range of -1.13 to -2.42 MPa for leaf turgor loss point for individual tree variation across 49 BCI tree species. This study.",
        "Ben Turner; pers. comm.",
        "Empirical",
        "To account for high macroporosity and direct flow paths in tropical soils that is not accounted for by small soil core samples [@Broedel:2017; @Kinner:2004; @Tomasella:1998]",
        "For Conrad catchment observed median Ksat (95% CI) for 12.5 cm depth is 38.3 mmhr^-1^ (25.4 - 51.2, n = 75), while for 60 cm depth 0.7 mmhr^-1^ (0.2 - 1.2, n = 40) [@Kinner:2004]. For reference, a storm with 12.5 mmhr^-1^ rainfall intensity has a 0.2 probability of occurring in any given rainfall event.",
        "Based on throughfall data from [@Zimmermann:2010kzl] and precipitation data for BCI from STRI Physical Monitoring program, defined for precipitation events greater than 10 mm."
        )
        )

