##--------------------
## Data by species
## Author: Rutuja
## Date April 14, 2020
##-------------------
rm(list = ls())
gc()
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, naniar)

### Demographic data -----
load("results/demo.sp.RData")
## keeping data only for adult
adult.demo.sp <- demo.sp %>% select(sp, grate.adult, mrate.adult, avg.abund) %>%
  subset(!(is.na(grate.adult) & is.na(mrate.adult))) %>% droplevels()
View(adult.demo.sp)
adult.demo.sp <- adult.demo.sp %>%
  mutate(mortrate.adult = ifelse(!is.na(mrate.adult), "1", NA))

### Growth data available for inversion- for adult trees----
load(file.path("results/4.1GLUEsetup_part2_BCI.RData"))
above10.sp.for.inversion <- growth_by_si.info$growth.meta %>% subset(size %in% c("medium"))
above30.sp.for.inversion <- growth_by_si.info$growth.meta %>% subset(size %in% c("large"))
growth.ts.for.inversion <- data.frame(sp = above10.sp.for.inversion$sp, growth.ts.above10cmDBH = "1") %>%
  full_join(data.frame(sp = above30.sp.for.inversion$sp, growth.ts.above30cmDBH = "1"), by= "sp")
growth.ts.for.inversion <- replace_with_na_all(growth.ts.for.inversion, condition = ~.x == "<NA>")

data.by.sp <- adult.demo.sp %>% select(-grate.adult, -mrate.adult, -grate.adult) %>%
  left_join(growth.ts.for.inversion, by = "sp")
### Canopy Species-----
## canopy or not?
# Metadata Panama Traits Definitions.doc:
# GRWFRM2$ - Values are as in GRWFRM1 except free-standing species can have multiple values if maximum size varies widely within Panama
# GRWFRM1$ - Values are Climber, HERB, S, U, M and T. S, U, M and T are free-standing species with maximum heights of 5, 10, 20 and > 30 m, respectively.

bci.traits <- read.csv("data-raw/traits/BCITRAITS_20101220.csv") %>%
  rename(form1 = GRWFRM1., form2 = GRWFRM2., sp = SP.) %>% mutate(sp = tolower(sp))

data.by.sp <- data.by.sp %>%
  left_join(bci.traits %>%
              mutate(max.height = recode(as.factor(form1), `S` = "5 m", `U` = "10 m", `M` = "20 m",  `T` = "> 30 m"),
                     canopy.sp = ifelse(form1 == "T", "1", NA)) %>%
              select(sp, max.height, form1, form2, canopy.sp), by = "sp")


### Xylem Deuterium Isotopic signature -----
iso.meinzer <- read.csv(file.path(paste0("data-raw/traits/isotopes/Meinzer1999_Table1_Xylem_Sap_deltaD_Fig4_sp_code.csv")),
         na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)
iso.jackson <- read.csv(file.path(paste0("data-raw/traits/isotopes/Oecologia 1995 Jackson _Fig3_Fig4.csv")),
                        na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)

data.by.sp <- data.by.sp %>%
  left_join(iso.meinzer %>%
              mutate(iso.Deut.Meinzer = ifelse(!is.na(Xylem_sap_deltaD_permil), "Yes;BCI", NA)) %>%
              select(sp, iso.Deut.Meinzer), by = "sp") %>%
  left_join(iso.jackson %>%
              mutate(iso.Deut.Jackson = ifelse(!is.na(Xylem_sap_deltaD_permil) & location == "Gigante", "Yes;Gigante",
                                             ifelse(!is.na(Xylem_sap_deltaD_permil) & location == "BCI", "Yes;BCI", NA))) %>%
              select(sp, iso.Deut.Jackson), by = "sp")

### Traits -----
traits.Kunert.indi <- read.csv("data-raw/traits/HydraulicTraits_Kunert/hydraulic_traits_panama_kunert.csv")
traits.Kunert.sp <- data.frame(sp = unique(traits.Kunert.indi$sp), traits.Kunert = "1")

tlp.Kunert.data <- read.csv("data-raw/traits/HydraulicTraits_Kunert/tlp_sp_mean.csv", na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)
lwp <- read_excel("data-raw/traits/2016ENSO_Panama_LWP_20170407181307/2016ENSO_Panama_LWP.xlsx",
                  sheet = "Panama_LWP")
hyd.Wolfe <- read.csv("data-raw/traits/HydraulicTraits_BrettWolfe/ht1_20200103.csv") # Brett's data
hyd.Wolfe.sp <- data.frame(sp = tolower(unique(hyd.Wolfe$sp)), hyd.traits.Wolfe = "1")

data.by.sp <- data.by.sp %>%
  left_join(traits.Kunert.sp, by = "sp") %>%
  left_join(tlp.Kunert.data %>%
              mutate(tlp.Kunert = ifelse(!is.na(tlp), "1", NA)) %>%
              select(sp, tlp.Kunert), by = "sp") %>%
  left_join(hyd.Wolfe.sp, by = "sp")

## Deciduousness-----ß
deci <- read_excel(file.path("data-raw/traits/nomenclature_R_20190524_Rready_Osvaldo Calderon & JoeWright_expert_opinion.xlsx"))

deci <- deci %>% mutate(sp = tolower(sp6)) %>%
  select(sp, deciduous) %>%
  mutate(deciduousness = recode(as.factor(deciduous), `DO` = "Obligate Deciduous", `DF` = "Facultative Deciduous",
                                       `DB` = "Brevideciduous",  `E` = "Evergreen")) %>%
  select(sp, deciduousness)
head(deci)

data.by.sp <- data.by.sp %>%
  left_join(deci, by = "sp")

## Inversely modelled Uptake Depth Index----

# Arrange by availability
data.by.sp <- data.by.sp %>%
  select(sp, canopy.sp,  max.height, deciduousness, avg.abund, growth.ts.above30cmDBH, growth.ts.above10cmDBH, everything()) %>%
  arrange(canopy.sp,growth.ts.above30cmDBH, hyd.traits.Wolfe,
          traits.Kunert, iso.Deut.Meinzer, iso.Deut.Jackson, mortrate.adult)
View(data.by.sp)
write.csv(data.by.sp, file.path(paste0("data-raw/traits/BCI_demo&traits_data_availability_by_sp.csv")), row.names = FALSE)

