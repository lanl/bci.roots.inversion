#______________________________
# Title: Phenology, demography, water availability
# Author : Rutuja Chitra-Tarak
# Original date: April 25, 2020
#_______________________________

## 1. How is species phenology located on the fast-slow continuum?
## 1.a Are deciduous species fast while evergreen species slow?
## 1.b demographically fast == high growth and mortality rates; slow == low growth and mortality rates
## 1.c fast traits high SLA, high Kmax, high TLP, deep roots
## 1.d fast species located on greater resource environments:
##        wetter sites within 50-ha plot, wet distributed along Panama gradient

#****************************
# ##   Libraries ######
#****************************

rm(list=ls())

if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, readxl, forcats, agricolae, gridExtra,
               scales, GGally, ggpmisc, Evapotranspiration, data.table)
# graphics info
theme_set(theme_bw())
theme_update(text = element_text(size = 14),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)

n.threshold = 50
figures.folder <- paste0("figures/PhenoDemoTraitsPsi")
if(!dir.exists(file.path(figures.folder))) {dir.create(file.path(figures.folder))}
results.folder <- paste0("results/PhenoDemoTraitsPsi")
if(!dir.exists(file.path(results.folder))) {dir.create(file.path(results.folder))}

#****************************
##   Custom Functions  ######
#****************************
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
indicator <- function(x, I.threshold, greater.than = TRUE) {
  if(greater.than == TRUE) {
    result <- ifelse((x > I.threshold), 1, 0)
  } else {
    result <- ifelse((x < I.threshold), 1, 0)
  }
  return(result)
}

rev_sqrt_trans <- function() {
  scales::trans_new(
    name = "rev_sqrt",
    transform = function(x) -sqrt(abs(x)),
    inverse = function(x) x^2);
}

## to calculate tlp based psi thresholds

# X1=sum(min(0, psi-psi_threshold)*PET), max/mean
# X2=sum(I(psi<psi_threshold)), max/mean
# X3=sum(min(0, psi-psi_threshold)), max/mean
# X4=sum(I(psi<psi_threshold)*PET),max/mean

psi.corr.fun.ls <- list(
  "gr.Psi" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := psi - df$psi_threshold][
        psi.mod < 0, psi.mod := 0][
        , keyby = .(depth, interval),
        .(gfac = mean(psi.mod, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    },
  "gr.Psi.Rad.VPD" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := psi - df$psi_threshold][
        psi.mod < 0, psi.mod := 0][
        , keyby = .(depth, interval),
        .(gfac = mean(psi.mod*std.Rs*std.VPD, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    } ,
  "gr.Psi.Rad.PET" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := psi - df$psi_threshold][
        psi.mod < 0, psi.mod := 0][
        , keyby = .(depth, interval),
        .(gfac = mean(psi.mod*std.Rs*std.pet.PM, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    } ,
  "mr.Psi" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := psi - df$psi_threshold][
                                   psi.mod > 0, psi.mod := 0][
                                   , keyby = .(depth, interval),
                                   .(gfac = mean(psi.mod, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    } ,
  "mr.Psi.PET" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := psi - df$psi_threshold][
          psi.mod > 0, psi.mod := 0][
            , keyby = .(depth, interval),
            .(gfac = mean(psi.mod*std.pet.PM, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    } ,
  "mr.Psi.VPD" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := psi - df$psi_threshold][
          psi.mod > 0, psi.mod := 0][
            , keyby = .(depth, interval),
            .(gfac = mean(psi.mod*std.VPD, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    } ,
  "mr.Psi.I" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := indicator(psi, df$psi_threshold, greater.than = FALSE)][
                                   , keyby = .(depth, interval),
                                   .(gfac = mean(psi.mod, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    } ,
  "mr.Psi.PET.I" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := indicator(psi, df$psi_threshold, greater.than = FALSE)][
                                   , keyby = .(depth, interval),
                                   .(gfac = mean(psi.mod*std.pet.PM, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    },
  "mr.Psi.VPD.I" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := indicator(psi, df$psi_threshold, greater.than = FALSE)][
          , keyby = .(depth, interval),
          .(gfac = mean(psi.mod*std.VPD, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    }
)
# Originally written by Sean
get.ts.lk <- function(df) {
  demo.ts <- df$demo.rate
  k <- length(demo.ts)
  psi.dat <- df[, grep("[0-9]", names(df))]
  loglk <- apply(psi.dat, 2, function(x) sum(dnorm(demo.ts, mean = scale(x), 1, log = TRUE)))
  R2 <- apply(psi.dat, 2,
              function(x)  summary(lm(demo.ts ~ x))$r.squared)
  corr <- apply(psi.dat, 2,
                function(x)  cor(demo.ts, x, method = "pearson", use = "complete.obs"))
  AIC.psi <- -2 * loglk + 2*k
  aic.sorted <- sort(AIC.psi, decreasing = FALSE, index.return = TRUE)
  depth <- names(aic.sorted$x)[which(aic.sorted$x <= aic.sorted$x[1] + 2)]
  return(list(depth = depth,
              aic.scores = aic.sorted$x[which(aic.sorted$x <= aic.sorted$x[1] + 2)],
              ml.scores = loglk[aic.sorted$ix[which(aic.sorted$x <= aic.sorted$x[1] + 2)]],
              R2 = R2[aic.sorted$ix[which(aic.sorted$x <= aic.sorted$x[1] + 2)]],
              corr = corr[aic.sorted$ix[which(aic.sorted$x <= aic.sorted$x[1] + 2)]]))
}

get.ml.max <- function(result) {
  result$ml.scores[1]
}

get.ml.depth.rsq <- function(result) {
  data.frame(depth = as.numeric(as.character(result$depth[which(result$corr >= 0)])),
             R2 = result$R2[which(result$corr >= 0)],
             corr = result$corr[which(result$corr >= 0)])
}

# Defines function to color according to correlation

cor_func <- function(data, mapping, method, symbol, ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)

  corr <- cor(x, y, method=method, use='complete.obs')
  colFn <- colorRampPalette(c("brown1", "white", "dodgerblue"),
                            interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length = 100))]

  ggally_text(
    label = paste(symbol, as.character(round(corr, 2))),
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    color = 'black',
    ...
  ) + #removed theme_void()
    theme(panel.background = element_rect(fill = fill))
}

ggpairs.theme <- theme(strip.background = element_rect(fill = "white"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 1, hjust = 1))

#******************************************************
## Load Deciduousness-----
#******************************************************

deci <- read_excel(file.path("data-raw/traits/nomenclature_R_20190524_Rready_Osvaldo Calderon & JoeWright_expert_opinion.xlsx"))
deci.level_key <- c("Evg" = "1", "DF" = "2", "DB" = "3", "DO" = "4") #c(a = "apple", b = "banana", c = "carrot")

deci <- deci %>% mutate(sp = tolower(sp6)) %>%
  select(sp4, sp, deciduous) %>%
  subset(deciduous %in% c("E", "DF", "DO", "DB")) %>% droplevels() %>%
  mutate(deciduousness = recode_factor(as.factor(deciduous), `E` = "Evergreen", `DB` = "Brevideciduous",
                                       `DF` = "Facultative Deciduous", `DO` = "Obligate Deciduous"),
         deciduousness.label = recode_factor(as.factor(deciduous), `E` = "E = Evergreen", `DB` = "DB = Brevideciduous",
                                       `DF` = "DF = Facultative Deciduous", `DO` = "DO = Obligate Deciduous")) %>%
  transform(deciduous = factor(deciduous,
                               levels = c("E", "DB", "DF", "DO"), ordered = TRUE),
            deciduousness = factor(deciduousness,
                                   levels = c("Evergreen", "Brevideciduous",
                                              "Facultative Deciduous", "Obligate Deciduous"), ordered = TRUE),
            deciduousness.label = factor(deciduousness.label,
                                   levels = c("E = Evergreen", "DB = Brevideciduous",
                                              "DF = Facultative Deciduous", "DO = Obligate Deciduous"), ordered = TRUE)) %>%
  mutate(DeciLvl = as.numeric(recode_factor(deciduous, !!!deci.level_key))) %>%
  select(sp4, sp, deciduous, deciduousness, deciduousness.label, DeciLvl)

head(deci)
# deci.2 <- read.csv("data-raw/traits/HydraulicTraits_Kunert/deciduous_species_Meakem.csv")
# deci.2 <- deci.2 %>% mutate(sp = as.character(Species.code), deciduousness = as.character(Deciduousness)) %>%
#   select(sp, deciduousness)

#******************************************************
### Load LWP -----
#******************************************************

lwp <- read_excel("data-raw/traits/2016ENSO_Panama_LWP_20170407181307/2016ENSO_Panama_LWP.xlsx",
                  sheet = "Panama_LWP")
lwp <- lwp %>%
  mutate(month = format(Date, "%b"),
         date.chr = as.character(Date),
         date = as.Date(Date),
         LWP_date.time.chr = paste(date.chr, LWP_coll_time),
         LWP_date.time = as.POSIXct(LWP_date.time.chr, tz = "America/Panama", format = "%Y-%m-%d %H:%M")) %>%
  subset(LWP_bar != -9999) %>% select(-date.chr, -LWP_date.time.chr, -Date)
head(lwp)
## same sp is not measured at two different locations
lwp.all <- lwp %>%
  mutate(sp = tolower(Species),
         # Brett's e-mail: ALBIED should be ALBIAD for Albizia adinocephala. We were using ALBIED mistakenly at the beginning
         sp = ifelse(sp == "albied", "albiad", sp)) %>%
  select(-Species) %>%
  rename(location = Location) %>%
  group_by(location, date, sp, Gasex_typeOrSequence, LWP_coll_time) %>%
  summarise(lwp.min = -min(LWP_bar/10, na.rm = TRUE)) %>% ## in MPa
  mutate(LWP_coll_time = as.numeric(LWP_coll_time)) %>%
  unite("sp_date", sp, date, remove = FALSE) %>%
  data.frame() %>%
  transform(location = factor(location, levels = c("PA-PNM", "PA-BCI", "PA-SLZ"), ordered = TRUE)) %>%
  rename(time = Gasex_typeOrSequence) %>%
  unite("sp_date_time", sp, date, time, remove = FALSE)

str(lwp.all)

## Most of the trees have multiple records per time step
lwp.all.indi <- lwp.all %>%
  group_by(location, sp, date, time, sp_date) %>%
  summarise(lwp.min = mean(lwp.min, na.rm = TRUE))

# First getting those species-(location)-wise dates and date-times with min.LWP
lwp.min.diurnal.dates <- lwp.all.indi %>%
  subset(!time %in% "Pre-dawn LWP") %>%
  group_by(sp, location, sp_date) %>%
  # Location not really needed as this data set does not have same sp measured at two different locations
  filter(lwp.min == min(lwp.min, na.rm = TRUE)) %>%
  unite("sp_date_time", sp, date, time, remove = FALSE)

## getting predawn LWPs for the same day as that of those minimum diurnal LWPs
## Retain multiple records per tree to calculate sd
lwp.min.diurnal <- lwp.all %>%
  subset(sp_date_time %in% lwp.min.diurnal.dates$sp_date_time) %>%
  group_by(sp, location) %>%
  # Location (sp_date_location) not really needed as this dat asets does not have same sp measured at two different locations
  bind_rows(lwp.all %>% subset(time == "Pre-dawn LWP" &
                                 # Location (sp_date_location) not really needed as this dat asets does not have same sp measured at two different locations
                                 sp_date %in% lwp.min.diurnal.dates$sp_date)) %>%
  mutate(time = ifelse(time == "Pre-dawn LWP", "Predawn", "Diurnal")) %>%
  mutate(time = recode_factor(time, `predawn` = "Predawn", `diurnal` = "Diurnal"), ordered = TRUE) %>%
  group_by(sp, time, location) %>%
  summarise(lwp.se = sd(lwp.min, na.rm = TRUE)/sqrt(n()),
            ## don't change this sequence
            lwp.min = mean(lwp.min, na.rm = TRUE)) %>%
  ungroup(sp, time, location) %>%
  arrange(sp, location, time) %>%
  group_by(sp, location) %>%
  mutate(lwp.diff = lwp.min - lag(lwp.min)) %>%
  ungroup(sp, location)

## also, to calculate difference between preadawn and min.dirunal measurement across all measurements
## first calculate minimun dirunal measurement for each day
# first get which time had the minimum water potential on each measurement day by species
lwp.diff.filter <- lwp.all %>%
  subset(!time %in% "Pre-dawn LWP") %>%
  group_by(sp_date) %>%
  # Location not really needed as this data set does not have same sp measured at two different locations
  slice(which.min(lwp.min)) %>% # selects only the match for the first  min
  arrange(sp, date, time)

lwp.diff <- lwp.all %>%
  subset(sp_date_time %in% lwp.diff.filter$sp_date_time) %>%
  bind_rows(lwp.all %>%
              subset(time == "Pre-dawn LWP")) %>%
  mutate(time = ifelse(time == "Pre-dawn LWP", "Predawn", "Diurnal")) %>%
  mutate(time = recode_factor(time, `predawn` = "Predawn", `diurnal` = "Diurnal"), ordered = TRUE) %>%
  arrange(sp, date, time) %>%
  group_by(sp, date, time, location) %>%
  summarise(lwp.min = mean(lwp.min, na.rm = TRUE)) %>%
  ungroup(sp, date, time, location) %>%
  group_by(sp, date, location) %>%
  arrange(sp, date, time) %>%
  mutate(lwp.diff = lwp.min - lag(lwp.min)) %>%
  ungroup(sp, date, location) %>%
  arrange(sp, date, time) %>%
  left_join(deci, by = "sp") %>%
  unite("sp_date", sp, date, remove = FALSE) %>%
  unite("deci_sp", deciduous, sp, remove = FALSE)
View(lwp.diff)
lwp.min.wide <- lwp.min.diurnal %>% pivot_wider(names_from = time, values_from = c(lwp.min, lwp.se, lwp.diff)) %>%
  mutate(lwp.diff = lwp.diff_Diurnal) %>% select(-lwp.diff_Diurnal, -lwp.diff_Predawn)

lwp.min <- lwp.min.diurnal %>%
  left_join(deci, by = "sp")

save(lwp.diff, file = file.path(results.folder, "lwp.diff.Rdata"))
save(lwp.min, file = file.path(results.folder, "lwp.min.Rdata"))
save(lwp.min.wide, file = file.path(results.folder, "lwp.min.wide.Rdata"))
load(file = file.path(results.folder, "lwp.min.Rdata"))

#******************************************************
## Load Panama rainfall gradient preference------
#******************************************************

moist.pref <- read.csv("data-raw/Condit_et_al_2013/TreeCommunityDrySeasonSpeciesResponse.csv")
moist.pref <- moist.pref %>% mutate(sp = paste0(tolower(str_sub(species, 1, 4)), str_sub(genus, 1, 2))) %>%
  rename(moist.pref = Moist, moist.pref.2 = Moist.2) %>% select(sp, moist.pref, moist.pref.2)
moist.pref <- moist.pref %>%
  group_by(sp) %>%
  summarise_all(list(~mean(., na.rm = TRUE)))
hab.swp <- read.csv(file.path("data-raw/sp.plot.hab.swp.csv"))
sp.hab <- moist.pref %>% full_join(hab.swp, by = "sp") %>%
  rename(Panama.moist.pref = moist.pref, Panama.moist.pref.2 = moist.pref.2) %>%
  mutate(Plot.swp.pref = med.swp.reg - mean(med.swp.reg, na.rm = TRUE),
         Plot.swp.ENSO = med.swp.dry - mean(med.swp.dry, na.rm = TRUE)) %>%
  select(-med.swp.reg, -med.swp.dry, -sd.swp.reg, -sd.swp.dry, -Panama.moist.pref.2, -Plot.swp.ENSO)

#******************************************************
## Load Demographic data----
#******************************************************

load("results/demo.sp.RData")
load("results/demo.sp_size.RData")
load("results/mrate.long.RData")
load(file.path("results/GLUEsetup_part2_BCI.RData"))
## growth rates when dbh.residuals = "on" are residuals from a dbh mixed effects model (for spp) of
## growth. A median residual for each sp_size is caluclated only when at least data from
# 3 trees are present across all census intervals.
# Medians within sp_size are then centered and scaled. {residual - E(residual)/sd(residual)}

if(growth_by_si.info$dbh.residuals == "on"){
  growth <- growth_by_si.info$growth
}

grate.long <- dplyr::bind_rows(growth, .id = 'sp_size')

census.years <- c(1982, 1985, 1990, 1995, 2000, 2005, 2010, 2015)

size.level_key <- c(tiny = "Tiny (1-5cm)", small = "Small (5-10cm)",
                    medium = "Medium (10-30cm)", large = "Large (>= 30cm)")
census.years.short <- format(strptime(census.years, "%Y"), "%y")

mrate.long <- mrate.long %>%
  separate(sp_size, c("sp", "size", sep = "_"), remove = FALSE, extra = "drop", fill = "right") %>%
  select(-"_") %>%
  left_join(deci, by = "sp") %>%
  mutate(censusint.m = recode(census, `1985` = "1982-85", `1990` = "1985-90", `1995` = "1990-95", `2000` = "1995-00", `2005` = "2000-05", `2010` = "2005-10", `2015` = "2010-15"),
         size.num = recode_factor(size, !!!size.level_key),
         size = factor(size, levels = c("tiny", "small", "medium", "large"))) %>%
  left_join(demo.sp_size %>% mutate(mean.mrate = mrate, mean.grate = grate) %>%
              select(sp_size, mean.mrate, mean.grate), by = "sp_size") %>%
  mutate(diff.mrate = mrate - mean.mrate)
grate.long <- grate.long %>%
  separate(sp_size, c("sp", "size", sep = "_"), remove = FALSE, extra = "drop", fill = "right") %>%
  select(-"_") %>%
  left_join(deci, by = "sp") %>%
  mutate(censusint.m = recode(interval, `1` = "1982-85", `2` = "1985-90",
                              `3` = "1990-95", `4` = "1995-00", `5` = "2000-05", `6` = "2005-10", `7` = "2010-15"),
         size.num = recode_factor(size, !!!size.level_key),
         size = factor(size, levels = c("tiny", "small", "medium", "large"))) %>%
  left_join(demo.sp_size %>% mutate(mean.mrate = mrate, mean.grate = grate) %>%
              select(sp_size, mean.mrate, mean.grate), by = "sp_size") %>%
  group_by(sp_size) %>%
  ungroup(sp_size)

save(grate.long, file = file.path("results/grate.long_by_species-size_decisuousness.Rdata"))
save(mrate.long, file = file.path("results/mrate.long_by_species-size_decisuousness.Rdata"))
#******************************************************
### Load Psi from ELM-FATES-------
#******************************************************

census.meds <- readr::read_rds("results/census.mediandates.rds")
census.beg <- census.meds[3: length(census.meds)]
cut.breaks <- census.beg
cut.breaks.2 <- as.Date(paste0(seq(1990, 2015, by = 5), "-01-01"))
cut.labels.2 <- paste0(seq(1990, 2010, by = 5), "-", seq(1995, 2015, by = 5))
cut.labels.interval <- 3: (length(census.meds)-1)
load(file.path("data-raw/psi.rda"))

## by depth panels

psi <- psi %>%
  mutate(interval.yrs = forcats::fct_explicit_na(cut(date, include.lowest = TRUE, breaks = cut.breaks,
                                                     labels = cut.labels.2, right = TRUE)))
#******************************************************
### Load climate data-------
#******************************************************

clim.dat <- read.csv("data-raw/BCI_1985_2018c_mod_2018substituted.csv")
str(clim.dat)
# converting to system time zone

clim.dat <- clim.dat %>%
  mutate(datetime = as.POSIXct(DateTime, format = "%m/%d/%y %H:%M",
                               tz = "America/Panama"),
         Year = format(datetime, "%Y")) %>%
  select(-DateTime) %>%
  mutate(date = as.Date(datetime, tz = "America/Panama"))
head(clim.dat$datetime)
# the above takes in the given data as in Panama format, but converts it to Sys timezone

clim.dat.yr <- clim.dat %>% subset(!is.na(Year) & Year != 2019) %>%
  group_by(Year) %>%
  summarise(rain = sum(Rainfall_mm_hr, na.rm = TRUE)) # m3/m2*1000 == L/m2 == mm
ggplot(clim.dat.yr, aes(x = Year, y = rain)) +
  geom_bar(stat = "identity")

clim.daily.rain <- clim.dat %>% select(-datetime) %>%
  select(date,  Rainfall_mm_hr) %>%
  group_by(date) %>%
  summarise(Precip = sum(Rainfall_mm_hr, na.rm = T))

clim.daily.min.max <- clim.dat %>%
  select(-datetime) %>%
  group_by(date) %>%
  summarise(Tmax = max(Temp_deg_C, na.rm = TRUE),
            Tmin = min(Temp_deg_C, na.rm = TRUE),
            RHmax = max(RH_prc, na.rm = TRUE),
            RHmin = min(RH_prc, na.rm = TRUE)) %>%
  mutate(RHmax = ifelse(RHmax > 100, 100, RHmax),
         RHmax = ifelse(RHmax < 0, 0, RHmin),
         RHmin = ifelse(RHmin > 100, 100, RHmax),
         RHmin = ifelse(RHmin < 0, 0, RHmin))

clim.daily <- clim.dat %>%
  select(-datetime, -Rainfall_mm_hr) %>%
  group_by(date, Year) %>%
  summarise(Rs = mean(SR_W_m2, na.rm = TRUE)* 0.0864, # from W/m2 MJ/m2/day
            uz = mean(Wind_Speed_km_hr, na.rm = TRUE)* 1000 / 3600, # to convert in m/s from km/hr
            Temp = mean(Temp_deg_C, na.rm = TRUE),
            RH = mean(RH_prc, na.rm = TRUE),
            Bp = mean(BP_Pa, na.rm = TRUE)) %>%
  ungroup(date, Year) %>%
  left_join(clim.daily.rain, by = "date") %>%
  left_join(clim.daily.min.max, by = "date") %>%
  mutate(Month = format(date, format = "%m"),
         Day = format(date, format = "%e"))

nrow(clim.daily[duplicated(clim.daily$date), ])
## there are no missing values

clim.daily.processed_data1 <- ReadInputs(varnames = c("Tmax", "Tmin", "RHmax", "RHmin", "uz", "Rs"),
                              subset(clim.daily, select = c(-date)),
                              constants.bci,
                              stopmissing = c(10, 10, 3),
                              timestep = "daily",
                              interp_missing_days = TRUE,
                              interp_missing_entries = TRUE,
                              interp_abnormal = TRUE,
                              missing_method = "DoY average",
                              abnormal_method = "DoY average"
)
#data("processeddata")
data("constants")
## Changing few constants to be BCI specific
constants.bci <- constants
constants.bci$z <- 2  # height of wind instrument
constants.bci$Elev <-
  150 # high 175 m, low 96m, so may be for now 150 m
constants.bci$lat <- 9.1666 # North, Long -79.85 W
constants.bci$lat_rad <- constants.bci$lat * pi / 180
constants.bci$PA <- 2650 # mm, average from 1929-2016, from this dataset
# Dont know these locally, so unchanged
# as fraction of extraterrestrial radiation reaching earth on sunless days = 0.23 for Australia (Roder- ick, 1999, page 181),
# bs difference between fracion of extraterrestrial radiation reaching full-sun days and that on sunless days = 0.5 for Australia (Roderick, 1999, page 181),

results_PM <-
  ET.PenmanMonteith(data = clim.daily.processed_data1,
                    constants = constants.bci, solar = "data", wind = "yes")

# Penman-Monteith FAO56 Reference Crop ET
# Evaporative surface: FAO-56 hypothetical short grass, albedo = 0.23 ; surface resistance = 70 sm^-1; crop height = 0.12  m; roughness height = 0.02 m
# Solar radiation data have been used directly for calculating evapotranspiration
# Wind data have been used for calculating the reference crop evapotranspiration
# Timestep: daily
# Units: mm
# Time duration: 1984-12-31 to 2019-01-01
# 12420 ET estimates obtained
# Basic stats
# Mean: 3.53
# Max: 7.16
# Min: 0.56
clim.daily <- clim.daily %>%
  left_join(data.frame(date = as.Date(index(results_PM$ET.Daily), tz = "America/Panama"),
                       pet.PM = coredata(results_PM$ET.Daily)), by = "date") %>%
  # Calculate the VPD = SVP x (1 – RH/100) = VPD
  mutate(SVP = 0.61121 * exp((18.678 - Temp/234.84) * (Temp/(257.14 + Temp))),
         VPD = SVP * (1 - RH/100),
         interval.yrs = forcats::fct_explicit_na(cut(date, include.lowest = TRUE, breaks = cut.breaks,
                                                     labels = cut.labels.2, right = TRUE)))

write.csv(clim.daily,
          file.path("results/clim.daily_with_pet.PM.csv"),
          row.names = FALSE)
save(clim.daily,
        file = file.path("results/clim.daily_with_pet.PM.Rdata"))
load(file = file.path("results/clim.daily_with_pet.PM.Rdata"))
#******************************************************
### Calculate Correlation of growth rates with psi by depth -------
#******************************************************
## 1. psi: interval mean psi by depth
## 2. psi.p50.g1 or g2: interval mean psi by depth: but only for the period when psi was above p50/p80
## 3. psi.p50.g1.n: number of days psi was above p50/p80 or group 1 p50 vs. group 2 p80

## Limiting to evergreen species for 2 & 3: For non-evergreen species periods in question will have to be a subset of days when leaves were on

growth.sub <- growth[grep("large", names(growth))]
mort.sub <- mrate.long %>% subset(size == "large") %>%
  mutate(interval = as.numeric(recode(census, `1985` = "1",
                              `1990` = "2", `1995` = "3",
                              `2000` = "4", `2005` = "5",
                              `2010` = "6", `2015` = "7"))) %>%
  subset(interval %in% cut.labels.interval) %>%
  select(sp_size, interval, mrate)

## combining psi, PET and VPD
soil.depths <- unique(psi$depth)

psi.m <- psi %>%
  mutate(interval = cut(date, include.lowest = TRUE, breaks = cut.breaks,
                        labels = cut.labels.interval, right = TRUE)) %>%
    full_join(clim.daily %>%
                mutate(std.pet.PM = range01(pet.PM),
                       std.VPD = range01(VPD),
                       std.Rs = Rs) %>%
                select(date, std.Rs, std.pet.PM, std.VPD), by = "date")
depth.breaks <- c(soil.depths[1], soil.depths[4], soil.depths[6],
                 soil.depths[7:length(soil.depths)])
depth.labels <- c(0.1, 0.4, signif(soil.depths[7:(length(soil.depths)-1)], 1),
                  signif(length(soil.depths), 2))

psi.study <- as.data.table(psi.m)[!is.na(interval),][,':='(doy = format(date, "%j"), year = format(date, "%Y"))][ doy < 120 & (!year %chin% c("1990")) &
                (!year %chin% c("1991") | depth == 2.9) &
                (!year %chin% c("1991", "1992") | depth < 2.9)][, c("doy", "year") := NULL][, ':='(depth = cut(depth, include.lowest = TRUE, breaks = depth.breaks,
              labels = depth.labels, right = TRUE))][, by = .(date, interval, interval.yrs, par.sam, depth),
                                                     lapply(.SD, mean, na.rm = TRUE)]
psi.study <- as.data.frame(psi.study)
tlp.sp <- read.csv("data-raw/traits/HydraulicTraits_Kunert/tlp_sp_mean.csv") # Nobby's data
n.tlp.sp <- length(tlp.sp$sp)
tlp.sp <- tlp.sp %>%
  mutate(psi_threshold = tlp*0.8)

tlp.sp.ls <- split(tlp.sp, f = list(tlp.sp$sp), drop = TRUE)

## for species with traits data
# growth.sub <- growth[names(growth) %in% paste0(unique(c(hyd$sp, traits$sp)), "_large")]
names.gfac <- names(psi.corr.fun.ls)#c("g.Psi", "g.Psi.Rad.VPD", "g.Psi.Rad.PET")
## Preparing PSI matrices to compare against
psi.interval <- vector(mode = "list", length = length(names.gfac))
names(psi.interval) <- names.gfac  # "psi.p50.g1", "psi.p50.g2"

for (i in 1: length(names.gfac)) {
    psi.interval[[i]] <- lapply(lapply(tlp.sp.ls, psi.corr.fun.ls[[i]]),
                             as.data.frame) %>%
      bind_rows(.id = "sp")
}

ml.ls <- vector("list", length(psi.interval))
ml.dens <- vector("list", length(psi.interval))
ml.rsq.ls <- vector("list", length(psi.interval))
for (i in names(psi.interval)) {
  if(grepl("gr.", i)) {
    demo <- lapply(growth.sub, as.data.frame) %>% bind_rows(.id = "sp_size") %>%
      rename(demo.rate = median)
  } else {
    demo <- mort.sub %>%
      rename(demo.rate = mrate)
  }
  demo.psi <- demo %>%
    mutate(interval = as.factor(interval)) %>%
    left_join(psi.interval[[i]] %>%
                mutate(sp_size = paste0(sp, "_large")), by = c("interval", "sp_size")) %>%
    subset(!is.na(`0.1`) & is.finite(`0.1`) & !is.na(demo.rate) & is.finite(demo.rate)) %>% droplevels()
  demo.psi.ls <- split(demo.psi,
                         f = list(demo.psi$sp_size), drop = TRUE)
  ml.ls[[i]] <- lapply(demo.psi.ls, get.ts.lk)
  ml.dens[[i]] <- sapply(ml.ls[[i]], get.ml.max)
  ml.rsq.ls[[i]] <- do.call(rbind, lapply(ml.ls[[i]], get.ml.depth.rsq))
}

# for(n in 4:length(ml.dens)) {
#   if(n == 1) plot(density(ml.dens[[n]][1:length(ml.dens[[n]])]))
#   lines(density(ml.dens[[n]]), col = terrain.colors(length(ml.dens))[n])
#   abline(v = median(ml.dens[[n]]))
# }

ml.rsq <- vector("list", length(psi.interval))
ml.rsq.best <- vector("list", length(psi.interval))
for (i in names(psi.interval)) {
  ml.rsq[[i]] <- ml.rsq.ls[[i]] %>%
    mutate(sp_size.depth = row.names(ml.rsq.ls[[i]])) %>%
    separate(sp_size.depth, c("sp", "size.depth", sep = "_"), remove = FALSE, extra = "drop", fill = "right") %>%
    select(-size.depth, -"_") %>%
    arrange(sp, depth) %>%
    left_join(deci %>% select(-sp4), by = "sp") %>%
    transform(deciduousness = factor(deciduousness,
                                     levels = c("Evergreen", "Brevideciduous",
                                                "Facultative Deciduous", "Obligate Deciduous"), ordered = TRUE)) %>%
    unite("deci_sp", deciduous, sp, remove = FALSE) %>%
    mutate(sp.plot = factor(sp, levels=unique(sp[order(deciduousness)]), ordered=TRUE),
           deci_sp.plot = factor(deci_sp, levels=unique(deci_sp[order(deciduousness)]), ordered=TRUE))

  # left_join(traits.long.hyd %>% select(deci_sp, deci_sp.plot, sp, sp.plot,  deciduousness), by = "sp") #%>%
  ml.rsq.best[[i]] <- ml.rsq[[i]] %>% group_by(sp) %>%
    mutate(R2.max = max(R2, na.rm = TRUE)) %>%
    subset(R2 == R2.max) %>% select(-R2.max) %>%
    ungroup(sp)

}
## species by depth by heat rsq
ml.rsq.combine <- dplyr::bind_rows(ml.rsq, .id = "variable") %>%
  transform(variable = factor(variable, levels = names.gfac))
ml.rsq.combine.best <- dplyr::bind_rows(ml.rsq.best, .id = "variable") %>%
  transform(variable = factor(variable, levels = names.gfac))

load(file = "data-raw/traits/isotopes/Oecologia 1995 Jackson_Fig3_Fig4_& Meinzer 1999_Fig4.Rdata")

ml.rsq.combine.best <- ml.rsq.combine.best %>%
  left_join(iso.1.3.join %>% select(sp, Xylem_sap_deltaD_permil, se, source), by = "sp")

ml.rsq.combine.best <- ml.rsq.combine.best %>%
  left_join(iso.1.3.join %>%
              group_by(sp) %>%
              summarise(Xylem_sap_deltaD_permil.mean = mean(Xylem_sap_deltaD_permil, na.rm = TRUE),
                        se.mean = mean(se, na.rm = TRUE)) %>%
              select(sp, Xylem_sap_deltaD_permil.mean, se.mean), by = "sp")

depth.rsq.isotopes <- ml.rsq.combine.best %>%
  group_by(variable, sp) %>%
  summarise_at(c("depth", "Xylem_sap_deltaD_permil", "se"), mean, na.rm = TRUE)

save(depth.rsq.isotopes, file = file.path(results.folder, "depth.rsq.isotopes.Rdata"))
save(ml.rsq.combine.best, file = file.path(results.folder, "ml.rsq.combine.best.Rdata"))
save(ml.rsq.combine, file = file.path(results.folder, "ml.rsq.combine.Rdata"))

#******************************************************
## Load BCI traits---
#******************************************************

bci.traits <- read.csv("data-raw/traits/BCITRAITS_20101220.csv") %>%
rename(form1 = GRWFRM1., sp = SP., WSG_Chave = WSG_CHAVE) %>% mutate(sp = tolower(sp))

#******************************************************
### Load Hydraulic traits by Brett Wolfe ---------
#******************************************************

hyd <- read.csv("data-raw/traits/HydraulicTraits_BrettWolfe/ht1_20200103.csv") #  # Brett's data
hyd <- hyd %>% select(-genus, -species, -deciduousness, -site) %>%
  rename(LMA = lma_gm2_m, WD = xylem_den_m, TLP = tlp_m,
         p50S = p50, p80S = p80, p88S = p88,
         Fcap_Xylem = cwr_xylem_elbow_origin_slope, CWR_Xylem = cwr_xylem_cwr_at_elbow, Felbow_Xylem = cwr_xylem_elbow,
         Fcap_Bark = cwr_bark_elbow_origin_slope, CWR_Bark = cwr_bark_cwr_at_elbow, Felbow_Bark = cwr_bark_elbow,
         BarkThick = barkthickness10mm) %>%
  ## Given Fcap slope values are -ve but they should be positive
  mutate(sp = tolower(sp), Fcap_Xylem = -Fcap_Xylem, Fcap_Bark = -Fcap_Bark) %>%
  left_join(lwp.min.wide %>% subset(location = "PA-BCI"), by = "sp") %>%
  mutate(HSMTLP.50S =  TLP - p50S, HSM50S = lwp.min_Diurnal - p50S,
         HSMTLP.88S =  TLP - p88S, HSM88S = lwp.min_Diurnal - p88S,
         HSMTLP = lwp.min_Diurnal - TLP,
         HSMFelbow_Xylem = lwp.min_Diurnal - Felbow_Xylem,
         HSMFelbow_Bark = lwp.min_Diurnal - Felbow_Bark,
         CWR_Total = CWR_Xylem + CWR_Bark) %>%
  left_join(deci %>% select(-sp4, -deciduousness.label), by = "sp") %>%
  left_join(sp.hab, by = "sp") %>%
  left_join(depth.rsq.isotopes, by = "sp")
length(unique(hyd$sp)) # 27 sp across BCI, PNM, San Lorenzo
traits.labels.table.1 <- data.frame(trait = factor(c("depth", "Xylem_sap_deltaD_permil",
                                              "lwp.min_Predawn", "lwp.min_Diurnal", "TLP",
                                              "p50S", "p88S",
                                              "HSMTLP", "HSM50S","HSM88S",
                                              "HSMTLP.50S", "HSMTLP.88S",
                                              "CWR_Total", "CWR_Xylem", "CWR_Bark",
                                              "Felbow_Xylem", "Felbow_Bark", "HSMFelbow_Xylem", "HSMFelbow_Bark",
                                              "Fcap_Xylem", "Fcap_Bark","WD",
                                              "Panama.moist.pref", "Plot.swp.pref", "LMA"),
                                    levels = c("depth", "Xylem_sap_deltaD_permil",
                                               "lwp.min_Predawn", "lwp.min_Diurnal", "TLP",
                                               "p50S", "p88S",
                                               "HSMTLP", "HSM50S","HSM88S",
                                               "HSMTLP.50S", "HSMTLP.88S",
                                               "CWR_Total", "CWR_Xylem", "CWR_Bark",
                                               "Felbow_Xylem", "Felbow_Bark", "HSMFelbow_Xylem", "HSMFelbow_Bark",
                                               "Fcap_Xylem", "Fcap_Bark","WD",
                                               "Panama.moist.pref", "Plot.swp.pref", "LMA"), ordered = TRUE)) %>%
  transform(trait.plot = factor(trait, labels = c(expression(Depth[italic('Rsq')]), expression(italic(delta)^2*H[Xylem]),
                                            expression(Psi[Predawn]), expression(Psi[min]), expression(Psi[TLP]),
                                            expression(italic('P')['50, Stem']),  expression(italic('P')['88, Stem']),
                                            expression(Psi[min]*' - '*Psi[TLP]),
                                            expression(Psi[min]*' - '*italic('P')['50, Stem']),
                                            expression(Psi[min]*' - '*italic('P')['88, Stem']),
                                            expression(Psi[TLP]*' - '*italic('P')['50, Stem']),
                                            expression(Psi[TLP]*' - '*italic('P')['88, Stem']),
                                            expression('CWR'['Total']), expression('CWR'['Xylem']), expression('CWR'['Bark']),
                                            expression(italic('F')['Elbow, Xylem']), expression(italic('F')['Elbow, Bark']),
                                            expression(Psi[min]*'-'*italic('F')['Elbow, Xylem']),
                                            expression(Psi[min]*'-'*italic('F')['Elbow, Bark']),
                                            expression(italic('F')['Cap, Xylem']), expression(italic('F')['Cap, Bark']),
                                            expression('WD'[stem]),
                                            expression('Panama'[wet]), expression('Plot'[wet]), "LMA")),
             trait.plot.chart = factor(trait, labels = c(expression(Depth[italic('Rsq')]), expression(italic(delta)^2*H[Xylem]),
                                                  expression(Psi[Predawn]), expression(Psi[min]), expression(Psi[TLP]),
                                                  expression(italic('P')['50,Stem']),  expression(italic('P')['88,Stem']),
                                                  expression(Psi[min]*'-'*Psi[TLP]),
                                                  expression(Psi[min]*'-'*italic('P')['50,Stem']),
                                                  expression(Psi[min]*'-'*italic('P')['88,Stem']),
                                                  expression(Psi[TLP]*'-'*italic('P')['50,Stem']),
                                                  expression(Psi[TLP]*'-'*italic('P')['88,Stem']),
                                                  expression('CWR'['Total']), expression('CWR'['Xylem']), expression('CWR'['Bark']),
                                                  expression(italic('F')['Elbow, Xylem']), expression(italic('F')['Elbow,Bark']),
                                                  expression(Psi[min]*'-'*italic('F')['Elbow,Xylem']),
                                                  expression(Psi[min]*'-'*italic('F')['Elbow,Bark']),
                                                  expression(italic('F')['Cap,Xylem']), expression(italic('F')['Cap,Bark']),
                                                  expression('WD'[stem]),
                                                  expression('Panama'[wet]), expression('Plot'[wet]), "LMA"))) %>% droplevels()



#******************************************************
####----Phenology by Brett Wolfe hydraulic traits-----
#******************************************************
hyd.long <- hyd %>% select(-DeciLvl) %>%
  select(sp, deciduousness, deciduous, location, TLP, p50S, p88S,
         CWR_Total, Fcap_Xylem, CWR_Xylem, Felbow_Xylem, Fcap_Bark, CWR_Bark, Felbow_Bark, WD, LMA,
         HSM50S, HSM88S, HSMTLP, HSMFelbow_Xylem, HSMFelbow_Bark, HSMTLP.50S, HSMTLP.88S,
         Panama.moist.pref, Plot.swp.pref, lwp.min_Diurnal, lwp.min_Predawn, depth, Xylem_sap_deltaD_permil, se) %>%
  gather(trait, value, -sp, -deciduousness, -deciduous,  -location, -se) %>%
  unite("deci_sp", deciduous, sp, remove = FALSE) %>%
  subset(deciduousness != "NA") %>%
  droplevels()
kruskal.list <- list()
for(i in unique(hyd.long$trait)) {
  xx <- hyd.long %>% subset(trait == i)
  kruskal.list[[i]] <- cbind(trait = i, kruskal(xx$value, xx$deciduousness, alpha = 0.1, group=TRUE, p.adj="bonferroni")$groups,
                             deciduousness = rownames(kruskal(xx$value, xx$deciduousness, alpha = 0.1, group=TRUE, p.adj="bonferroni")$groups))
}
hyd.kruskal.labels <- do.call(rbind.data.frame, kruskal.list)
head(hyd.kruskal.labels)
hyd.labels <- hyd.kruskal.labels
unique(hyd.labels$trait)
hyd.labels.data <- hyd.labels %>%
  left_join(hyd.long %>% group_by(trait) %>%
              summarise(value = max(value, na.rm = TRUE)), by = c("trait")) %>%
  subset(deciduousness != "NA") %>%
  droplevels() %>%
  left_join(traits.labels.table.1 %>% select(trait, trait.plot), by = "trait")

hyd.long <- hyd.long %>%
  left_join(traits.labels.table.1, by = "trait")

save(hyd, file = file.path(results.folder, "hyd.traits.all.RData"))
save(hyd.long, file = file.path(results.folder, "hyd.traits.key.long.RData"))

#******************************************************
## Load Nobert Kunert traits --------
#******************************************************

traits.indi <- read.csv("data-raw/traits/HydraulicTraits_Kunert/hydraulic_traits_panama_kunert.csv") # Nobby's data
traits <- traits.indi %>%
  rename(TLP = mean_TLP_Mpa, Chl = Chl_m2_per_g) %>%
  select(sp, TLP, Chl) %>% group_by(sp) %>%
  summarise_all(mean, na.rm = TRUE)
traits <- traits %>%
  left_join(lwp.min.wide %>% select(sp, lwp.min_Predawn, lwp.min_Diurnal) %>%
              subset(location = "PA-BCI"), by = "sp")

leaf_cond.models <- read.csv("data-raw/traits/HydraulicTraits_Kunert/Panama_fits_leaf_K_p50_Kunert.csv")
leaf.k.p80 <- leaf_cond.models %>% subset(model_type == "Exponential") %>%
  mutate(sp = data.type, p50L = -psi_kl50, p80L = -psi_kl80, KmaxL = Kmax) %>% # these are Kmax that are extrapolated from the exponential curve
  select(sp, KmaxL, p50L, p80L) # 21 species
traits <- traits %>%
  full_join(leaf.k.p80 %>% mutate(sp = as.character(sp)), by = "sp") %>%
  mutate(HSMTLP.50L = TLP - p50L, HSMLWP.50L = lwp.min_Diurnal - p50L,
         HSMTLP.80L = TLP - p80L, HSMLWP.80L = lwp.min_Diurnal - p80L,
         HSMLWP.TLP = lwp.min_Diurnal - TLP) %>%
  left_join(sp.hab, by = "sp") %>%
  left_join(deci %>% select(-sp4, -deciduousness.label), by = "sp") %>%
  left_join(bci.traits %>% select(sp, form1, WSG_Chave), by = "sp") %>%
  left_join(depth.rsq.isotopes, by = "sp")

traits.labels.table.2 <- data.frame(trait = factor(c("depth", "Xylem_sap_deltaD_permil",
                                              "KmaxL", "lwp.min_Predawn", "lwp.min_Diurnal", "TLP", "p50L", "p80L",
                                              "HSMLWP.TLP", "HSMLWP.50L", "HSMTLP.50L",
                                              "HSMLWP.80L", "HSMTLP.80L",
                                              "Panama.moist.pref", "Plot.swp.pref", "WSG_Chave", "Chl"),
                                              levels = c("depth", "Xylem_sap_deltaD_permil",
                                                           "KmaxL", "lwp.min_Predawn", "lwp.min_Diurnal", "TLP", "p50L", "p80L",
                                                           "HSMLWP.TLP", "HSMLWP.50L", "HSMTLP.50L",
                                                           "HSMLWP.80L", "HSMTLP.80L",
                                                           "Panama.moist.pref", "Plot.swp.pref", "WSG_Chave", "Chl"), ordered = TRUE)) %>%
  transform(trait.plot = factor(trait, labels = c(expression(Depth[italic('Rsq')]), expression(italic(delta)^2*H[Xylem]),
                                           expression(italic(K)[max]), expression(Psi[predawn]), expression(Psi[min]),
                                           expression(Psi[TLP]), expression(italic('P')['50, Leaf']), expression(italic('P')['80, Leaf']),
                                           expression(Psi[min]*' - '*Psi[TLP]),
                                           expression(Psi[min]*' - '*italic('P')['50, Leaf']),
                                           expression(Psi[TLP]*' - '*italic('P')['50, Leaf']),
                                           expression(Psi[min]*' - '*italic('P')['80, Leaf']),
                                           expression(Psi[TLP]*' - '*italic('P')['80, Leaf']),
                                           expression('Panama'[wet]), expression('Plot'[wet]), expression('WSG'[Chave]), "LMA")),
            trait.plot.chart = factor(trait, labels = c(expression(Depth[italic('Rsq')]), expression(italic(delta)^2*H[Xylem]),
                                                 expression(italic(K)[max]), expression(Psi[predawn]), expression(Psi[min]),
                                                 expression(Psi[TLP]), expression(italic('P')['50,Leaf']), expression(italic('P')['80,Leaf']),
                                                 expression(Psi[min]*'-'*Psi[TLP]),
                                                 expression(Psi[min]*'-'*italic('P')['50,Leaf']),
                                                 expression(Psi[TLP]*'-'*italic('P')['50,Leaf']),
                                                 expression(Psi[min]*'-'*italic('P')['80,Leaf']),
                                                 expression(Psi[TLP]*'-'*italic('P')['80,Leaf']),
                                                 expression('Panama'[wet]), expression('Plot'[wet]), expression('WSG'[Chave]), "LMA")))

# > with(traits, table(deciduousness))
# Evergreen    Obligate Deciduous Facultative Deciduous        Brevideciduous
# 26                     3                    11                     8

#******************************************************
####----Phenology by Kunert hydraulic traits-----
#******************************************************

traits.long <- traits %>% select(-DeciLvl) %>%
  gather(trait, value, -sp, -deciduousness, -deciduous, -form1) %>%
  subset(deciduousness != "NA") %>%
  droplevels()

kruskal.list <- list()
for(i in unique(traits.long$trait)) {
  xx <- traits.long %>% subset(trait == i)
  kruskal.list[[i]] <- cbind(trait = i, kruskal(xx$value, xx$deciduousness, alpha = 0.1, group=TRUE, p.adj="bonferroni")$groups,
                             deciduousness = rownames(kruskal(xx$value, xx$deciduousness, alpha = 0.1, group=TRUE, p.adj="bonferroni")$groups))
}
traits.kruskal.labels <- do.call(rbind.data.frame, kruskal.list)
head(traits.kruskal.labels)

traits.labels <- traits.kruskal.labels
traits.labels.data <- traits.labels %>%
  left_join(traits.long %>% group_by(trait) %>%
              summarise(value = max(value, na.rm = TRUE)), by = c("trait")) %>%
  subset(deciduousness != "NA") %>%
  droplevels() %>%
  transform(deciduousness = factor(deciduousness,
                                   levels = c("Evergreen", "Brevideciduous",
                                              "Facultative Deciduous", "Obligate Deciduous"), ordered = TRUE)) %>%
  left_join(traits.labels.table.2 %>% select(trait, trait.plot), by = "trait")

traits.long <- traits.long %>%
  left_join(traits.labels.table.2, by = "trait")

## Kunert traits species wise for sp in hyd.traits----
select.traits <- c("lwp.min_Predawn", "lwp.min_Diurnal", "TLP", "p50S", "p88S",
                   "HSMTLP", "HSM50S", "HSM88S", "HSMTLP.50S", "HSMTLP.88S")

traits.long <- traits.long %>%
  unite("deci_sp", deciduous, sp, remove = FALSE) %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(deciduousness)]), ordered=TRUE),
         deci_sp.plot = factor(deci_sp, levels=unique(deci_sp[order(deciduousness)]), ordered=TRUE))
# just for sp with hyd.traits, but traits.long does not have all those sp, and hab preference and WSG traits will be missed
## so beginning with those other traits
traits.wide <- traits.long %>% select(-trait.plot, -trait.plot.chart) %>% pivot_wider(names_from = trait, values_from = value)
traits.long.hyd <- sp.hab %>%
  full_join(bci.traits %>% select(sp, form1, WSG_Chave), by = "sp") %>%
  full_join(deci %>% select(-sp4), by = "sp") %>%
  subset(sp %in% unique(c(depth.rsq.isotopes$sp, hyd$sp, traits$sp))) %>%
  left_join(traits.wide %>% select(-form1, -deciduous, -deciduousness,
                                   -WSG_Chave, -Panama.moist.pref, -Plot.swp.pref,
                                   -se), by = "sp") %>%
  pivot_longer(cols = c(-sp, -form1, -deciduous, -deciduousness, -DeciLvl),
               names_to = "trait", values_to = "value") %>%
  unite("deci_sp", deciduous, sp, remove = FALSE) %>%
  left_join(traits.labels.table.2, by = "trait") %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(deciduousness)]), ordered=TRUE),
         deci_sp.plot = factor(deci_sp, levels=unique(deci_sp[order(deciduousness)]), ordered=TRUE)) %>%
  droplevels()

save(traits, file = file.path(results.folder, "kunert.traits.all.RData"))
save(traits.long, file = file.path(results.folder, "kunert.traits.key.long.RData"))
save(traits.long.hyd, file = file.path(results.folder, "kunert.traits.key.long_in_Wolfe_traits_species_list.RData"))
load(file = file.path(results.folder, "kunert.traits.key.long.RData"))

### Plot best correlated depth against isotopic data and traits-----

# ml.rsq.combine.best <- ml.rsq.combine.best %>%
#   left_join(traits.wide.hyd %>% select(sp, KmaxL, Panama.moist.pref:HSMLWP.TLP), by = "sp")

heat.rsq <- ggplot(ml.rsq.combine,
                   aes(y = deci_sp, x = as.factor(depth))) +
  geom_tile(aes(fill = R2)) +
  ylab("Species") + xlab("Depth (m)") +
  facet_wrap(. ~ variable, nrow = 1) +
  scale_fill_viridis_c("R-squared", direction = -1, option = "plasma") #+
ggsave("psi.corr_all.depths_phenology_heat_by_variable.jpeg",
       plot = heat.rsq, file.path(figures.folder), device = "jpeg", height = 5, width = 12, units='in')

# theme(axis.text.y = element_text(angle = 90, vjust = 0.5)) +
# scale_x_continuous(breaks = soil.depths[c(1,8:13)])
heat.best.rsq <- heat.rsq %+% subset(ml.rsq.combine.best, R2 >= 0.3)
ggsave("psi.corr_best.depth_phenology_heat_by_variable.jpeg",
       plot = heat.best.rsq, file.path(figures.folder), device = "jpeg", height = 5, width = 12, units='in')


xylem.label <- expression('Xylem Sap '*delta~""^2*"H (\u2030)"*'')
ml.rsq.combine.best <- ml.rsq.combine.best %>% left_join(bci.traits %>% select(sp, form1), by = "sp")

formula = y~x
p0 <- ggplot(ml.rsq.combine.best %>% subset(R2 >= 0.3 & form1 == "T"),
             aes(x = Xylem_sap_deltaD_permil.mean, y = depth)) + #HSMTLP.80L)) +
  geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil.mean + se.mean,
                     xmin = Xylem_sap_deltaD_permil.mean - se.mean, color = deciduousness),
                 size = 0.2) +
  facet_wrap( ~ variable, nrow = 2) +
  geom_text(aes(x =  Xylem_sap_deltaD_permil.mean, y = depth, label = sp, color = deciduousness), nudge_y = 0.1, nudge_x = 0.2,
            size = 3, show.legend = FALSE) +
  ylab(expression("Water Uptake Depth (m)")) + xlab(xylem.label) +
  scale_y_continuous(trans=reverselog_trans(10), breaks = ml.rsq.combine.best$depth) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.6, npcy = 0.2, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.6, npcy = 0.1, size = 4) +
  geom_point(size = 1, show.legend = TRUE, aes(color = deciduousness)) + theme(legend.position = "top")
ggsave("psi.corr_best.depth_xylem_sap_deltaD_phenology_mean_isotope_source.jpeg",
       plot = p0, file.path(figures.folder), device = "jpeg", height = 6, width = 9, units = 'in')

p0 <- ggplot(ml.rsq.combine.best %>% subset(R2 >= 0.3 & form1 == "T"),
             aes(x = Xylem_sap_deltaD_permil, y = depth)) + #HSMTLP.80L)) +
  geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + se,
                     xmin = Xylem_sap_deltaD_permil - se, color = deciduousness),
                 size = 0.2) +
  facet_wrap( ~ variable, nrow = 2) +
  geom_text(aes(x =  Xylem_sap_deltaD_permil, y = depth, label = sp, color = deciduousness), nudge_y = 0.1, nudge_x = 0.2,
            size = 3, show.legend = FALSE) +
  ylab(expression("Water Uptake Depth (m)")) + xlab(xylem.label) +
  scale_y_continuous(trans=reverselog_trans(10), breaks = ml.rsq.combine.best$depth) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.6, npcy = 0.2, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.6, npcy = 0.1, size = 4) +
  geom_point(size = 1, show.legend = TRUE, aes(color = deciduousness)) + theme(legend.position = "top")
ggsave("psi.corr_best.depth_xylem_sap_deltaD_phenology_two_isotope_sources.jpeg",
       plot = p0, file.path(figures.folder), device = "jpeg", height = 6, width = 9, units = 'in')

g3 <- ggplot(ml.rsq.combine.best %>% subset(R2 >= 0.3 & !duplicated(sp) & !is.na(depth)),
             aes(x = deci_sp.plot, y = depth)) +
  facet_wrap(. ~ variable) +
  geom_col(aes(fill = deciduousness)) +
  guides(fill = guide_legend(title = "")) +
  theme(legend.position = "top") +
  ylab("Depth (m)") + xlab("Species") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(trans=reverselog_trans(10), breaks = ml.rsq.combine.best$depth)
g4 <- g3 + coord_flip() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))
ggsave("psi.corr_best.depth_phenology.jpeg",
       plot = g4, file.path(figures.folder), device = "jpeg", height = 5, width = 6.5, units='in')

### Plot against hydraulic traits

d1 <- ggplot(traits.long.hyd %>% subset(!is.na(deciduousness))) +
  facet_wrap(. ~  trait.plot, scales = "free_x", labeller = label_parsed, nrow = 2) +
  geom_col(aes(x = deci_sp.plot, y = value,
               fill = deciduousness),
           position = position_dodge2(width = 0.9, preserve = "single")) +
  theme(axis.text.y = element_text(face = "plain", size = 8)) +
  coord_flip() +
  guides(fill = guide_legend(title = "Deciduousness")) +
  xlab("Deciduousness") + ylab("Value")
ggsave(file.path(figures.folder, paste0("Kunert_traits_vs_deciduousness_sp_bar_depth_isotopes_all.sp.with.depth.jpeg")),
       plot = d1, height = 7, width = 14, units ='in')
d2 <- d1 %+% subset(traits.long.hyd, !is.na(deciduousness) & sp %in% traits$sp)
ggsave(file.path(figures.folder, paste0("Kunert_traits_vs_deciduousness_sp_bar_depth_isotopes.jpeg")),
       plot = d2, height = 7, width = 14, units ='in')

## Correlation chart

# Check correlations (as scatterplots), distribution and print correlation coefficient
select.traits.1 <- c("TLP", "p50S", "p88S", "HSMTLP", "HSM50S",
                     "lwp.min_Predawn", "lwp.min_Diurnal",
                     "HSM88S", "HSMTLP.50S", "HSMTLP.88S")
select.traits.2 <- c("depth", "Xylem_sap_deltaD_permil", "lwp.min_Predawn", "lwp.min_Diurnal",
                     "TLP",  "p88S", "HSMTLP", "HSM88S", "HSMTLP.88S")

hyd.pairs.1 <- hyd.long %>%
  subset(trait %in% select.traits.1) %>%
  select(sp, deciduousness, trait.plot.chart, value) %>%
  pivot_wider(names_from = trait.plot.chart, values_from = value)

hyd.pairs.2 <- hyd.long %>%
  subset(trait %in% select.traits.2) %>%
  select(sp, deciduousness, trait.plot.chart, value) %>%
  pivot_wider(names_from = trait.plot.chart, values_from = value)

select.traits.3 <- c("KmaxL", "lwp.min_Predawn",
                     "lwp.min_Diurnal", "TLP", "p50L",
  "HSMLWP.TLP", "HSMLWP.50L", "HSMTLP.50L")

select.traits.4 <- c("depth", "Xylem_sap_deltaD_permil",
                     "KmaxL", "lwp.min_Predawn", "TLP", "p50L",
                     "HSMLWP.TLP", "Panama.moist.pref", "Plot.swp.pref", "WSG_Chave", "Chl")

traits.pairs.1 <- traits.long %>%
  subset(trait %in% select.traits.3) %>%
  select(sp, deciduousness, trait.plot.chart, value) %>%
  pivot_wider(names_from = trait.plot.chart, values_from = value)
traits.pairs.2 <- traits.long %>%
  subset(trait %in% select.traits.4) %>%
  select(sp, deciduousness, trait.plot.chart, value) %>%
  pivot_wider(names_from = trait.plot.chart, values_from = value)

chart.hyd.1 <- ggpairs(hyd.pairs.1 %>% select(-sp, -deciduousness),
              upper = list(continuous = wrap(cor_func,
                                             method = 'spearman', symbol = expression('\u03C1 ='))),
              lower = list(continuous = function(data, mapping, ...) {
                ggally_smooth_lm(data = data, mapping = mapping) +
                  theme(panel.background = element_blank())}),
              diag = list(continuous = function(data, mapping, ...) {
                ggally_densityDiag(data = data, mapping = mapping) +
                  theme(panel.background = element_blank())}
              ), labeller = "label_parsed")
chart.hyd.2 <- ggpairs(hyd.pairs.2 %>% select(-sp, -deciduousness),
                       upper = list(continuous = wrap(cor_func,
                                                      method = 'spearman', symbol = expression('\u03C1 ='))),
                       lower = list(continuous = function(data, mapping, ...) {
                         ggally_smooth_lm(data = data, mapping = mapping) +
                           theme(panel.background = element_blank())}),
                       diag = list(continuous = function(data, mapping, ...) {
                         ggally_densityDiag(data = data, mapping = mapping) +
                           theme(panel.background = element_blank())}
                       ), labeller = "label_parsed")
# chart.hyd.2 <- chart.hyd.1 %+% hyd.pairs.2
# chart.traits.1 <- chart.hyd.1 %+% traits.pairs.1
# chart.traits.2 <- chart.hyd.1 %+% traits.pairs.2
chart.traits.1 <- ggpairs(traits.pairs.1 %>% select(-sp, -deciduousness),
                          upper = list(continuous = wrap(cor_func,
                                                         method = 'spearman', symbol = expression('\u03C1 ='))),
                          lower = list(continuous = function(data, mapping, ...) {
                            ggally_smooth_lm(data = data, mapping = mapping) +
                              theme(panel.background = element_blank())}),
                          diag = list(continuous = function(data, mapping, ...) {
                            ggally_densityDiag(data = data, mapping = mapping) +
                              theme(panel.background = element_blank())}
                          ), labeller = "label_parsed")
chart.traits.2 <- ggpairs(traits.pairs.2 %>% select(-sp, -deciduousness),
        upper = list(continuous = wrap(cor_func,
                                       method = 'spearman', symbol = expression('\u03C1 ='))),
        lower = list(continuous = function(data, mapping, ...) {
          ggally_smooth_lm(data = data, mapping = mapping) +
            theme(panel.background = element_blank())}),
        diag = list(continuous = function(data, mapping, ...) {
          ggally_densityDiag(data = data, mapping = mapping) +
            theme(panel.background = element_blank())}
        ), labeller = "label_parsed")

ggsave(file.path(figures.folder, paste0("Brett_Wolfe_traits_cor.chart.jpeg")),
       plot = chart.hyd.1 + ggpairs.theme, height = 10, width = 10, units ='in')
ggsave(file.path(figures.folder, paste0("Brett_Wolfe_traits_depth_isotopes_cor.chart.jpeg")),
       plot = chart.hyd.2 + ggpairs.theme, height = 10, width = 10, units ='in')
ggsave(file.path(figures.folder, paste0("Kunert_traits_cor.chart.jpeg")),
       plot = chart.traits.1 + ggpairs.theme, height = 8, width = 8, units ='in')
ggsave(file.path(figures.folder, paste0("Kunert_traits_depth_isotopes_cor.chart.jpeg")),
       plot = chart.traits.2 + ggpairs.theme, height = 9, width = 10, units ='in')

#******************************************************
### Plot LWP -----
#******************************************************

ggplot(lwp.all, aes(x = LWP_coll_time, y = lwp.min)) +
  facet_grid(location ~ .) +
  geom_point(aes(color = as.factor(date))) +
  geom_line(aes(group = c(sp_date), color = as.factor(date)), show.legend = FALSE) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE) +
  guides(color = guide_legend(title = "Date"))
ggsave(file.path(figures.folder, paste0("LWP_time_series_all_Data.jpeg")), height = 6, width = 5, units ='in')

# Across all the dryseason measurements which diurnal measurement was minimum by species and location
## also get the predawn for the same day
# https://stackoverflow.com/questions/46971945/how-can-i-have-a-greek-symbol-within-one-of-my-facet-labels
plot.lwp.base <- ggplot(lwp.min %>% subset(!is.na(deciduousness)), aes(x = time, y = lwp.min)) +
  facet_grid(. ~ location) +
  geom_point(aes(color = deciduousness)) +
  ylab(expression(psi[min])) + xlab("Time") +
  guides(colour = guide_legend(order = 1, title = "Deciduousness")) +
  scale_x_discrete(name="",
                   breaks = c("Predawn", "Diurnal"),
                   labels = c("Predawn", expression('Diurnal'['min'])))
plot.lwp <- plot.lwp.base +
  geom_line(aes(group = sp, color = deciduousness)) +
  geom_errorbar(aes(ymax = lwp.min + lwp.se, ymin = lwp.min - lwp.se, color = deciduousness), width = 0.05) +
  ggsave(file.path(figures.folder, paste0("LWP_min.jpeg")), plot = plot.lwp, height = 4, width = 7.5, units ='in')

plot.lwp.diff <- plot.lwp.base %+%
  subset(lwp.diff, !is.na(deciduousness)) +
  geom_line(aes(group = sp_date, color = deciduousness), size = 0.2)
ggsave(file.path(figures.folder, paste0("LWP_min_all_days.jpeg")), plot = plot.lwp.diff, height = 4, width = 7.5, units ='in')

plot.lwp.diff.ts.base <- ggplot(lwp.diff %>% subset(!is.na(deciduousness) & !is.na(lwp.diff) & location != "PA-BCI"),
                                aes(x = date, y = lwp.diff)) +
  facet_grid(. ~ location) +
  ylab(expression(psi[Predawn] - psi['Diurnal'['min']])) + xlab("Time")
plot.lwp.diff.ts <- plot.lwp.diff.ts.base +
  geom_point(aes(color = deciduousness)) +
  geom_line(aes(group = sp, color = deciduousness), size = 0.2) +
  guides(colour = guide_legend(order = 1, title = "Deciduousness"))
ggsave(file.path(figures.folder, paste0("LWP_min_all_days.jpeg")), plot = plot.lwp.diff.ts, height = 4, width = 7.5, units ='in')
plot.lwp.diff.ts.sp <- plot.lwp.diff.ts.base +
  geom_point(aes(shape = deciduousness, color = deci_sp), size = 3) +
  geom_line(aes(group = deci_sp, color = deci_sp), size = 0.2) +
  guides(colour = guide_legend(order = 1, title = "Sp"),
         shape = guide_legend(order = 2, title = "Deciduousness")) +
  theme(legend.direction = "vertical", legend.box = "horizontal")
ggsave(file.path(figures.folder, paste0("LWP_min_all_days_sp_color.jpeg")), plot = plot.lwp.diff.ts.sp, height = 4, width = 8, units ='in')

plot.lwp.ts.sp.pnm <- ggplot(lwp.diff %>% subset(!is.na(deciduousness) & location == "PA-PNM") %>% droplevels(),
                             aes(x = time, y = lwp.min)) +
  facet_grid(location ~ as.factor(date)) +
  ylab(expression(psi[min])) + xlab("Time") +
  ylim(c(-3.5, 0)) +
  scale_x_discrete(name="",
                   breaks = c("Predawn", "Diurnal"),
                   labels = c("Predawn", expression('Diurnal'['min']))) +
  geom_point(aes(color = deci_sp, shape = deciduousness), size = 2) +
  geom_line(aes(group = deci_sp, color = deci_sp), size = 0.2) +
  guides(colour = guide_legend(order = 1, title = "Sp"),
         shape = guide_legend(order = 2, title = "Deciduousness")) +
  theme(legend.direction = "vertical", legend.box = "horizontal")
plot.lwp.ts.sp.snl <- plot.lwp.ts.sp.pnm %+%
  droplevels(subset(lwp.diff, !is.na(deciduousness) & location == "PA-SLZ")) +
  scale_color_brewer(palette = "Set1")
ggsave(file.path(figures.folder, paste0("LWP_min_all_days_diurnal_sp_color.jpeg")),
       plot = arrangeGrob(plot.lwp.ts.sp.pnm, plot.lwp.ts.sp.snl), height = 5.5, width = 12, units ='in')


#******************************************************
####----Plot Phenology by Wolfe hydraulic traits-----
#******************************************************

h1 <- ggplot(hyd.labels.data, aes(x = deciduousness, y = value)) +
  geom_boxplot(data = hyd.long, aes(fill = deciduousness), stat = "boxplot", notch = TRUE) +
  geom_jitter(data = hyd.long, width = 0.05, shape = 21, fill = "darkgray", color = "black", show.legend = FALSE, alpha = 0.7) +
  facet_wrap(. ~  trait.plot, scales = "free_y", labeller = label_parsed) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 1, hjust = 1)) +
  scale_color_brewer(palette = "Greens", direction = -1) +
  scale_fill_brewer(name = "Deciduousness", palette = "Greens", direction = -1) +
  xlab("Deciduousness") + ylab("Value")
ggsave(file.path(figures.folder, paste0("BrettWolfe_traits_vs_deciduousness_isotopes_depth.jpeg")),
       plot = h1, height = 7, width = 10, units ='in')
h1.1 <- h1 + geom_text(aes(label = groups), vjust = 1, hjust = 0, show.legend = FALSE)
ggsave(file.path(figures.folder, paste0("BrettWolfe_traits_vs_deciduousness_kruskal.labels.jpeg")),
       plot = h1.1, height = 7, width = 10, units ='in')

h1.2 <- ggplot(hyd.labels.data %>% subset(!trait %in% c("depth",  "Xylem_sap_deltaD_permil")),
               aes(x = deciduousness, y = value)) +
  geom_boxplot(data = hyd.long %>% subset(!trait %in% c("depth",  "Xylem_sap_deltaD_permil")),
               aes(fill = deciduousness), stat = "boxplot", notch = TRUE) +
  geom_jitter(data = hyd.long %>% subset(!trait %in% c("depth",  "Xylem_sap_deltaD_permil")),
              width = 0.05, shape = 21, fill = "darkgray", color = "black", show.legend = FALSE, alpha = 0.7) +
  facet_wrap(. ~  trait.plot, scales = "free_y", labeller = label_parsed) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 1, hjust = 1)) +
  scale_color_brewer(palette = "Greens", direction = -1) +
  scale_fill_brewer(name = "Deciduousness", palette = "Greens", direction = -1) +
  xlab("Deciduousness") + ylab("Value")
ggsave(file.path(figures.folder, paste0("BrettWolfe_traits_vs_deciduousness.jpeg")),
       plot = h1.2, height = 7, width = 10, units ='in')
h1.3 <- h1.2 + geom_text(aes(label = groups), vjust = 1, hjust = 0, show.legend = FALSE)
ggsave(file.path(figures.folder, paste0("BrettWolfe_traits_vs_deciduousness_kruskal.labels.jpeg")),
       plot = h1.3, height = 7, width = 10, units ='in')

select.traits <- c("depth", "Xylem_sap_deltaD_permil", "lwp.min_Predawn", "lwp.min_Diurnal", "TLP", "p50S", "p88S",
                   "HSMTLP", "HSM50S","HSM88S", "HSMTLP.50S", "HSMTLP.88S")
hyd.long <- hyd.long %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(deciduousness)]), ordered=TRUE),
         deci_sp.plot = factor(deci_sp, levels=unique(deci_sp[order(deciduousness)]), ordered=TRUE))
h2 <- ggplot(hyd.long %>% subset(trait %in% select.traits)) +
  facet_wrap(. ~  trait.plot, scales = "free_x", labeller = label_parsed, nrow = 2) +
  geom_col(aes(x = deci_sp.plot, y = value,
               fill = deciduousness),
           position = position_dodge2(width = 0.9, preserve = "single")) +
  theme(axis.text.y = element_text(face = "plain", size = 8)) +
  coord_flip() +
  guides(fill = guide_legend(title = "Deciduousness")) +
  xlab("Deciduousness") + ylab("Value")
ggsave(file.path(figures.folder, paste0("BrettWolfe_traits_vs_deciduousness_sp_bar_HSM.jpeg")),
       plot = h2, height = 7, width = 10, units ='in')
h3 <- h2 %+% subset(hyd.long, !trait %in% select.traits)
ggsave(file.path(figures.folder, paste0("BrettWolfe_traits_vs_deciduousness_sp_bar_CWR.jpeg")),
       plot = h3, height = 7, width = 12, units ='in')

#******************************************************
####----Plot Phenology by Kunert hydraulic traits-----
#******************************************************

ggplot(traits.labels.data %>% subset(trait != "se"), aes(x = deciduousness, y = value)) +
  facet_wrap(. ~  trait.plot, scales = "free_y", labeller = label_parsed) +
  geom_text(aes(label = groups), vjust = 1, hjust = 0, show.legend = FALSE) +
  geom_boxplot(data = traits.long %>% subset(trait != "se"), aes(fill = deciduousness), stat = "boxplot", notch = TRUE) +
  geom_jitter(data = traits.long %>% subset(trait != "se"), width = 0.05, shape = 21, fill = "darkgray",
              color = "black", show.legend = FALSE, alpha = 0.7) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 1, hjust = 1)) +
  scale_color_brewer(palette = "Greens", direction = -1) +
  scale_fill_brewer(name = "Deciduousness", palette = "Greens", direction = -1)
ggsave(file.path(figures.folder, paste0("kunert_traits_vs_deciduousness.jpeg")), height = 8, width = 11, units ='in')

## species wise for sp in hyd.traits

t2 <- ggplot(traits.long.hyd %>% subset(sp %in% unique(hyd$sp) & trait != "se")) +
  facet_wrap(. ~  trait.plot, scales = "free_x", labeller = label_parsed, nrow = 2) +
  geom_col(aes(x = deci_sp.plot, y = value,
               fill = deciduousness),
           position = position_dodge2(width = 0.9, preserve = "single")) +
  theme(axis.text.y = element_text(face = "plain", size = 8)) +
  coord_flip() +
  guides(fill = guide_legend(title = "Deciduousness")) +
  xlab("Deciduousness") + ylab("Value")
ggsave(file.path(figures.folder, paste0("Kunert_traits_vs_deciduousness_sp_bar.jpeg")),
       plot = t2, height = 7, width = 14, units ='in')

## Thus evergreen species have:
# significantly higher median TLP than F. Deci,
# significantly lower Leaf Kmax than F Deci and Obligate deci
# significantly higher SPAD than F. Deci or O. Deci
# higher lwp.min, greater moist site preference
## Evg, BDeci & FDeci have significantly higher Chl and WD than O Deci

#******************************************************
## Mortality vs growth rate by phenology-----
#******************************************************

growth.type <- "med"
formula = y ~ x
demo.sp_size <- demo.sp_size %>%
  left_join(deci, by = "sp")
p.d0 <- ggplot(demo.sp_size %>% subset(size != "NA"),
               aes(y = mrate, x = grate)) +
  geom_point() +
  facet_wrap(. ~ size, scales = "free_y") +
  geom_smooth(method = "lm") +
  scale_x_log10() + scale_y_log10() +
  xlab(expression("Mean Growth Rate (mm/yr)")) +  ylab(expression("Mean Mortality Rate (% per year)"))
p.d0.1 <- p.d0 +  stat_poly_eq(npcx = 0.1, npcy = 0.2, size = 4, aes(label = paste(..rr.label..)), rr.digits = 2, formula = formula, parse = TRUE) +
  stat_fit_glance(npcx = 0.1, npcy = 0.1, size = 4, method = 'lm',  method.args = list(formula = formula), geom = 'text_npc', aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")))
ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_mean_Growth_vs_mrate.jpeg")),
       plot = p.d0.1, height = 6, width = 6, units='in')

p.d1 <-  p.d0 + facet_wrap(size ~ deciduousness, scales = "free_y") + theme(axis.text.x = element_text(angle = 90)) +
  stat_poly_eq(npcx = 0.05, npcy = 0.2, size = 4, aes(label = paste(..rr.label..)), rr.digits = 2, formula = formula, parse = TRUE) +
  stat_fit_glance( npcx = 0.05, npcy = 0.1, size = 4, method = 'lm',  method.args = list(formula = formula), geom = 'text_npc', aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")))
ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_mean_Growth_vs_mrate_Deci.jpeg")),
       plot = p.d1, height = 9, width = 10, units='in')
p.d2 <-  p.d0 %+% subset(demo.sp_size, avg.abund >= n.threshold  &
                           !deciduousness %in% c("Obligate Deciduous")) +
  facet_wrap(. ~ deciduousness, scales = "free_y") + theme(axis.text.x = element_text(angle = 90)) +
  stat_poly_eq(npcx = 0.85, npcy = 0.2, size = 4, aes(label = paste(..rr.label..)), rr.digits = 2, formula = formula, parse = TRUE) +
  stat_fit_glance( npcx = 0.85, npcy = 0.1, size = 4, method = 'lm',  method.args = list(formula = formula), geom = 'text_npc', aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")))+
  ggtitle("Large Size class (>= 30 cm DBH)") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_mean_Growth_vs_mrate_Deci_avg.abund_above", n.threshold,"_large.jpeg")),
       plot = p.d2, height = 6, width = 6, units='in')

#******************************************************
### Is leaf phenology linked to vulnerability to different drought intensity and duration?------
#******************************************************
y.label.1 <- expression(atop(Mortality[Interval], '-'~Mean[Mortality]~('%'*yr^{-1})))

m1 <- ggplot(mrate.long %>%
         subset(!is.na(size) & avg.abund >= n.threshold & !is.na(deciduousness)),
       aes(x = deciduous, y = diff.mrate, color = avg.abund)) +
  scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
                       low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
  facet_grid(size.num ~ censusint.m, scales = "free_y") +
  geom_hline(aes(yintercept = 0), color = "blue", size = 0.5) +
  geom_boxplot(aes(fill = deciduousness.label), stat = "boxplot", notch = TRUE) +
  scale_fill_brewer(name = "", palette = "Greens", direction = -1) +
  theme(legend.position = "top") +
  ylab(y.label.1) + xlab("Deciduousness") +
  ggtitle("Mortality Anomaly by Leaf Phenology") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_deci_by_size.jpeg")), plot = m1, height = 8, width = 9, units='in')
m2 <- m1 %+% subset(mrate.long, size == "large" & avg.abund >= n.threshold & !is.na(deciduousness))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_deci_by_size_large.jpeg")), plot = m2, height = 3, width = 9, units='in')
m2 <- m1 %+% subset(mrate.long, size == "large" & avg.abund >= n.threshold & !is.na(deciduousness)) +
  geom_jitter(width = 0.05, shape = 21, fill = "darkgray", color = "black", show.legend = FALSE, alpha = 0.7)
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_deci_by_size_large_points.jpeg")), plot = m2, height = 3, width = 9, units='in')

mrate.long.hyd <- subset(mrate.long, sp %in% hyd$sp) %>%
  left_join(hyd.wide %>% select(sp, p88S, HSMTLP.88S, HSM88S, HSM50S), by = "sp") %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(HSMTLP.88S)]), ordered=TRUE))
# show_col(viridis_pal()(4))
m3.base <- ggplot(subset(mrate.long.hyd, size == "large" & sp %in% c("cordal", "luehse", "tab1ro"))) +
  facet_grid(size.num ~ censusint.m, scales = "free_y") +
  theme(legend.position = "top") +
  ylab(y.label.1) + xlab("Deciduousness") +
  guides(fill = guide_legend(title = "Deciduousness")) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y = element_text(face = "plain", size = 8))
m3.1 <- m3.base + geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Species leafless in early wet season, with increasing HSM '*Psi[TLP]*' - '*italic('P')['88,Stem']))
  #  does not work: fill = "Facultative Deciduous"
  # guides(fill = guide_legend(title = "Deciduousness",
  #                            override.aes = list(fill = c("Facultative Deciduous" = "#35B779FF"))))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_deci_HSMTLP.88S_increasing_spp_leafless in early wet season.jpeg")),
       plot = m3.1, height = 3, width = 9, units='in')

m4.1 <- m3.base %+% subset(mrate.long.hyd, sp %in% hyd$sp & size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Mortality for Evergreen Species with increasing HSM '*Psi[TLP]*' - '*italic('P')['88,Stem']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_HSMTLP.88S_increasing_evergreens.jpeg")),
       plot = m4.1, height = 3, width = 9, units='in')

mrate.long.hyd <- mrate.long.hyd %>% mutate(sp.plot = factor(sp, levels=unique(sp[order(-p88S)]), ordered=TRUE))
m3.2 <- m3.base + geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Species leafless in early wet season, with increasingly more negative '*italic('P')['88,Stem']))
ggsave(file.path(paste0(figures.folder, "/sp_Mortality_rate_by_period_deci_p88S_increasing_spp_leafless in early wet season.jpeg")),
       plot = m3.2, height = 3, width = 9, units='in')

m4.2 <- m3.base %+% subset(mrate.long.hyd, sp %in% hyd$sp & size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Mortality for Evergreen Species with increasingly more negative '*italic('P')['88,Stem']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_p88S_increasing_evergreens.jpeg")),
       plot = m4.2, height = 3, width = 9, units='in')
mrate.long.hyd <- mrate.long.hyd %>% mutate(sp.plot = factor(sp, levels=unique(sp[order(-HSM50S)]), ordered=TRUE))
m4.3 <- m3.base %+% subset(mrate.long.hyd, sp %in% hyd$sp & size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Mortality for Evergreen Species with increasing HSM '*Psi[min]*' - '*italic('P')['50,Stem']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_HSMLWP.50S_increasing_evergreens.jpeg")),
       plot = m4.3, height = 3, width = 9, units='in')

mrate.long.traits <- subset(mrate.long, sp %in% traits$sp) %>%
  left_join(traits.wide %>% select(-deciduous, -deciduousness), by = "sp") %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(KmaxL)]), ordered=TRUE))

m4.4 <- m3.base %+% subset(mrate.long.traits, size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Mortality for Evergreen Species with increasing '*italic('K')['max, Leaf']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_KmaxL_increasing_evergreens.jpeg")),
       plot = m4.4, height = 3, width = 9, units='in')

mrate.long.traits <- mrate.long.traits %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(-HSMLWP.50L)]), ordered=TRUE))
m4.5 <- m3.base %+% subset(mrate.long.traits, size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Mortality for Evergreen Species with increasing HSM '*Psi[min]*' - '*italic('P')['50, Leaf']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_HSMLWP.50L_increasing_evergreens.jpeg")),
       plot = m4.5, height = 3, width = 9, units='in')

## all evergreens
mrate.long <- mrate.long %>% mutate(sp.plot = factor(sp, levels=unique(sp[order(mean.mrate)]), ordered=TRUE))
m5 <- m3.base %+% subset(mrate.long, size == "large" & deciduous == "E") +
  geom_col(aes(x = sp, y = diff.mrate, fill = deciduousness)) +
  ggtitle("Mortality for Evergreen Species") +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 0.5, hjust = 1, size = 4))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_evergreens.jpeg")),
       plot = m5, height = 3, width = 15, units='in')
mrate.long <- mrate.long %>% mutate(sp.plot = factor(sp, levels=unique(sp[order(mean.mrate)]), ordered=TRUE))

m6 <- m3.base %+% subset(mrate.long, size == "large" & deciduous == "DF") +
  geom_col(aes(x = sp, y = diff.mrate, fill = deciduousness)) +
  ggtitle("Mortality for Facultative Deciduous Species") +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 0.5, hjust = 1, size = 4))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_Facultative Deciduous.jpeg")),
       plot = m6, height = 3, width = 15, units='in')

#******************************************************
### Is leaf phenology linked to growth vulnerability to different drought intensity and duration?------
# Indeed facultative deciduous species show greater reduction in growth in 2005-2010 period two successive early wet seasons were dry
# a large chunk of their limited growing period
#******************************************************
y.label.2 <- expression(atop(Std.~Growth[interval], '-'~Mean[Std.~Growth]))

g1 <- ggplot(grate.long %>%
               subset(!is.na(size) & !is.na(deciduousness)),
             aes(x = deciduous, y = median)) +
  scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
                       low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
  facet_grid(size.num ~ censusint.m, scales = "free_y") +
  geom_hline(aes(yintercept = 0), color = "blue", size = 0.5) +
  geom_boxplot(aes(fill = deciduousness.label), stat = "boxplot", notch = TRUE) +
  scale_fill_brewer(name = "", palette = "Greens", direction = -1) +
  theme(legend.position = "top") +
  ylab(y.label.2) + xlab("Deciduousness") +
  ggtitle("Growth Rates by Leaf Phenology") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_deci_by_size.jpeg")), plot = g1, height = 8, width = 9, units='in')
g2 <- g1 %+% subset(grate.long, size == "large" & !is.na(deciduousness))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_deci_by_size_large.jpeg")), plot = g2, height = 3, width = 9, units='in')
g2 <- g1 %+% subset(grate.long, size == "large" & !is.na(deciduousness)) +
  geom_jitter(width = 0.05, shape = 21, fill = "darkgray", color = "black", show.legend = FALSE, alpha = 0.7)
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_deci_by_size_large_points.jpeg")), plot = g2, height = 3, width = 9, units='in')

hyd.wide <- hyd.long %>% pivot_wider(names_from = trait, values_from = value, -trait.plot)
grate.long.hyd <- subset(grate.long, sp %in% hyd$sp) %>%
  left_join(hyd.wide %>% select(sp, p88S, HSMTLP.88S, HSM88S), by = "sp") %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(HSMTLP.88S)]), ordered=TRUE))

# show_col(viridis_pal()(4))
g3.base <- ggplot(subset(grate.long.hyd, size == "large" & sp %in% c("cordal", "luehse", "tab1ro"))) +
  facet_grid(size.num ~ censusint.m, scales = "free_y") +
  theme(legend.position = "top") +
  ylab(expression("Growth Rate - Mean (% per year)")) + xlab("Deciduousness") +
  # ylab(expression('Growth Rate'[interval]*' - '*italic('E')'[Growth Rate'[interval]'] (% yr'^-1')')) +
  # xlab("Deciduousness") +
  guides(fill = guide_legend(title = "Deciduousness")) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 0.5, hjust = 1))
g3.1 <- g3.base + geom_col(aes(x = sp.plot, y = median, fill = deciduousness)) +
  ggtitle(expression('Species leafless in early wet season, with increasing HSM '*Psi[TLP]*' - '*italic('P')['88,Stem']))
#  does not work: fill = "Facultative Deciduous"
# guides(fill = guide_legend(title = "Deciduousness",
#                            override.aes = list(fill = c("Facultative Deciduous" = "#35B779FF"))))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_deci_HSMTLP.88S_increasing_spp_leafless in early wet season.jpeg")),
       plot = g3.1, height = 4, width = 9, units='in')

g4.1 <- g3.base %+% subset(grate.long.hyd, sp %in% hyd$sp & size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = median, fill = deciduousness)) +
  ggtitle(expression('Growth for Evergreen Species with increasing HSM '*Psi[TLP]*' - '*italic('P')['88,Stem']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_HSMTLP.88S_increasing_evergreens.jpeg")),
       plot = g4.1, height = 4, width = 9, units='in')

grate.long.hyd <- grate.long.hyd %>% mutate(sp.plot = factor(sp, levels=unique(sp[order(-p88S)]), ordered=TRUE))
g3.2 <- g3.base + geom_col(aes(x = sp.plot, y = median, fill = deciduousness)) +
  ggtitle(expression('Species leafless in early wet season, with increasingly more negative '*italic('P')['88,Stem']))
ggsave(file.path(paste0(figures.folder, "/sp_Growth_rate_by_period_deci_p88S_increasing_spp_leafless in early wet season.jpeg")),
       plot = g3.2, height = 4, width = 9, units='in')

g4.2 <- g3.base %+% subset(grate.long.hyd, sp %in% hyd$sp & size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = median, fill = deciduousness)) +
  ggtitle(expression('Growth for Evergreen Species with increasingly more negative '*italic('P')['88,Stem']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_p88S_increasing_evergreens.jpeg")),
       plot = g4.2, height = 4, width = 9, units='in')

## all evergreens
grate.long <- grate.long %>%
  left_join(traits.wide %>% select(sp, TLP, KmaxL), by = "sp") %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(KmaxL)]), ordered=TRUE))

g5 <- g3.base %+% subset(grate.long, size == "large" & deciduous == "E") +
  geom_col(aes(x = sp, y = median, fill = deciduousness)) +
  ggtitle("Growth rates (DBH residuals) for Evergreen Species") +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 0.5, hjust = 1, size = 8))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_evergreens.jpeg")),
       plot = g5, height = 4, width = 15, units='in')

g6 <- g3.base %+% subset(grate.long, size == "large" & deciduous == "DF") +
  geom_col(aes(x = sp, y = median, fill = deciduousness)) +
  ggtitle("Growth rates (DBH residuals) for Facultative Deciduous Species") +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 0.5, hjust = 1, size = 12))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_facultative_deciduous.jpeg")),
       plot = g6, height = 4, width = 15, units='in')

## Plot additional mortality by period mean swp-------

#******************************************************
### Yearly psi dynamics versus climatology-------
#******************************************************
psi.stat.1 <- psi %>%
  group_by(interval.yrs, date, depth) %>%
  summarise(mean = -mean(psi, na.rm = TRUE),
            upper.CI = -quantile(psi, probs = 0.975),
            lower.CI = -quantile(psi, probs = 0.025))
psi.stat.1 <- psi.stat.1 %>%
  group_by(interval.yrs, depth) %>%
  mutate(days = 1:n())
rectangles <- data.frame(
  xmin = as.Date(paste0(c(1990:2018), "-04-01")),
  xmax = as.Date(paste0(c(1990:2018), "-11-01")),
  ymin = 0,
  ymax = 2.5
)
plot.psi.stat.1 <- ggplot(psi.stat.1 %>%
                            subset(depth %in% c(0.21,  0.37,  1.00,  1.70)) %>% droplevels()) +
  geom_rect(data=rectangles, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill='gray80', alpha=0.8) +
  # geom_ribbon(aes(x = date, ymin = lower.CI, ymax = upper.CI), alpha = 0.3) +
  theme(panel.grid.major.y = element_line()) +
  geom_line(aes(x = date, y = mean, group = depth, color = depth), size = 0.3) +
  scale_color_continuous(trans = "reverse", guide = guide_colorbar(title = "Depth\n(cm)")) +
  geom_vline(xintercept = census.beg) + # as.Date(paste0(c(1990:2015), "-01-01")) +
  scale_y_reverse(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 10, 15), limits = c(2.5, 0)) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%Y")) +
  coord_cartesian(xlim = c(as.Date("1990-01-01"), as.Date("2018-12-31"))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  # facet_grid(interval.yrs ~ .) +
  ylab("-Soil Water Potential [MPa]") +
  theme(text = element_text(size = 12)) + xlab("")
  # ggtitle("PSI:Violins of interval means of best-fit ensembles")
ggsave("psi_model_daily_bestfit_params.top.few_CI_full.jpeg", plot = plot.psi.stat.1,
       file.path(figures.folder), device = "jpeg", height = 3, width = 20, units='in')

psi.2 <- psi %>%
  mutate(interval.yrs.to.plot = forcats::fct_explicit_na(cut(date, include.lowest = TRUE, breaks = cut.breaks.2,
                                    labels = cut.labels.2, right = TRUE)))

psi.stat.2 <- psi.2 %>%
  group_by(interval.yrs.to.plot, date, depth) %>%
  summarise(mean = mean(psi, na.rm = TRUE),
            upper.CI = quantile(psi, probs = 0.975),
            lower.CI = quantile(psi, probs = 0.025))
psi.stat.2 <- psi.stat.2 %>%
  group_by(interval.yrs.to.plot, depth) %>%
  mutate(days = 1:n())
rectangles.2 <- data.frame(
  xmin = c(0:4)*365 + 120,
  xmax = c(0:4)*365 + 335,
  ymin = 0,
  ymax = -2.5
)

plot.psi.stat.2 <- ggplot(psi.stat.2 %>%
                            subset(depth %in% c(0.1, 0.21,  0.37,  1.00,  1.70)) %>% droplevels()) +
  geom_rect(data=rectangles.2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill='gray80', alpha=0.8) +
  geom_ribbon(aes(x = days, ymin = lower.CI, ymax = upper.CI), alpha = 0.3, fill = "red") +
  geom_line(aes(x = days, y = mean, group = depth, color = depth), size = 0.3) +
  scale_color_continuous(trans = "reverse", guide = guide_colorbar(title = "Depth\n(cm)")) +
  # scale_y_continuous(breaks = -c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 10, 15)) +
  scale_x_continuous(breaks = c(1:5)*365, labels = paste0("Yr.", c(1:5))) +
  coord_cartesian(ylim = c(-2.5, 0)) +
  theme(panel.grid.major.y = element_line()) +
  facet_grid(interval.yrs.to.plot ~ .) +
  ylab("Soil Water Potential (MPa)") +
  theme(text = element_text(size = 12)) + xlab("Census Year")
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_censuspanels.jpeg", plot = plot.psi.stat.2,
       file.path(figures.folder), device = "jpeg", height = 5, width = 7, units='in')

## by depth panels
plot.psi <- ggplot(psi %>% subset(date < as.Date("1992-12-31") & depth >= 1), aes(x = date, y = -psi)) +
  scale_y_continuous(trans="rev_sqrt") +
  geom_line(aes(group = par.sam, color = as.factor(par.sam)), show.legend = F, size = 0.2) +
  geom_vline(xintercept = census.beg, color = "gray") +
  facet_grid(depth ~ ., scales = "free_y") +
  ylab("-Soil Water Potential [MPa]") + xlab("Date") +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%Y")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Soil Water Potential for best-fit ensembles")
ggsave("psi_model_daily_all_depths_params.top.few_full_until_1992-12-31.jpeg", plot = plot.psi, path =
         file.path(figures.folder), height = 15, width = 8.94, units='in')

psi.stat.3 <- psi %>%
  mutate(doy = format(date, "%j")) %>%
  group_by(doy, depth) %>%
  summarise(mean = mean(psi, na.rm = TRUE),
            sd = sd(psi, na.rm = TRUE),
            upper.CI = quantile(psi, probs = 0.975),
            lower.CI = quantile(psi, probs = 0.025)) %>%
  ungroup(doy, depth) %>%
  mutate(doy = as.numeric(doy))
psi.stat.4 <- psi %>%
  mutate(doy = format(date, "%j"),
         year = format(date, "%Y")) %>%
  ## due to time nneded for model initialisation soil layers deeper than 2.9 m
  ## are not fully recharged until the end of 1992, layer 2.9 m recharges by the end of 1992
  subset(!year %in% c("1990")) %>%
  subset(!year %in% c("1990", "1991") | depth == 2.9) %>%
  subset(!year %in% c("1990", "1991", "1992") | depth < 2.9) %>%
  group_by(date, doy, year, interval.yrs, depth) %>%
  summarise(mean = mean(psi, na.rm = TRUE),
            median = median(psi, na.rm = TRUE),
            sd = sd(psi, na.rm = TRUE),
            upper.CI = quantile(psi, probs = 0.975),
            lower.CI = quantile(psi, probs = 0.025)) %>%
  ungroup(doy, year, depth) %>%
  mutate(doy = as.numeric(doy))
psi.stat.5 <- psi.stat.4 %>%
  group_by(doy, depth) %>%
  summarise(mean.clim = mean(mean, na.rm = TRUE),
            sd.clim = sd(mean, na.rm = TRUE),
            upper.CI.clim = quantile(mean, probs = 0.975),
            lower.CI.clim = quantile(mean, probs = 0.025)) %>%
  ungroup(doy, depth) %>%
  mutate(doy = as.numeric(doy))
rectangles.3 <- data.frame(
  xmin = 120,
  xmax = 335,
  ymin = 0,
  ymax = -2.5
)
source("code/Utilities/plot.ticks.R")
plot.psi.stat.5.base <- ggplot(psi.stat.5 %>% droplevels()) +
  # geom_rect(data=rectangles.3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
  #           fill='gray80', alpha=0.8) +
  scale_x_continuous(breaks = c(seq(0, 360, by = 60))) +
  coord_cartesian(ylim = c(-2.5, 0)) +
  theme(panel.grid.major.y = element_line()) +
  ylab("Soil Water Potential (MPa)") + xlab("Day of the Year") +
  theme(text = element_text(size = 12))
plot.psi.stat.5 <- plot.psi.stat.5.base +
  geom_ribbon(aes(x = doy, ymin = lower.CI.clim, ymax = upper.CI.clim), alpha = 0.3, fill = "grey20") +
  geom_line(aes(x = doy, y = mean.clim, group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
  guides(color = guide_legend(title = "Depth(m)", order = 2, override.aes = list(size = 3)))
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_censuspanels_climatology.jpeg",
       plot = plot.ticks(plot.psi.stat.5),
       file.path(figures.folder), device = "jpeg", height = 4, width = 7, units='in')

plot.psi.stat.5.over <- plot.psi.stat.5.base %+%
  subset(psi.stat.5, depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)) +
  geom_line(aes(x = doy, y = mean.clim, linetype = "climatology", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
  geom_line(data = subset(psi.stat.4, year == "2016" & depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)),
            aes(x = doy, y = mean, linetype = "2016", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
  guides(color = guide_legend(title = "Depth(m)", legend.position = "right", order = 1, override.aes = list(size = 3)),
         linetype = guide_legend(order = 2, title = NULL, legend.position = "top", override.aes =
                                   list(linetype = c("climatology" = "dashed", "2016" = "solid")))) +
  coord_cartesian(ylim = c(-3, 0)) + ggtitle("2016")
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_censuspanels_climatology_over.jpeg",
       plot = plot.ticks(plot.psi.stat.5.over),
       file.path(figures.folder), device = "jpeg", height = 4, width = 7, units='in')

pdf(paste0(figures.folder, "/psi_model_daily_bestfit_params.top.few_CI_full_censuspanels_climatology_over_by_year.pdf"), height = 4, width = 7)
for (i in unique(psi.stat.4$year)) {
  plot.psi.stat.5.yr <- plot.psi.stat.5.base %+% subset(psi.stat.5, depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)) +
    geom_line(aes(x = doy, y = mean.clim, linetype = "climatology", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
    geom_line(data = subset(psi.stat.4, year == i & depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)),
              aes(x = doy, y = mean, linetype = "Year", group = as.factor(depth), color = as.factor(depth)), size = 0.5) +
    guides(color = guide_legend(title = "Depth(m)", order = 1, override.aes = list(size = 3)),
           linetype = guide_legend(order = 2, title = NULL, override.aes =
                                     list(linetype = c("climatology" = "solid", "Year" = "dashed")))) +
    coord_cartesian(ylim = c(-3, 0), xlim = c(0, 200)) + ggtitle(i)
  print(plot.psi.stat.5.yr)
}
dev.off()
#******************************************************
### Census Interval psi dynamics versus climatology-------
#******************************************************
psi.stat.4 <- psi.stat.4 %>%
  mutate(interval.yrs.2 = forcats::fct_explicit_na(cut(date, include.lowest = TRUE, breaks = c(cut.breaks, max(psi.stat.4$date, na.rm = TRUE)),
    labels = c(cut.labels.2, "2015-2018"), right = TRUE))) %>%
  left_join(psi.stat.5, by = c("doy", "depth")) %>%
  mutate(mean.below.sd1 = ifelse(mean < c(mean.clim - sd.clim), mean, NA),
         mean.below.sd2 = ifelse(mean < c(mean.clim - 2*sd.clim), mean, NA),
         mean.below.lower.CI = ifelse(mean < lower.CI.clim, mean, NA))

pct.drought.days <- psi.stat.4 %>%
  mutate(season =  ifelse(doy < 120, "Dry Season", "Wet Season")) %>%
  group_by(depth, interval.yrs, interval.yrs.2, season) %>%
  summarise(pct.days.below.sd1.clim = 100*round(sum(!is.na(mean.below.sd1))/n(), 3),
         pct.days.below.sd2.clim = 100*round(sum(!is.na(mean.below.sd2))/n(), 3),
         pct.days.below.lower.CI.clim = 100*round(sum(!is.na(mean.below.lower.CI))/n(), 3),
         pct.days.above.sd1.clim = 100 - pct.days.below.sd1.clim,
         pct.days.above.sd2.clim = 100 - pct.days.below.sd2.clim,
         pct.days.above.lower.CI.clim = 100 - pct.days.below.lower.CI.clim,
         pct.days.below.lower.CI.clim.0.5 = 100*round(sum(!is.na(mean.below.lower.CI) & mean.below.lower.CI < -0.5)/n(), 3),
         pct.days.above.lower.CI.clim.0.5 = 100 - pct.days.below.lower.CI.clim.0.5) %>%
  ungroup(depth, interval.yrs) %>% #subset(depth %in% c(0.06, 0.62, 1))
  mutate(depth.fac = factor(depth, levels = sort(unique(psi.stat.4$depth), decreasing = TRUE))) %>%
  mutate(int.ssn = paste(interval.yrs, season, sep = "_"))

heat.fr.drought.days.base <- ggplot(pct.drought.days %>% subset(interval.yrs.2 != "(Missing)"),
                   aes(x = interval.yrs.2, y = depth.fac)) +
  facet_wrap(~season, nrow = 2) +
  xlab("Census Interval") + ylab("Depth (m)") +
  scale_fill_viridis_c(expression('% Days '*Psi['Soil, DOY']*' < '*italic(Q)[italic(p)*'=0.975, '*Psi['Soil, DOY']]),
                       direction = -1, option = "plasma") +
  theme(legend.position = "top", legend.direction = "horizontal")
heat.fr.drought.days <- heat.fr.drought.days.base +
  geom_tile(aes(fill = pct.days.below.lower.CI.clim))
ggsave("pct.days.below.lower.CI.clim_by depth & intervalyrs&ssn_full.jpeg",
       plot = heat.fr.drought.days, file.path(figures.folder), device = "jpeg", height = 8, width = 6, units='in')

heat.fr.drought.days.study <- heat.fr.drought.days %+% subset(pct.drought.days, interval.yrs != "(Missing)") +
  geom_tile(aes(fill = pct.days.below.lower.CI.clim))
ggsave("pct.days.below.lower.CI.clim_by depth & intervalyrs&ssn_study_period.jpeg",
       plot = heat.fr.drought.days.study, file.path(figures.folder), device = "jpeg", height = 8, width = 5, units='in')

heat.fr.drought.days.0.5 <- heat.fr.drought.days.base +
  geom_tile(aes(fill = pct.days.below.lower.CI.clim.0.5))
ggsave("pct.days.below.lower.CI.clim.0.5_by depth & intervalyrs&ssn_full.jpeg",
       plot = heat.fr.drought.days.0.5, file.path(figures.folder), device = "jpeg", height = 8, width = 6, units='in')

heat.fr.drought.days.study.0.5 <- heat.fr.drought.days %+% subset(pct.drought.days, interval.yrs != "(Missing)") +
  geom_tile(aes(fill = pct.days.below.lower.CI.clim.0.5))
ggsave("pct.days.below.lower.CI.clim.0.5_by depth & intervalyrs&ssn_study_period.jpeg",
       plot = heat.fr.drought.days.study.0.5, file.path(figures.folder), device = "jpeg", height = 8, width = 5, units='in')

rectangles.4 <- data.frame(
  xmin = 120,
  xmax = 335,
  ymin = 0,
  ymax = -3.0
)

## add legends for ribbon areas for drought rarity, lines for climatology,
# Could add a line for sd2
plot.psi.stat.6.interval.base <- plot.psi.stat.5.base %+%
  subset(psi.stat.4, depth %in% c(0.06, 0.62, 1)) +
  facet_wrap(. ~ interval.yrs.2) +
  geom_rect(data=rectangles.4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill='gray80', alpha=0.8) +
  geom_line(aes(x = doy, y = mean.clim, group = as.factor(depth), color = as.factor(depth)), size = 0.3, linetype = "solid") +
  geom_ribbon(aes(x = doy, ymin = mean.below.sd2, ymax = mean.clim, group = as.factor(depth_year),
                  fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE) +
    theme(panel.grid.major.y = element_line(size = 0.1)) +
  guides(color = guide_legend(title = "Depth(m)", order = 1, override.aes = list(size = 3))) +
  coord_cartesian(ylim = c(-3, 0), xlim = c(0, 200)) +
  theme(legend.position = "top")
plot.psi.stat.6.interval.sd1 <- plot.psi.stat.6.interval.base  +
  geom_ribbon(aes(x = doy, ymin = mean.below.sd1, ymax = mean.clim, group = as.factor(depth_year),
                fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE)
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_sd1.jpeg",
       plot = plot.psi.stat.6.interval.sd1, file.path(figures.folder), device = "jpeg", height = 4, width = 7, units='in')

plot.psi.stat.7.interval.sd1 <- plot.psi.stat.6.interval.base %+%
  subset(psi.stat.4, depth %in% c(0.06, 0.62, 1) & interval.yrs != "(Missing)") +
  facet_wrap(. ~ interval.yrs, nrow = 1) +
  geom_ribbon(aes(x = doy, ymin = mean.below.sd1, ymax = mean.clim, group = as.factor(depth_year),
                  fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE)
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_study_period_sd1.jpeg",
       plot = plot.psi.stat.7.interval.sd1, file.path(figures.folder), device = "jpeg", height = 2.5, width = 7, units='in')

plot.psi.stat.6.interval.sd2 <- plot.psi.stat.6.interval.base  +
  geom_ribbon(aes(x = doy, ymin = mean.below.sd2, ymax = mean.clim, group = as.factor(depth_year),
                  fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE)
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_sd2.jpeg",
       plot = plot.psi.stat.6.interval.sd2, file.path(figures.folder), device = "jpeg", height = 4, width = 7, units='in')
plot.psi.stat.7.interval.sd2 <- plot.psi.stat.6.interval.base %+%
  subset(psi.stat.4, depth %in% c(0.06, 0.62, 1) & interval.yrs != "(Missing)") +
  facet_wrap(. ~ interval.yrs, nrow = 1) +
  geom_ribbon(aes(x = doy, ymin = mean.below.sd2, ymax = mean.clim, group = as.factor(depth_year),
                  fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE)
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_study_period_sd2.jpeg",
       plot = plot.psi.stat.7.interval.sd2, file.path(figures.folder), device = "jpeg", height = 2.5, width = 7, units='in')

## plot over observed
plot.psi.stat.7.interval.mean <- plot.psi.stat.6.interval.base %+%
  subset(psi.stat.4, depth %in% c(0.06, 0.62, 1) & interval.yrs.2 != "2015-2020") +
  facet_wrap(. ~ interval.yrs, nrow = 1) +
  geom_ribbon(aes(x = doy, ymin = mean, ymax = mean.clim, group = as.factor(depth_year),
                  fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE)
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_study_period_all_observed.jpeg",
       plot = plot.psi.stat.7.interval.mean, file.path(figures.folder), device = "jpeg", height = 2.5, width = 7, units='in')


plot.list.grid <- list()
for (i in 1:length(unique(psi.stat.4$year))) {
  year.on = unique(psi.stat.4$year)[i]
  plot.psi.stat.5.yr <-  ggplot(psi.stat.5 %>% subset(depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)) %>%
                                  droplevels()) +
    # geom_rect(data=rectangles.3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
    #           fill='gray80', alpha=0.8) +
    scale_x_continuous(breaks = c(seq(0, 360, by = 60))) +
    coord_cartesian(ylim = c(-2.5, 0)) +
    theme(panel.grid.major.y = element_line()) +
    ylab("Soil Water Potential (MPa)") + xlab("Day of the Year") +
    theme(text = element_text(size = 12)) +
    # geom_ribbon(aes(x = doy, ymin = lower.CI, ymax = upper.CI), alpha = 0.3, fill = "grey20") +
    geom_line(aes(x = doy, y = mean.clim, linetype = "climatology", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
    geom_line(data = subset(psi.stat.4, year == year.on & depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)),
              aes(x = doy, y = mean, linetype = "Year", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
    guides(color = guide_legend(title = "Depth(m)", order = 1, override.aes = list(size = 3)),
           linetype = guide_legend(order = 2, title = NULL, override.aes =
                                     list(linetype = c("climatology" = "solid", "Year" = "dashed")))) +
    coord_cartesian(ylim = c(-3, 0), xlim = c(0, 200)) +
      ggtitle(year.on) + xlab("") + ylab("")
  if(year.on != "2018") {
    plot.psi.stat.5.yr <- plot.psi.stat.5.yr + theme(legend.position = "none")
  }
  # if (as.numeric(year.on)%%5 != 0) {
  #   plot.psi.stat.5.yr <- plot.psi.stat.5.yr + ylab("")
  # }
  # if (as.numeric(year.on)[i] %in% c(2015:2018)) {
  #   plot.psi.stat.5.yr <- plot.psi.stat.5.yr + xlab("")
  # }
  plot.list.grid[[i]] <- plot.psi.stat.5.yr
}

plot.list <- lapply(plot.list.grid, plot.ticks)
# pdf("all.pdf")
# invisible(lapply(plot.list, print))
# dev.off()

ggsave(file.path(figures.folder, "arrange.pdf"), arrangeGrob(grobs = plot.list.grid, ncol = 5), width = 15, height = 10, units = "in")

### Does % of the growing season < p50 or p80 explain inter-census growth or mortality?
### Use Meinzer data for rooting depths, and assign the depth whose water potential to track


head(psi.stat.4)
