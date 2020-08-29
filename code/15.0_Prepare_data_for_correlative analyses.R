
rm(list=ls())

if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, readxl, forcats, agricolae, gridExtra,
               scales, GGally, ggpmisc, Evapotranspiration,
               data.table, bci.elm.fates.hydro, mgcv, lubridate, smooth)
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
range01 <- function(x){(x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))}
range01.p <- function(x){(x - pmin(x, na.rm = TRUE))/(pmax(x, na.rm = TRUE)-pmin(x, na.rm = TRUE))}

indicator <- function(x, I.threshold, greater.than = TRUE) {
  if(greater.than == TRUE) {
    result <- ifelse((x > I.threshold), 1, 0)
  } else {
    result <- ifelse((x < I.threshold), 1, 0)
  }
  return(result)
}
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

rev_sqrt_trans <- function() {
  scales::trans_new(
    name = "rev_sqrt",
    transform = function(x) -sqrt(abs(x)),
    inverse = function(x) x^2);
}

reverselog_trans <- function(base = exp(1)) {
  scales::trans_new(name = paste0("reverselog-", format(base)),
                    log_breaks(base = base),
                    domain = c(1e-100, Inf),
                    transform = function(x) -log(x, base),
                    inverse = function(x) base^(-x))
}
## to calculate tlp based psi thresholds

# X1=sum(min(0, psi-psi_threshold)*PET), max/mean
# X2=sum(I(psi<psi_threshold)), max/mean
# X3=sum(min(0, psi-psi_threshold)), max/mean
# X4=sum(I(psi<psi_threshold)*PET),max/mean

## funciton used to estimate K_leaf from Psi and species specific parameters
Exponential <- function (A, B, psi) {
  A * exp(-B * psi)
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
## Load Isotopic data-----
#******************************************************

load(file = "data-raw/traits/isotopes/Oecologia 1995 Jackson_Fig3_Fig4_& Meinzer 1999_Fig4.Rdata")
leafless_mar.apr <- read.csv("data-raw/traits/isotopes/Meinzer_1999_isotope_sp_leafless_in_mar_april.csv")
# load(file = "results/all_isotopic_record.Rdata")
iso.2.raw <- read.csv("data-raw/traits/isotopes/Meinzer1999_Xylem_Sap_deltaD_March97_DBH_Fig5B.csv", na.strings = c("NA", ""), header = T, row.names = NULL, check.names = F)

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
save(deci, file = file.path(results.folder, "deci.prepped.Rdata"))

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
# View(lwp.diff)
lwp.min.wide <- lwp.min.diurnal %>% pivot_wider(names_from = time, values_from = c(lwp.min, lwp.se, lwp.diff)) %>%
  mutate(lwp.diff = lwp.diff_Diurnal) %>% select(-lwp.diff_Diurnal, -lwp.diff_Predawn)

lwp.min <- lwp.min.diurnal %>%
  left_join(deci, by = "sp")

save(lwp.diff, file = file.path(results.folder, "lwp.diff.Rdata"))
save(lwp.min, file = file.path(results.folder, "lwp.min.Rdata"))
save(lwp.min.wide, file = file.path(results.folder, "lwp.min.wide.Rdata"))

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
load("results/adult.mrate.long.RData")
load(file.path("results/GLUEsetup_part2_BCI.RData"))
## growth rates when dbh.residuals = "on" are residuals from a dbh mixed effects model (for spp) of
## growth. A median residual for each sp_size is calculated only when at least data from
# 3 trees are present across all census intervals.
# Medians within sp_size are then centered and scaled. {residual - E(residual)/sd(residual)}

if(growth_by_si.info$dbh.residuals == "on"){
  growth <- growth_by_si.info$growth
}

grate.long <- dplyr::bind_rows(growth, .id = 'sp_size')
## to add absolute growth
intervals <- 5
growth.selection <- "size_class_predefined_cc_scaled"
load(file = paste0("results/gro.long.cc.med_", intervals, "_", growth.selection, ".Rdata"))

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
adult.mrate.long <- adult.mrate.long %>%
  left_join(deci, by = "sp") %>%
  mutate(censusint.m = recode(census, `1985` = "1982-85", `1990` = "1985-90", `1995` = "1990-95", `2000` = "1995-00", `2005` = "2000-05", `2010` = "2005-10", `2015` = "2010-15")) %>%
  # left_join(demo.sp_size %>%
  #             subset(size %in% c("medium", "large")) %>%
  #             group_by(sp) %>%
  #             summarize(mean.mrate = mean(mrate, na.rm = TRUE),
  #                       mean.grate = mean(grate, na.rm = TRUE)), by = "sp") %>%
  left_join(demo.sp %>%
              mutate(mean.mrate = mrate.adult, na.rm = TRUE,
                     mean.grate = grate.adult, na.rm = TRUE,
                     grate.se = g.adult.se, na.rm = TRUE) %>%
              select(sp, mean.mrate, mean.grate, grate.se), by = "sp") %>%
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
  ungroup(sp_size) %>%
  left_join(gro.long.cc.med %>% select(sp_size, interval, med.dbh.resid, med.growth),
            by = c("sp_size", "interval"))

save(grate.long, file = file.path(results.folder, "grate.long_by_species-size_deciduousness.Rdata"))
save(mrate.long, file = file.path(results.folder, "mrate.long_by_species-size_deciduousness.Rdata"))
save(adult.mrate.long, file = file.path(results.folder, "adult.mrate.long_by_species-size_deciduousness.Rdata"))

#******************************************************
### Load Psi from ELM-FATES-------
#******************************************************

census.meds <- readr::read_rds("results/census.mediandates.rds")
census.beg <- census.meds[3: length(census.meds)]
cut.breaks <- census.beg
cut.breaks.2 <- as.Date(paste0(seq(1990, 2015, by = 5), "-01-01"))
cut.labels.2 <- paste0(seq(1990, 2010, by = 5), "-", seq(1995, 2015, by = 5))
cut.labels.interval <- 3: (length(census.meds)-1)
psi <- bci.elm.fates.hydro::psi
gpp <- bci.elm.fates.hydro::gpp

psi <- psi %>%
  mutate(interval.yrs = forcats::fct_explicit_na(cut(date, include.lowest = TRUE, breaks = cut.breaks,
                                                     labels = cut.labels.2, right = TRUE)))

# new.depths = data.frame(depth =
#                           seq(from = 0.1, to =
#                                 psi$depth[length(psi$depth)], by = 0.1))
#
# psi.date <- split(psi %>% select(-interval.yrs), f = list(psi$date, psi$par.sam), drop = TRUE)
#
# psi.interp.approx <- function(df) {
#   x <- df$depth
#   y <- df$psi
#   xout <- new.depths$depth
#   yout <- predict(interpSpline(x, y), xout)
#   # yout <- approx(x, y, xout, method = "linear")
#   df.1 <- data.frame(date = df$date[1], par.sam = df$par.sam[1],
#                                depth = yout$x, psi = yout$y)
#   # plot(depth ~ psi, data = df.1)
#   return(df.1)
# }
#
# psi.int <- lapply(lapply(psi.date, psi.interp.approx),
#                         as.data.frame) %>%
#   bind_rows()
# save(psi.int,
#      file = file.path("results/psi.interpolated.depths.Rdata"))
# load(file = file.path("results/psi.interpolated.depths.Rdata"))
#
# depth.breaks.1 <- c(0, 0.5, 1, c(1:13 + 0.5))
# depth.labels.1 <- c(0.25, 0.75, 1.25, 2:13)
# psi.depths <- as.data.table(psi.int)[,
#   depths := as.numeric(as.character(cut(depth, breaks = depth.breaks.1, labels = depth.labels.1)))][,
#   keyby = .(date, depths, par.sam), .(psi = mean(psi, na.rm = TRUE))][,]
# psi.depths <- psi.depths %>% as.data.frame() %>% rename(depth = depths) %>%
#   mutate(interval.yrs = forcats::fct_explicit_na(cut(date, include.lowest = TRUE, breaks = cut.breaks,
#                                                      labels = cut.labels.2, right = TRUE)))
#
# save(psi.depths,
#      file = file.path(results.folder, "psi.interpolated.depths_layers_combined.Rdata"))

load(file = file.path(results.folder, "psi.interpolated.depths_layers_combined.Rdata"))
save(psi,
     file = file.path(results.folder, "psi.prepped.Rdata"))

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

clim.dat.yr <- clim.dat %>% subset(!is.na(Year) & !Year %in% c(1985, 2019)) %>%
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
          file.path(results.folder, "clim.daily_with_pet.PM.csv"),
          row.names = FALSE)
save(clim.daily,
     file = file.path(results.folder, "clim.daily_with_pet.PM.Rdata"))

#******************************************************
### GPP vs. RAD VPD or PET models------
#******************************************************

# # GPP from ELM.FATES output
# gpp.rel.daily <- gpp %>% group_by(date) %>%
#   summarise(gpp = median(value)) %>%
#   left_join(clim.daily %>% select(date, VPD, pet.PM, Rs), by = "date")

# GPP Observed
figures.folder.gpp <- paste0("figures/PhenoDemoTraitsPsi/gpp_by_env_variables")
if(!dir.exists(file.path(figures.folder.gpp))) {dir.create(file.path(figures.folder.gpp))}

bci.tower <- read.csv("data-raw/BCI_v3.1.csv", header = TRUE)
bci.tower <- as.data.frame(bci.tower[-1, ])
bci.tower$datetime <- strptime(bci.tower$date, format = "%m/%d/%Y %H:%M")
bci.tower$datetime <- as.POSIXct(bci.tower$datetime, format = "%Y-%m-%d %H:%M:%S", tz = "")
str(bci.tower$datetime)

obs.gpp.d <- bci.tower %>% select(datetime, gpp, vpd) %>% # in mumol/m2/2: units must be mumol/m2/s
  mutate(date = as.Date(format(datetime, "%Y-%m-%d")),
         gpp.mumol = as.numeric(as.character(gpp)),
         vpd = as.numeric(as.character(vpd))) %>%
  group_by(date) %>% summarise(gpp.tower = mean(gpp.mumol, na.rm = T),
                               VPD.tower = mean(vpd, na.rm = T)) %>%
  mutate(gpp.tower = gpp.tower*12*1e-06*24*60*60) %>% # gC/m2/d %>%
  subset(date != is.na(date)) %>%
  droplevels

rectangles <- data.frame(
  xmin = as.Date(paste0(c(2012:2017), "-05-01")),
  xmax = as.Date(paste0(c(2012:2017), "-12-01")),
  ymin = 0,
  ymax = 15
)

gpp.daily <- ggplot(obs.gpp.d) +
  geom_rect(data=rectangles, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill='gray80', alpha=0.8) +
  geom_point(aes(y = gpp.tower, x = date)) +
  ylab(expression('GPP (gC'*m^-2*day^-1*')')) + xlab("Date") +
  scale_x_date(breaks = c(rectangles$xmin, rectangles$xmax), date_labels = "%b-%y") +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 1, hjust = 1))
ggsave(("gpp.daily.jpeg"),
       plot = gpp.daily, file.path(figures.folder.gpp), device = "jpeg", height = 3, width = 9, units='in')

obs.qet.d <- bci.hydromet::forcings %>% select(date, AET, AET.flag.day) %>% # in mm/day
  rename(AET.gap.tower = AET.flag.day,
         AET.tower = AET) %>% # this is gap filled AET
  subset(date > "2012-07-02" & date < "2017-09-01") %>% ## date after which ET data is present
  droplevels()
head(obs.qet.d)

gpp.rel.daily <- obs.gpp.d %>%
  left_join(clim.daily %>% select(date, VPD, pet.PM, Rs), by = "date") %>%
  left_join(obs.qet.d, by = "date") %>%
  left_join(bci.hydromet::met.tower %>% select(date, Rs) %>% rename(Rs.tower = Rs), by = "date") %>%
  left_join(bci.hydromet::bci.met %>% select(date, PET_man) %>% rename(Pan.Evap = PET_man), by = "date") %>%
  left_join(bci.hydromet::met.petPM %>% select(date, pet.PM) %>% rename(pet.PM.tower = pet.PM), by = "date") %>%
  subset(date < max(obs.gpp.d$date, na.rm = TRUE)) %>%
  droplevels()

gpp.rel.daily.long <- gpp.rel.daily %>% pivot_longer(-date, values_to = "value", names_to = "variable")

gpp.rel.daily.plot <- ggplot(gpp.rel.daily.long) +
  # geom_rect(data=rectangles, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
  #           fill='gray80', alpha=0.8) +
  geom_vline(data = rectangles, aes(xintercept = xmin)) +
  geom_vline(data = rectangles, aes(xintercept = xmax)) +
  geom_line(aes(y = value, x = date, color = variable), show.legend = FALSE) +
  geom_point(size = 0.1, aes(y = value, x = date, color = variable), show.legend = FALSE) +
  facet_wrap(variable ~ . , scales = "free_y", ncol = 1) +
  ylab(expression('Value')) + xlab("Date") +
  scale_x_date(breaks = c(rectangles$xmin, rectangles$xmax), date_labels = "%b-%y") +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 1, hjust = 1))
ggsave(("gpp.re.daily.jpeg"),
       plot = gpp.rel.daily.plot, file.path(figures.folder.gpp), device = "jpeg", height = 9, width = 9, units='in')

gpp.rel.monthly <- gpp.rel.daily %>%
  mutate(month = format(date, format = "%b%Y")) %>%
  select(-date) %>%
  group_by(month) %>%
  summarise_all(list(~mean(., na.rm = TRUE)))
gpp.rel.weekly <- gpp.rel.daily %>%
  mutate(week = format(date, format = "%U%Y")) %>%
  select(-date) %>%
  group_by(week) %>%
  summarise_all(list(~mean(., na.rm = TRUE)))

data.scale.on <- "monthly"

if (data.scale.on == "monthly"){
  point.size <- 2
  gpp.rel <- gpp.rel.monthly
} else if (data.scale.on == "weekly"){
  point.size <- 2
  gpp.rel <- gpp.rel.weekly
  } else if (data.scale.on == "daily"){
  point.size <- 1
  gpp.rel <- gpp.rel.daily
}

formula = y ~ poly(x, 2, raw=TRUE)
# https://stackoverflow.com/questions/21748598/add-or-override-aes-in-the-existing-mapping-object
# add_modify_aes <- function(mapping, ...) {
#   ggplot2:::rename_aes(modifyList(mapping, ...))
# }
# ggplot(df, add_modify_aes(mapping, aes(color=new_col, y=new_y)))

gpp.vpd <- ggplot(gpp.rel, aes(y = gpp.tower, x = VPD.tower)) +
  geom_point(size = point.size, alpha = 0.7) +
  stat_smooth(method="lm", se=TRUE, fill=NA,
              formula = formula, colour = "red") +
  stat_poly_eq(aes(label = stat(eq.label)),
               npcx = 0.95, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4, colour = "red") +
  stat_poly_eq(aes(label = stat(adj.rr.label)),
               npcx = 0.95, npcy = 0.87, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4, colour = "red") +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.95, npcy = 0.77, size = 4, colour = "red") +
  ylab(expression('GPP (gC'*m^-2*day^-1*')')) + xlab("VPD (kPa)")
ggsave(paste0(data.scale.on,"_gpp.vpd.jpeg"),
       plot = gpp.vpd, file.path(figures.folder.gpp), device = "jpeg", height = 3, width = 3, units='in')
ggsave(paste0(data.scale.on,"_gpp.vpd.tiff"),
       plot = gpp.vpd, file.path(figures.folder.gpp), device = "tiff", height = 3, width = 3, units='in')

if(data.scale.on == "monthly") {
  mapping.gpp.aet <- aes(y = gpp.tower, x = AET.tower)
} else {
  mapping.gpp.aet <- aes(y = gpp.tower, x = AET.gap.tower)
}
gpp.aet <- ggplot(gpp.rel, mapping.gpp.aet) +
  geom_point(size = point.size, alpha = 0.7) +
  stat_smooth(method="lm", se=TRUE, fill=NA,
              formula=formula, colour = "red") +
  stat_poly_eq(aes(label = paste(stat(eq.label), stat(adj.rr.label), sep = "~~~~")),
               npcx = 0.95, npcy = 0.98, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4, color = "red") +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.95, npcy = 0.90, size = 4, color = "red") +
  ylab(expression('GPP (gC'*m^-2*day^-1*')')) + xlab("AET (mm)")#xlab("Penman-Monteith PET (mm)")
ggsave(paste0(data.scale.on,"_gpp.aet.jpeg"),
       plot = gpp.aet, file.path(figures.folder.gpp), device = "jpeg", height = 3, width = 3, units='in')

gpp.pet <- ggplot(gpp.rel, aes(y = gpp.tower, x = pet.PM.tower)) +
  geom_point(size = point.size, alpha = 0.7) +
  stat_smooth(method="lm", se=TRUE, fill=NA,
              formula=formula, colour = "red") +
  stat_poly_eq(aes(label = paste(stat(eq.label), stat(adj.rr.label), sep = "~~~~")),
               npcx = 0.95, npcy = 0.98, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4, color = "red") +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.95, npcy = 0.90, size = 4, color = "red") +
  ylab(expression('GPP (gC'*m^-2*day^-1*')')) + xlab("Penman-Monteith PET (mm)")
ggsave(paste0(data.scale.on,"_gpp.pet.jpeg"),
       plot = gpp.pet, file.path(figures.folder.gpp), device = "jpeg", height = 3, width = 3, units='in')

gpp.pan <- ggplot(gpp.rel, aes(y = gpp.tower, x = Pan.Evap)) +
  geom_point(size = point.size, alpha = 0.7) +
  stat_smooth(method="lm", se=TRUE, fill=NA,
              formula=formula, colour = "red") +
  stat_poly_eq(aes(label = paste(stat(eq.label), stat(adj.rr.label), sep = "~~~~")),
               npcx = 0.95, npcy = 0.98, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4, color = "red") +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.95, npcy = 0.90, size = 4, color = "red") +
  ylab(expression('GPP (gC'*m^-2*day^-1*')')) + xlab("Pan Evaporation (mm)")#xlab("Penman-Monteith PET (mm)")
ggsave(paste0(data.scale.on,"_gpp.pan.jpeg"),
       plot = gpp.pan, file.path(figures.folder.gpp), device = "jpeg", height = 3, width = 3, units='in')

gpp.rad <- ggplot(gpp.rel, aes(y = gpp.tower, x = Rs.tower)) +
  geom_point(size = point.size, alpha = 0.7) +
  stat_smooth(method="lm", se=TRUE, fill=NA,
              formula=formula, colour = "red") +
  stat_poly_eq(aes(label = paste(stat(eq.label), stat(adj.rr.label), sep = "~~~~")),
               npcx = 0.95, npcy = 0.98, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4, color = "red") +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.95, npcy = 0.90, size = 4, color = "red") +
  ylab(expression('GPP (gC'*m^-2*day^-1*')')) + xlab(expression('Solar Radiation (MJ'*m^-2*day^-1*')'))
ggsave(paste0(data.scale.on,"_gpp.rad.jpeg"),
       plot = gpp.rad, file.path(figures.folder.gpp), device = "jpeg", height = 3, width = 3, units='in')

if(data.scale.on == "monthly") {
  mapping.aet.pet <- aes(y = AET.tower, x = pet.PM.tower)
} else {
  mapping.aet.pet <- aes(y = AET.gap.tower, x = pet.PM.tower)
}
aet.pet <- ggplot(gpp.rel, mapping.aet.pet) +
  geom_abline(intercept = 0, slope = 1, lty = "dotted") +
  geom_point(size = point.size, alpha = 0.7) +
  xlim(c(2, 5)) + ylim(c(2,5)) +
  stat_smooth(method="lm", se=TRUE, fill=NA,
              formula=formula, colour = "red") +
  stat_poly_eq(aes(label = paste(stat(eq.label), stat(adj.rr.label), sep = "~~~~")),
               npcx = 0.95, npcy = 0.98, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4, color = "red") +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.95, npcy = 0.90, size = 4, color = "red") +
  ylab("AET (mm)") +xlab("Penman-Monteith PET (mm)") #xlab("Penman-Monteith PET (mm)")
ggsave(paste0(data.scale.on,"_aet.pet.jpeg"),
       plot = aet.pet, file.path(figures.folder.gpp), device = "jpeg", height = 3, width = 3, units='in')

eq.gpp.vpd.tower <- lm(gpp.tower ~ poly(VPD.tower, 2, raw = TRUE), data = gpp.rel)
eq.gpp.aet.tower <- lm(gpp.tower ~ poly(AET.gap.tower, 2, raw = TRUE), data = gpp.rel)
eq.gpp.rad.tower <- lm(gpp.tower ~ poly(Rs.tower, 2, raw = TRUE), data = gpp.rel)
eq.gpp.vpd <- lm(gpp.tower ~ poly(VPD, 2, raw = TRUE), data = gpp.rel)
eq.gpp.pet <- lm(gpp.tower ~ poly(pet.PM, 2, raw = TRUE), data = gpp.rel)
eq.gpp.rad <- lm(gpp.tower ~ poly(Rs, 2, raw = TRUE), data = gpp.rel)
eq.gpp.rad.vpd.tower.gam <- gam(gpp.tower ~ s(Rs.tower, k = 5) + s(VPD.tower, k = 5), data = gpp.rel)
eq.gpp.rad.vpd.gam <- gam(gpp.tower ~ s(Rs, k = 5) + s(VPD, k = 5), data = gpp.rel)
eq.gpp.rad.aet.gam <- gam(gpp.tower ~ s(Rs.tower, k = 5) + s(AET.gap.tower, k = 8), data = gpp.rel)
eq.gpp.rad.pan.gam <- gam(gpp.tower ~ s(Rs.tower, k = 5) + s(Pan.Evap, k = 8), data = gpp.rel)
eq.gpp.rad.pet.tower.gam <- gam(gpp.tower ~ s(Rs.tower, k = 5) + s(pet.PM.tower, k = 8), data = gpp.rel)
eq.gpp.rad.pet.gam <- gam(gpp.tower ~ s(Rs, k = 5) + s(pet.PM, k = 8), data = gpp.rel)
gam.check(eq.gpp.rad.vpd.gam)
gam.check(eq.gpp.rad.aet.gam)
gam.check(eq.gpp.rad.pan.gam)
gam.check(eq.gpp.rad.pet.gam)
summary(eq.gpp.rad.vpd.tower.gam)
summary(eq.gpp.rad.vpd.gam)
summary(eq.gpp.rad.aet.gam)
summary(eq.gpp.rad.pan.gam)
summary(eq.gpp.rad.pet.tower.gam)
summary(eq.gpp.rad.pet.gam)

# plot(eq.gpp.rad.vpd.gam$fitted.values ~ eq.gpp.rad.vpd.gam$)
# eq.gpp.vpd <- lm(gpp ~ poly(VPD, 2, raw=TRUE), data = gpp.rel.daily)
# eq.gpp.pet <- lm(gpp ~ poly(pet.PM, 2, raw=TRUE), data = gpp.rel.daily)
# eq.gpp.rad <- lm(gpp ~ poly(Rs, 2, raw=TRUE), data = gpp.rel.daily)
gpp.models <- list(eq.gpp.vpd = eq.gpp.vpd,
                   eq.gpp.pet = eq.gpp.pet,
                   eq.gpp.rad = eq.gpp.rad,
                   eq.gpp.rad.vpd.tower.gam = eq.gpp.rad.vpd.tower.gam,
                   eq.gpp.rad.pet.tower.gam = eq.gpp.rad.pet.tower.gam,
                   eq.gpp.rad.vpd.gam = eq.gpp.rad.vpd.gam,
                   eq.gpp.rad.pet.gam = eq.gpp.rad.pet.gam)
save(gpp.models, file = file.path(results.folder, "gpp.models.Rdata"))

#******************************************************
## Load BCI traits----
#******************************************************

bci.traits <- read.csv("data-raw/traits/BCITRAITS_20101220.csv") %>%
  rename(form1 = GRWFRM1., sp = SP., SG100C_AVG = SG100C_AVG) %>% mutate(sp = tolower(sp))

#******************************************************
### Load Hydraulic traits by Brett Wolfe ---------
#******************************************************
# nlrq(ks~m*(exp(-(-psych_psi/b)^a))
hyd <- read.csv("data-raw/traits/HydraulicTraits_BrettWolfe/ht1_20200103.csv") # Brett's data
# vc_a = vc_a, vc_b = vc_b, vc_a_se = vc_a_se, vc_b_se = vc_b_se,
hyd <- hyd %>% select(-genus, -species, -deciduousness, -site) %>%
  rename(KmaxS = vc_m, KmaxS_se = vc_m_se,
         LMA = lma_gm2_m, WD = xylem_den_m, TLP = tlp_m,
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
  left_join(sp.hab, by = "sp")

length(unique(hyd$sp)) # 27 sp across BCI, PNM, San Lorenzo
traits.labels.table.1 <- data.frame(trait = factor(c("lwp.min_Predawn", "lwp.min_Diurnal", "TLP",
                                                     "KmaxS", "p50S", "p88S",
                                                     "HSMTLP", "HSM50S","HSM88S",
                                                     "HSMTLP.50S", "HSMTLP.88S",
                                                     "CWR_Total", "CWR_Xylem", "CWR_Bark",
                                                     "Felbow_Xylem", "Felbow_Bark", "HSMFelbow_Xylem", "HSMFelbow_Bark",
                                                     "Fcap_Xylem", "Fcap_Bark","WD",
                                                     "Panama.moist.pref", "Plot.swp.pref", "LMA", "vc_a", "vc_b"),
                                                   levels = c("lwp.min_Predawn", "lwp.min_Diurnal", "TLP",
                                                              "KmaxS", "p50S", "p88S",
                                                              "HSMTLP", "HSM50S","HSM88S",
                                                              "HSMTLP.50S", "HSMTLP.88S",
                                                              "CWR_Total", "CWR_Xylem", "CWR_Bark",
                                                              "Felbow_Xylem", "Felbow_Bark", "HSMFelbow_Xylem", "HSMFelbow_Bark",
                                                              "Fcap_Xylem", "Fcap_Bark","WD",
                                                              "Panama.moist.pref", "Plot.swp.pref", "LMA", "vc_a", "vc_b"), ordered = TRUE)) %>%
  transform(trait.plot = factor(trait, labels = c(expression(Psi[Predawn]), expression(Psi[min]), expression(Psi[tlp]),
                                                  expression(italic('K')['max, stem']),  expression(italic('P')['50, stem']),  expression(italic('P')['88, stem']),
                                                  expression(Psi[min]*' - '*Psi[tlp]),
                                                  expression(Psi[min]*' - '*italic('P')['50, stem']),
                                                  expression(Psi[min]*' - '*italic('P')['88, stem']),
                                                  expression(Psi[tlp]*' - '*italic('P')['50, stem']),
                                                  expression(Psi[tlp]*' - '*italic('P')['88, stem']),
                                                  expression('CWR'['total']), expression('CWR'['xylem']), expression('CWR'['bark']),
                                                  expression(italic('F')['elbow, xylem']), expression(italic('F')['elbow, bark']),
                                                  expression(Psi[min]*'-'*italic('F')['Elbow, Xylem']),
                                                  expression(Psi[min]*'-'*italic('F')['Elbow, Bark']),
                                                  expression(italic('F')['cap, xylem']), expression(italic('F')['cap, bark']),
                                                  expression('WD'[stem]),
                                                  expression('Panama'[wet]), expression('Plot'[wet]), "LMA",
                                                  expression('A'['vc, stem']), expression('B'['vc, stem']))),
            trait.plot.chart = factor(trait, labels = c(expression(Psi[Predawn]), expression(Psi[min]), expression(Psi[tlp]),
                                                        expression(italic('K')['max, stem']),  expression(italic('P')['50,stem']),  expression(italic('P')['88,stem']),
                                                        expression(Psi[min]*'-'*Psi[tlp]),
                                                        expression(Psi[min]*'-'*italic('P')['50,stem']),
                                                        expression(Psi[min]*'-'*italic('P')['88,stem']),
                                                        expression(Psi[tlp]*'-'*italic('P')['50,stem']),
                                                        expression(Psi[tlp]*'-'*italic('P')['88,stem']),
                                                        expression('CWR'['total']), expression('CWR'['xylem']), expression('CWR'['bark']),
                                                        expression(italic('F')['elbow, xylem']), expression(italic('F')['elbow,bark']),
                                                        expression(Psi[min]*'-'*italic('F')['elbow,Xylem']),
                                                        expression(Psi[min]*'-'*italic('F')['elbow,Bark']),
                                                        expression(italic('F')['cap,xylem']), expression(italic('F')['cap,bark']),
                                                        expression('WD'[stem]),
                                                        expression('Panama'[wet]), expression('Plot'[wet]), "LMA",
                                                        expression('A'['vc, stem']), expression('B'['vc, stem'])))) %>% droplevels()


#******************************************************
####----Phenology by Brett Wolfe hydraulic traits-----
#******************************************************

hyd.long <- hyd %>% select(-DeciLvl) %>%
  select(sp, deciduousness, deciduous, location, TLP, p50S, p88S, vc_a, vc_b,
         CWR_Total, Fcap_Xylem, CWR_Xylem, Felbow_Xylem, Fcap_Bark, CWR_Bark, Felbow_Bark, WD, LMA,
         HSM50S, HSM88S, HSMTLP, HSMFelbow_Xylem, HSMFelbow_Bark, HSMTLP.50S, HSMTLP.88S,
         Panama.moist.pref, Plot.swp.pref, lwp.min_Diurnal, lwp.min_Predawn) %>%
  gather(trait, value, -sp, -deciduousness, -deciduous,  -location) %>%
  unite("deci_sp", deciduous, sp, remove = FALSE) %>%
  subset(deciduousness != "NA") %>%
  droplevels()
kruskal.list <- list()
for(i in unique(hyd.long$trait)) {
  xx <- hyd.long %>% subset(trait == i)
  kruskal.list[[i]] <- cbind(trait = i, kruskal(xx$value, xx$deciduousness, alpha = 0.1, group = TRUE, p.adj = "bonferroni")$groups,
                             deciduousness = rownames(kruskal(xx$value, xx$deciduousness, alpha = 0.1, group = TRUE, p.adj = "bonferroni")$groups))
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

hyd.error <- hyd %>% select(sp, KmaxS_se, vc_b_se, vc_a_se, tlp_sd) %>%
  rename(KmaxS = KmaxS_se, vc_b = vc_b_se, vc_a = vc_a_se, TLP = tlp_sd) %>%
  gather(trait, se, -sp) %>%
  ## But for TLP it's not se but sd
  mutate(se = ifelse(trait == "TLP", NA, se),
         sd = ifelse(trait == "TLP", se, NA))
hyd.long <- hyd.long %>% left_join(hyd.error, by = c("sp", "trait"))
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
  left_join(bci.traits %>% select(sp, form1, SG100C_AVG), by = "sp")

traits.labels.table.2 <- data.frame(trait = factor(c("KmaxL", "lwp.min_Predawn", "lwp.min_Diurnal", "TLP", "p50L", "p80L",
                                                     "HSMLWP.TLP", "HSMLWP.50L", "HSMTLP.50L",
                                                     "HSMLWP.80L", "HSMTLP.80L",
                                                     "Panama.moist.pref", "Plot.swp.pref", "SG100C_AVG", "Chl"),
                                                   levels = c("depth", "Xylem_sap_deltaD_permil",
                                                              "KmaxL", "lwp.min_Predawn", "lwp.min_Diurnal", "TLP", "p50L", "p80L",
                                                              "HSMLWP.TLP", "HSMLWP.50L", "HSMTLP.50L",
                                                              "HSMLWP.80L", "HSMTLP.80L",
                                                              "Panama.moist.pref", "Plot.swp.pref", "SG100C_AVG", "Chl"), ordered = TRUE)) %>%
  transform(trait.plot = factor(trait, labels = c(expression(italic(K)[max]), expression(Psi[predawn]), expression(Psi[min]),
                                                  expression(Psi[tlp]), expression(italic('P')['50, leaf']), expression(italic('P')['80, leaf']),
                                                  expression(Psi[min]*' - '*Psi[tlp]),
                                                  expression(Psi[min]*' - '*italic('P')['50, leaf']),
                                                  expression(Psi[tlp]*' - '*italic('P')['50, leaf']),
                                                  expression(Psi[min]*' - '*italic('P')['80, leaf']),
                                                  expression(Psi[tlp]*' - '*italic('P')['80, leaf']),
                                                  expression('Panama'[wet]), expression('Plot'[wet]), expression('SG'[100*~degree*C]), "LMA")),
            trait.plot.chart = factor(trait, labels = c(expression(italic(K)[max]), expression(Psi[predawn]), expression(Psi[min]),
                                                        expression(Psi[tlp]), expression(italic('P')['50,leaf']), expression(italic('P')['80,leaf']),
                                                        expression(Psi[min]*'-'*Psi[tlp]),
                                                        expression(Psi[min]*'-'*italic('P')['50,leaf']),
                                                        expression(Psi[tlp]*'-'*italic('P')['50,leaf']),
                                                        expression(Psi[min]*'-'*italic('P')['80,leaf']),
                                                        expression(Psi[tlp]*'-'*italic('P')['80,leaf']),
                                                        expression('Panama'[wet]), expression('Plot'[wet]), expression('SG'[100*~degree*C]), "LMA")))

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
  full_join(bci.traits %>% select(sp, form1, SG100C_AVG), by = "sp") %>%
  full_join(deci %>% select(-sp4), by = "sp") %>%
  subset(sp %in% unique(c(hyd$sp, traits$sp))) %>%
  left_join(traits.wide %>% select(-form1, -deciduous, -deciduousness,
                                   -SG100C_AVG, -Panama.moist.pref, -Plot.swp.pref), by = "sp") %>%
  pivot_longer(cols = c(-sp, -form1, -deciduous, -deciduousness, -DeciLvl,
                        -deciduousness.label, -sp.plot, -deci_sp, -deci_sp.plot),
               names_to = "trait", values_to = "value") %>%
  unite("deci_sp", deciduous, sp, remove = FALSE) %>%
  left_join(traits.labels.table.2, by = "trait") %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(deciduousness)]), ordered=TRUE),
         deci_sp.plot = factor(deci_sp, levels=unique(deci_sp[order(deciduousness)]), ordered=TRUE)) %>%
  droplevels()

save(traits, file = file.path(results.folder, "kunert.traits.all.RData"))
save(traits.long, file = file.path(results.folder, "kunert.traits.key.long.RData"))
save(traits.long.hyd, file = file.path(results.folder, "kunert.traits.key.long_in_Wolfe_traits_species_list.RData"))


#******************************************************
## Predictors of KmaxL-LWP relationship --------
#******************************************************
# WD is measured by Nobby at 60 deg C
traits.for.kcurves <- traits.indi %>%
  rename(TLP = mean_TLP_Mpa, Chl = Chl_m2_per_g, LMA = LMA_g_per_m2,
         SPAD = mean_SPAD, WD = WD_g_per_cm3, LD = LD_g_per_cm3) %>%
  select(sp, TLP, Chl, LMA, SPAD, WD, LD) %>% group_by(sp) %>%
  summarise_all(mean, na.rm = TRUE)
sp.exp.param <- leaf_cond.models %>% subset(model_type == "Exponential") %>%
  mutate(sp = data.type) %>% select(sp, A, B, Kmax, Kmax_at_0.1) %>%
  left_join(traits.for.kcurves, by = "sp") %>%
  left_join(deci %>% select(-sp4, -deciduousness.label), by = "sp") %>%
  left_join(bci.traits %>%
              select(sp, SG100C_AVG, SG60C_AVG, LMALEAF_AVD, LMALAM_AVD, LMADISC_AVD, ## LMADISC_AVD is not available when LMALAM_AVD isn't
                     HEIGHT_AVG, WSG_CHAVE), by = "sp")
range(sp.exp.param$Kmax, na.rm = TRUE)
# 0.9679718 9.9476435
range(sp.exp.param$Kmax_at_0.1, na.rm = TRUE)
# 0.8453948 7.1963332

## available data
nrow(sp.exp.param %>% subset(!is.na(WD))) # 21 spp
nrow(sp.exp.param %>% subset(!is.na(SG100C_AVG))) # 19 spp
nrow(sp.exp.param %>% subset(!is.na(SG60C_AVG))) # 19 spp
nrow(sp.exp.param %>% subset(!is.na(WSG_CHAVE))) # 20 spp
nrow(sp.exp.param %>% subset(!is.na(LMA))) # 21 spp
nrow(sp.exp.param %>% subset(!is.na(LMALEAF_AVD))) # 16 spp
nrow(sp.exp.param %>% subset(!is.na(LMADISC_AVD))) # 19 spp
nrow(sp.exp.param %>% subset(!is.na(LMALAM_AVD))) # 16 spp

## available data in bci.traits
nrow(bci.traits %>% subset(!is.na(SG100C_AVG))) # 517 spp
nrow(bci.traits %>% subset(!is.na(SG60C_AVG))) # 478 spp
nrow(bci.traits %>% subset(!is.na(WSG_CHAVE))) # 152 spp
nrow(bci.traits %>% subset(!is.na(LMALEAF_AVD))) # 253 spp
nrow(bci.traits %>% subset(!is.na(LMADISC_AVD))) # 214 spp
nrow(bci.traits %>% subset(!is.na(LMALAM_AVD))) # 203 spp
singlepoly <- function(vec.x, vec.y, common.spp = TRUE, given.df) {
  result.list <- list()
  for(i in 1 : length(vec.x)) {
    rn <- i
    if (common.spp == TRUE) {
      ## to compare with data present on all variables on same species
      df <- given.df[complete.cases(given.df), ]
    } else {
      df <- given.df[!is.na(given.df[, vec.y[i]]),]
    }
    y <- df[, vec.y[i]]
    x <- df[, vec.x[i]]
    model1 <- lm(y ~ poly(x, degree = 2, raw = TRUE), data = df)
    ## rounded coefficients for better output
    cf <- round(coef(model1), 2)
    ## sign check to avoid having plus followed by minus for negative coefficients
    eq <- paste0("y = ", cf[1],
                 ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), "x",
                 ifelse(sign(cf[3])==1, " + ", " - "), abs(cf[3]), "x2")
    r2 <- paste0("R-squared = ", round(summary(model1)$r.squared, 2),
                 ", P-val = ", round(lmp(model1), 2))
    result.df <- data.frame(y =  vec.y[i], x =  vec.x[i],
                            eq = eq, R2 = round(summary(model1)$r.squared, 2),
                            P = round(lmp(model1), 2), n = nrow(df))
    result.list[[rn]] <- result.df
  }
  return(do.call("rbind", result.list))
}

comm.singlevar.1 <- singlepoly(given.df = bci.traits, vec.x = c("SG60C_AVG", "SG100C_AVG", "WSG_CHAVE", "WSG_CHAVE"),
                               vec.y = c("SG100C_AVG", "SG60C_AVG","SG100C_AVG", "SG60C_AVG"), common.spp = FALSE) %>% arrange(desc(R2))
comm.singlevar.2 <- singlepoly(given.df = bci.traits, vec.x = c("LMALAM_AVD", "LMADISC_AVD",  "LMALEAF_AVD", "LMADISC_AVD"),
                               vec.y = c("LMALEAF_AVD", "LMALEAF_AVD","LMALAM_AVD", "LMALAM_AVD"), common.spp = FALSE) %>% arrange(desc(R2))

comm.singlevar <- rbind(comm.singlevar.1, comm.singlevar.2)
save(comm.singlevar, file = file.path(results.folder, "comm.singlevar_common.spp_FALSE.Rdata"))
# y           x                         eq   R2 P   n
# 1  SG100C_AVG   SG60C_AVG  y = 0.02 + 0.92x + 0.04x2 0.99 0 517
# 2   SG60C_AVG  SG100C_AVG y = -0.01 + 1.06x - 0.04x2 0.99 0 478
# 3   WSG_CHAVE   SG60C_AVG  y = 0.02 + 0.91x + 0.03x2 0.77 0 152
# 4   WSG_CHAVE  SG100C_AVG  y = 0.01 + 0.96x + 0.01x2 0.76 0 152
# 5 LMALEAF_AVD  LMALAM_AVD    y = 13.05 + 0.84x - 0x2 0.94 0 253
# 6  LMALAM_AVD LMALEAF_AVD    y = -2.89 + 0.99x - 0x2 0.93 0 203
# 7 LMADISC_AVD  LMALAM_AVD    y = -8.25 + 1.28x - 0x2 0.90 0 214
# 8 LMADISC_AVD LMALEAF_AVD   y = -12.15 + 1.26x - 0x2 0.83 0 214

## Do Nobby's measures predict larger common measures?

wd.lma.singlevar.1 <- singlepoly(given.df = sp.exp.param, vec.x = rep("WD", 3), vec.y = c("SG100C_AVG", "SG60C_AVG", "WSG_CHAVE")) %>% arrange(desc(R2))
wd.lma.singlevar.2 <- singlepoly(given.df = sp.exp.param, vec.x = rep("LMA", 3), vec.y = c("LMALEAF_AVD", "LMADISC_AVD", "LMALAM_AVD")) %>% arrange(desc(R2))
wd.lma.singlevar <- rbind(wd.lma.singlevar.1, wd.lma.singlevar.2)
save(wd.lma.singlevar, file = file.path(results.folder, "wd.lma.singlevar_common.spp_TRUE.Rdata"))
# y   x                          eq   R2    P  n
# 1   SG60C_AVG  WD  y = -0.25 + 2.52x - 1.64x2 0.58 0.01 14
# 2  SG100C_AVG  WD  y = -0.14 + 2.06x - 1.23x2 0.57 0.01 14
# 3   WSG_CHAVE  WD  y = -0.17 + 2.11x - 1.39x2 0.48 0.03 14
# 4 LMADISC_AVD LMA    y = 88.4 - 0.8x + 0.01x2 0.93 0.00 14
# 5  LMALAM_AVD LMA y = 112.19 - 1.43x + 0.01x2 0.86 0.00 14
# 6 LMALEAF_AVD LMA y = 132.26 - 1.72x + 0.02x2 0.79 0.00 14
# wd.lma.singlevar
wd.lma.singlevar.1 <- singlepoly(given.df = sp.exp.param, vec.x = rep("WD", 3), vec.y = c("SG100C_AVG", "SG60C_AVG", "WSG_CHAVE"), common.spp = FALSE) %>% arrange(desc(R2))
wd.lma.singlevar.2 <- singlepoly(given.df = sp.exp.param, vec.x = rep("LMA", 3), vec.y = c("LMALEAF_AVD", "LMADISC_AVD", "LMALAM_AVD"), common.spp = FALSE) %>% arrange(desc(R2))
wd.lma.singlevar <- rbind(wd.lma.singlevar.1, wd.lma.singlevar.2)
save(wd.lma.singlevar, file = file.path(results.folder, "wd.lma.singlevar_common.spp_FALSE.Rdata"))

## With common species, Nobby's LMA is well related LMADISC followed by LMALAM and less so with LMALEAF
## but WD is only moderately related to SG100C_AVG and SG60C_AVG
# y   x                          eq             R2    P  n
# 1  SG100C_AVG  WD    y = 0.04 + 1.4x - 0.71x2 0.45 0.01 19
# 2   SG60C_AVG  WD  y = -0.04 + 1.74x - 1.02x2 0.45 0.01 19
# 3   WSG_CHAVE  WD  y = -0.04 + 1.63x - 0.99x2 0.41 0.01 20
# 4 LMALEAF_AVD LMA y = 137.51 - 1.86x + 0.02x2 0.78 0.00 16
# 5 LMADISC_AVD LMA  y = 90.27 - 0.56x + 0.01x2 0.58 0.00 19
# 6  LMALAM_AVD LMA y = 117.58 - 1.57x + 0.01x2 0.86 0.00 16

## Nobby's LMA is well related to LMALAM and followed by LMALEAF and less so with LMADISC
## but WD is only moderately related to SG100C_AVG and SG60C_AVG


figures.folder.kleaf <- paste0("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf")
if(!dir.exists(file.path(figures.folder.kleaf))) {dir.create(file.path(figures.folder.kleaf))}

## ******
### PCA to check which variables strongly correlate with Parameters A & B
## ******
df.pca <- sp.exp.param %>% remove_rownames %>% column_to_rownames(var = "sp") %>% select_if(is.numeric) ## more species without Kmax data %>% select(-SafetyMargin.p50, -p50_Kmax, -Kmax)
# if (diff(range(df.pca$root.95, na.rm = TRUE)) == 0) {
df.pca <- df.pca %>% subset(complete.cases(df.pca))
result.pca <- prcomp(df.pca, center = TRUE, scale = TRUE)

pdf(file.path(figures.folder.kleaf ,"kamx_by_psi_params_pca.pdf"), height = 6, width = 6)
biplot(result.pca, choices = 1:2, pc.biplot = TRUE, main = "")
dev.off()

## single explanatory variable
var.y <- c("A", "A", "A", "A", "A", "A", "A", "A",
           "B", "B", "B", "B", "B", "B", "B")
var.x <- c("SG100C_AVG", "SG60C_AVG", "WSG_CHAVE", "WD", "LMA", "LMALAM_AVD", "LMALEAF_AVD", "B",
           "SG100C_AVG", "SG60C_AVG",  "WSG_CHAVE", "WD", "LMA", "LMALAM_AVD", "LMALEAF_AVD")
ab.singlevar <- data.frame(y = var.y, x = var.x, eq = NA, R2 = NA, P = NA, n = NA)

for(i in 1:length(var.x)){
  jpeg(file.path(figures.folder.kleaf, paste0("kmax_by_psi_", var.y[i], "_by_", var.x[i], ".jpeg")),
       width = 4, height = 4, units = "in", pointsize = 10,
       quality = 100, res = 300)
  df <- sp.exp.param[!is.na(sp.exp.param[, var.x[i]]),]
  df$x <- df[, var.x[i]]
  df$y <- df[, var.y[i]]
  model1 <- lm(y ~ poly(x, 2), data = df, raw = TRUE)
  plot(y ~ x, data = df, xlab = var.x[i], ylab = var.y[i])
  x0 <- seq(min(df$x, na.rm = TRUE), max(df$x, na.rm = TRUE), length = 100)  ## prediction grid
  y0 <- predict(model1, newdata = list(x = x0))  ## predicted values
  lines(x0, y0, col = 2)
  ## rounded coefficients for better output
  cf <- round(coef(model1), 2)
  ## sign check to avoid having plus followed by minus for negative coefficients
  eq <- paste0("y = ", cf[1],
               ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), "x",
               ifelse(sign(cf[3])==1, " + ", " - "), abs(cf[3]), "x2")
  r2 <- paste0("R-squared = ", round(summary(model1)$r.squared, 2),
               ", P-val = ", round(lmp(model1), 2))
  mtext(eq, 3, line = -2)
  mtext(r2, 3, line = -4)
  dev.off()
  ab.singlevar$eq[i] <- eq
  ab.singlevar$R2[i] <- round(summary(model1)$r.squared, 2)
  ab.singlevar$P[i] <- round(lmp(model1), 2)
  ab.singlevar$n[i]<- nrow(sp.exp.param[!is.na(sp.exp.param[, var.x[i]]),])
}
# ab.singlevar
# y           x                        eq   R2    P  n
# 1  A  SG100C_AVG y = 3.04 + 5.31x + 1.64x2 0.30 0.05 19
# 2  A   SG60C_AVG y = 3.04 + 5.08x + 1.73x2 0.28 0.07 19
# 3  A   WSG_CHAVE y = 3.06 + 5.03x + 1.43x2 0.27 0.07 20
# 4  A          WD y = 3.38 + 3.36x - 0.67x2 0.08 0.47 21
# 5  A         LMA y = 3.38 + 4.78x - 2.42x2 0.20 0.14 21
# 6  A  LMALAM_AVD y = 3.75 - 0.19x - 4.08x2 0.14 0.38 16
# 7  A LMALEAF_AVD y = 3.75 - 0.45x - 3.63x2 0.11 0.47 16
# 8  A           B y = 3.38 + 8.74x - 2.07x2 0.55 0.00 21
# 9  B  SG100C_AVG    y = 2.13 + 2.84x + 2x2 0.60 0.00 19
# 10 B   SG60C_AVG y = 2.13 + 2.73x + 2.09x2 0.58 0.00 19
# 11 B   WSG_CHAVE y = 2.25 + 2.37x + 2.07x2 0.39 0.02 20
# 12 B          WD y = 2.36 + 1.41x - 0.11x2 0.07 0.55 21
# 13 B         LMA y = 2.36 + 1.34x + 0.52x2 0.07 0.53 21
# 14 B  LMALAM_AVD y = 2.49 + 0.96x - 0.84x2 0.06 0.67 16
# 15 B LMALEAF_AVD y = 2.49 + 1.18x - 1.34x2 0.12 0.45 16

# So SG100C_AVG, SG60C_AVG,  WSG_CHAVE, are fairly ok predictors of A, but not WD.
# LMA is a weak predictor of A. LMALAM_AVD LMALEAF_AVD are worse
# B is fairly well predictor of A

# SG100C_AVG, SG60C_AVG are good predictors of B, WSG_CHAVE is OK, but not WD
## LMA and LMALAM_AVD are not, LMALEAF_AVD is weakly related

## Best R2 have two variables

multipoly <- function(vec.z, vec.x, vec.y, common.spp = TRUE, given.df) {
  result.list <- list()
    for (i in 1 : length(vec.y)) {
      for(j in 1 : length(vec.x)) {
        rn <- (i - 1)*length(vec.x) + j
        if (common.spp == TRUE) {
          ## to compare with data present on all variables on same species
          df.select <- given.df %>% select(vec.z, vec.x, vec.y)
          df <- df.select[complete.cases(df.select), ]
        } else {
          df <- given.df[!is.na(given.df[, vec.y[i]]) & !is.na(given.df[, vec.x[j]]), ]
        }
        z <- df[, vec.z]
        y <- df[, vec.y[i]]
        x <- df[, vec.x[j]]
        model1 <- lm(z ~ polym(x, y, degree = 2, raw = TRUE), data = df)
        ## rounded coefficients for better output
        cf <- round(coef(model1), 2)
        ## sign check to avoid having plus followed by minus for negative coefficients
        eq <- paste0("z = ", cf[1],
                     ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), "x",
                     ifelse(sign(cf[3])==1, " + ", " - "), abs(cf[3]), "x2",
                     ifelse(sign(cf[4])==1, " + ", " - "), abs(cf[4]), "xy",
                     ifelse(sign(cf[5])==1, " + ", " - "), abs(cf[5]), "y",
                     ifelse(sign(cf[6])==1, " + ", " - "), abs(cf[6]), "y2")
        r2 <- paste0("R-squared = ", round(summary(model1)$r.squared, 2),
                     ", P-val = ", round(lmp(model1), 2))
        result.df <- data.frame(z = vec.z, x =  vec.x[j], y =  vec.y[i],
                                eq = eq, R2 = round(summary(model1)$r.squared, 2),
                                P = round(lmp(model1), 2), n = nrow(df))
        result.list[[rn]] <- result.df
      }
    }
  return(do.call("rbind", result.list))
}
var.z <- c("B")
var.x <- c("SG100C_AVG", "SG60C_AVG", "WSG_CHAVE", "WD") #, "SG100C_AVG", "SG60C_AVG", "WSG_CHAVE", "WD")
var.y <- c("LMA", "LMALAM_AVD", "LMALEAF_AVD", "LMADISC_AVD")  #"SG100C_AVG", "SG60C_AVG", "WSG_CHAVE", "WD", "LMA", "LMALAM_AVD", "LMALEAF_AVD")

B.multivar <- multipoly(var.z, var.x, var.y, common.spp = TRUE, given.df = sp.exp.param) %>%
  arrange(desc(R2))
save(B.multivar, file = file.path(results.folder, "B.multivar_common.spp_TRUE.Rdata"))
B.multivar <- multipoly(var.z, var.x, var.y, common.spp = FALSE, given.df = sp.exp.param) %>%
  arrange(desc(R2))
save(B.multivar, file = file.path(results.folder, "B.multivar_common.spp_FALSE.Rdata"))
# When species are common (14), SG100C_AVG & SG60C_AVG with any LMAs are the best and nearly equal (R2 0.88 - 0.85)
# Worst fits with WD

A.multivar.1 <- multipoly("A", c("B", "SG100C_AVG", "SG60C_AVG", "WSG_CHAVE", "WD"),
                          c("LMA", "LMALAM_AVD", "LMALEAF_AVD", "LMADISC_AVD"),
                          common.spp = TRUE, given.df = sp.exp.param)

A.multivar.2 <- multipoly("A", c("B"), c("SG100C_AVG", "SG60C_AVG", "WSG_CHAVE", "WD"),
                          common.spp = TRUE, given.df = sp.exp.param)
A.multivar <- rbind(A.multivar.1, A.multivar.2) %>% arrange(desc(R2))
save(A.multivar, file = file.path(results.folder, "A.multivar_common.spp_TRUE.Rdata"))

A.multivar.1 <- multipoly("A", c("B", "SG100C_AVG", "SG60C_AVG", "WSG_CHAVE", "WD"),
                          c("LMA", "LMALAM_AVD", "LMALEAF_AVD", "LMADISC_AVD") ,
                          common.spp = FALSE, given.df = sp.exp.param)

A.multivar.2 <- multipoly("A", c("B"), c("SG100C_AVG", "SG60C_AVG", "WSG_CHAVE", "WD"),
                          common.spp = FALSE, given.df = sp.exp.param)
A.multivar <- rbind(A.multivar.1, A.multivar.2) %>% arrange(desc(R2))
save(A.multivar, file = file.path(results.folder, "A.multivar_common.spp_FALSE.Rdata"))

## When species are common (14), B with SG60C_AVG followed by SG100C_AVG or WD
## But all the data present seen B with LMA (21spp) R2 = 0.75; B with SG60C_AVG (19 spp) R2 = 0.65
## B with LMALAM (16spp) R2 = 0.78; and LMADISC (19 spp) R2 = 0.75
## LMA predicts LMALAM well, so using B with LMALAM

## gap-filling LMALAM with LMADISC or LMA or, LMALEAF
## gap-filling WSG100 with WSG60 or WSG_CHAVE, finally, and *may be* WD becuase it doesn't predict well
gap.models <- vector(mode = "list", length = 6)
names(gap.models) <- c("LMA.LAM.LEAF", "LMA.LAM.DISC", "LMA.LAM.LMA",
                            "WSG.100.60", "WSG.100.CHAVE", "WSG.100.WD")
gap.models$LMA.LAM.LEAF <- lm(LMALAM_AVD ~ LMALEAF_AVD, data = bci.traits)
gap.models$LMA.LAM.DISC <- lm(LMALAM_AVD ~ LMADISC_AVD, data = bci.traits)
gap.models$LMA.LAM.LMA <- lm(LMALAM_AVD ~ LMA, data = sp.exp.param)
gap.models$WSG.100.60 <- lm(SG100C_AVG ~ SG60C_AVG, data = bci.traits)
gap.models$WSG.100.CHAVE <- lm(SG100C_AVG ~ WSG_CHAVE, data = bci.traits)
gap.models$WSG.100.WD <- lm(SG100C_AVG ~ WD, data = sp.exp.param)

sp.exp.param$LMALAM_AVD.filled <- sp.exp.param$LMALAM_AVD;
sp.exp.param$SG100C_AVG.filled <- sp.exp.param$SG100C_AVG
## Out of five gaps in LMALAM_AVD three gaps filled with LMA.LAM.DISC, two by LMA, LMALEAF has same gaps, so can't be filled
## Two gaps in	SG100C_AVG, but SG60C_AVG has the same gaps. WSG_CHAVE may have filled one, WD two, but they are not strongly related with WSG100C and worsen the fit

sp.exp.param <- sp.exp.param %>%
  mutate(LMALAM_AVD.filled =
           ifelse(is.na(LMALAM_AVD.filled), predict(gap.models$LMA.LAM.DISC, newdata = sp.exp.param), LMALAM_AVD.filled),
         LMALAM_AVD.filled =
           ifelse(is.na(LMALAM_AVD.filled), predict(gap.models$LMA.LAM.LMA, newdata = sp.exp.param), LMALAM_AVD.filled))

var.z <- c("B")
var.x <- c("SG100C_AVG.filled", "WD") #"SG100C_AVG", "SG60C_AVG", "WSG_CHAVE", "WD")
var.y <- c("LMA", "LMALAM_AVD.filled")  #"SG100C_AVG", "SG60C_AVG", "WSG_CHAVE", "WD", "LMA", "LMALAM_AVD", "LMALEAF_AVD")

B.multivar <- multipoly(var.z, var.x, var.y, common.spp = TRUE, given.df = sp.exp.param) %>%
  arrange(desc(R2))
B.multivar
save(B.multivar, file = file.path(results.folder, "B.multivar_filled_common.spp_TRUE.Rdata"))
A.multivar.1 <- multipoly("A", c("B"),
                          c("LMALAM_AVD.filled", "LMA"), common.spp = TRUE, given.df = sp.exp.param)
                          # ## matching species with those in B slightly reduces R2
                          # given.df = sp.exp.param[complete.cases(sp.exp.param %>%
                          #                                          select(c("SG100C_AVG.filled", "WD", "LMALAM_AVD.filled", "LMA"))),])
A.multivar.1
save(A.multivar, file = file.path(results.folder, "A.multivar_filled_common.spp_TRUE.Rdata"))

## Gap filling in the community-level data
sp.soft.filled <- list()
## But LMALMA_AVD are fewer than LMALEAF_AVD, so filling up some gaps in them
sp.soft.filled$LMALEAF_AVD <- length(bci.traits$LMALAM_AVD[is.na(bci.traits$LMALAM_AVD) & !is.na(bci.traits$LMALEAF_AVD)])
# 50 sp can be filled
summary(gap.models$LMA.LAM.LEAF)
sp.soft.filled$LMADISC_AVD <- length(bci.traits$LMALAM_AVD[is.na(bci.traits$LMALAM_AVD) &
                                                             is.na(bci.traits$LMALEAF_AVD) &
                                                             !is.na(bci.traits$LMADISC_AVD)])
# Only 9 sp can be filled
# Adjusted R-squared:  0.8703 p-value: < 2.2e-16
sp.soft.filled$SG60C_AVG <- length(bci.traits$SG100C_AVG[is.na(bci.traits$SG100C_AVG) & !is.na(bci.traits$SG60C_AVG)])
## So only 1 sp can be substituted by SG60, it is not filled by WSG_CHAVE, so filling up


k_by_psi.models <- vector(mode = "list", length = 7)
names(k_by_psi.models) <- c("B.WSG100.LMA", "A.B.LMA.LAM",
                            "B.WSG100", "B.LMA.LAM",
                            "A.WSG100", "A.LMA.LAM", "A.B")
#models to predict from community level data
### B
k_by_psi.models$B.WSG100.LMA <- lm(B ~ polym(SG100C_AVG, LMALAM_AVD.filled, degree = 2, raw = TRUE),
                                   data = sp.exp.param)
summary(k_by_psi.models$B.WSG100.LMA)
## With gap-filled data
# Multiple R-squared:  0.7778,	Adjusted R-squared:  0.6924
# F-statistic: 9.104 on 5 and 13 DF,  p-value: 0.0006666
## with original data
# nrow(sp.exp.param %>% subset(!is.na(B) & !is.na(LMALAM_AVD) & !is.na(SG100C_AVG)))
# 14 spp
# Adjusted R-squared:  0.8122;  p-value: 0.001392
## with only SG100C_AVG, R-squared:  0.547; p-value: 0.0006908
nrow(sp.exp.param %>% subset(!is.na(B) & !is.na(LMALAM_AVD) & !is.na(SG100C_AVG)))
k_by_psi.models$B.WSG100 <- lm(B ~ poly(SG100C_AVG, degree = 2, raw = TRUE),
                                   data = sp.exp.param)
k_by_psi.models$B.LMA.LAM <- lm(B ~ poly(LMALAM_AVD.filled, degree = 2, raw = TRUE),
                                   data = sp.exp.param)

k_by_psi.models$A.B.LMA.LAM <- lm(A ~ polym(B, LMALAM_AVD.filled, degree = 2, raw = TRUE), data = sp.exp.param)
summary(k_by_psi.models$A.B.LMA.LAM)
## Using gap-filled data: 21 spp
## Adjusted R-squared:  0.7371; p-value: 7.509e-05
## Using original data
nrow(sp.exp.param %>% subset(!is.na(B) & !is.na(LMALAM_AVD)))
# 16
# Adjusted R-squared:  0.667, p-value: 0.004658
## so using either LMALAM_AVD is better than LMALEAF_AVD in predicting A
k_by_psi.models$A.WSG100 <- lm(A ~ poly(SG100C_AVG, degree = 2, raw = TRUE),
                               data = sp.exp.param)
k_by_psi.models$A.LMA.LAM <- lm(A ~ poly(LMALAM_AVD.filled, degree = 2, raw = TRUE),
                                data = sp.exp.param)
k_by_psi.models$A.B <- lm(A ~ poly(B, degree = 2, raw = TRUE),
                                data = sp.exp.param)
save(k_by_psi.models, file = file.path(results.folder, "k_by_psi.models.Rdata"))

###**********
## Predict A & B from soft traits and the above models
###**********
bci.AB <- bci.traits %>% select(sp, SG100C_AVG, SG60C_AVG, WSG_CHAVE,
                                LMALAM_AVD, LMALEAF_AVD, LMADISC_AVD) %>%
  left_join(sp.exp.param %>% select(sp, LMA), by = "sp")

bci.AB <- bci.AB %>% mutate(sp.LMA.sub = ifelse(!is.na(LMALAM_AVD), "original", NA),
                            sp.WSG.sub = ifelse(!is.na(SG100C_AVG), "original", NA),
                            SG100C_AVG = ifelse(is.na(SG100C_AVG),
                                     predict(gap.models$WSG.100.60, newdata = bci.AB),
                                     SG100C_AVG),
                            sp.WSG.sub = ifelse(!is.na(SG100C_AVG) & is.na(sp.WSG.sub),
                                                "filled-WSG60", sp.WSG.sub),
                            SG100C_AVG = ifelse(is.na(SG100C_AVG),
                                                predict(gap.models$WSG.100.CHAVE, newdata = bci.AB),
                                                SG100C_AVG),
                            sp.WSG.sub = ifelse(!is.na(SG100C_AVG) & is.na(sp.WSG.sub),
                                                "filled-CHAVE", sp.WSG.sub),
                            LMALAM_AVD =
                              ifelse(is.na(LMALAM_AVD),
                                     predict(gap.models$LMA.LAM.DISC, newdata = bci.AB),
                                     LMALAM_AVD),
                            sp.LMA.sub = ifelse(!is.na(LMALAM_AVD) & is.na(sp.LMA.sub),
                                                "filled-LMA-DISC", sp.LMA.sub),
                            LMALAM_AVD =
                              ifelse(is.na(LMALAM_AVD),
                                     predict(gap.models$LMA.LAM.LMA, newdata = bci.AB),
                                     LMALAM_AVD),
                            sp.LMA.sub = ifelse(!is.na(LMALAM_AVD) & is.na(sp.LMA.sub),
                                                "filled-LMA-Nobby", sp.LMA.sub),
                            LMALAM_AVD = ifelse(is.na(LMALAM_AVD),
                                     predict(gap.models$LMA.LAM.LEAF, newdata = bci.AB),
                                     LMALAM_AVD),
                            sp.LMA.sub = ifelse(!is.na(LMALAM_AVD) & is.na(sp.LMA.sub),
                                                "filled-LEAF", sp.LMA.sub)) %>%
  mutate(LMALAM_AVD.filled = LMALAM_AVD, SG100C_AVG.filled = SG100C_AVG)

bci.AB <- bci.AB %>% mutate(B = predict(k_by_psi.models$B.WSG100.LMA, newdata = bci.AB))
bci.AB <- bci.AB %>% mutate(A = predict(k_by_psi.models$A.B.LMA.LAM, newdata = bci.AB))
sp.soft.filled$SG100C_AVG.gaps.left <- length(bci.AB$SG100C_AVG[is.na(bci.AB$SG100C_AVG)])
## 95 species
sp.soft.filled$LMALEAF_AVD.gaps.left <- length(bci.AB$LMALAM_AVD[is.na(bci.AB$LMALAM_AVD)])
## 359 left
sp.exp.param <- sp.exp.param %>% mutate(model.B = predict(k_by_psi.models$B.WSG100.LMA, newdata = sp.exp.param))
sp.exp.param <- sp.exp.param %>% mutate(model.A = predict(k_by_psi.models$A.B.LMA.LAM, newdata = sp.exp.param))

###**********
## Do params A & B predicted from soft traits match those fitted to data?
###**********
data.model.AB <- bci.AB %>%
  subset(B >= 0 & A >= 0) %>%
  rename(model.A = A, model.B = B) %>%
  left_join(sp.exp.param %>% select(sp, A, B) %>%
              rename(data.A = A, data.B = B), by = "sp") %>%
  droplevels()
# species with AB data and modeled A B >= 0
sp.soft.filled$AB.data.model <- length(unique(data.model.AB$sp[!is.na(data.model.AB$data.A) &
                                 !is.na(data.model.AB$model.A) & !is.na(data.model.AB$model.B)]))
comm.predict.list <- vector("list", nrow(data.model.AB))
names(comm.predict.list) <- data.model.AB$sp

for (i in 1:nrow(data.model.AB)) {
  params <- data.model.AB[i,]
  df <- data.frame(psi = seq(0, 5, length.out = 1000)) %>%
    mutate(sp = params$sp,
           k.predict = Exponential(A = params$model.A, B = params$model.B, psi = psi))
  comm.predict.list[[i]] <- df

  Kmax <- max(df$k.predict, na.rm = TRUE)
  Kmax_at_0.1 <- df$k.predict[which(round(df$psi, 3) == 0.1)]

  data.model.AB$psi_kl20[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.8*Kmax)$y
  data.model.AB$psi_kl50[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.5*Kmax)$y
  data.model.AB$psi_kl80[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.2*Kmax)$y
  data.model.AB$psi_kl88[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.12*Kmax)$y
  data.model.AB$psi_kl95[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.05*Kmax)$y

  data.model.AB$Kmax[i] <- Kmax
  data.model.AB$Kmax_at_0.1[i] <- Kmax_at_0.1
  data.model.AB$psi_at_0.1_kl50[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.5*Kmax_at_0.1)$y
  data.model.AB$psi_at_0.1_kl80[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.2*Kmax_at_0.1)$y
  data.model.AB$psi_at_0.1_kl88[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.12*Kmax_at_0.1)$y

}
comm.sp.vcurves <- do.call(rbind.data.frame, comm.predict.list)

## based on parameters with fits to observed data
obs.data.model.AB <- data.model.AB %>% subset(!is.na(data.A) & !is.na(data.B)) %>% droplevels()
obs.predict.list <- vector("list", nrow(data.model.AB))
names(obs.predict.list) <- data.model.AB$sp

for (i in 1:nrow(obs.data.model.AB)) {
  params <- obs.data.model.AB[i,]
  df <- data.frame(psi = seq(0, 5, length.out = 1000)) %>%
    mutate(sp = params$sp,
           k.predict = Exponential(A = params$data.A, B = params$data.B, psi = psi))
  obs.predict.list[[i]] <- df

  Kmax <- max(df$k.predict, na.rm = TRUE)
  Kmax_at_0.1 <- df$k.predict[which(round(df$psi, 3) == 0.1)]

  obs.data.model.AB$psi_kl20[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.8*Kmax)$y
  obs.data.model.AB$psi_kl50[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.5*Kmax)$y
  obs.data.model.AB$psi_kl80[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.2*Kmax)$y
  obs.data.model.AB$psi_kl88[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.12*Kmax)$y
  obs.data.model.AB$psi_kl95[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.05*Kmax)$y

  obs.data.model.AB$Kmax[i] <- Kmax
  obs.data.model.AB$Kmax_at_0.1[i] <- Kmax_at_0.1
  obs.data.model.AB$psi_at_0.1_kl50[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.5*Kmax_at_0.1)$y
  obs.data.model.AB$psi_at_0.1_kl80[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.2*Kmax_at_0.1)$y
  obs.data.model.AB$psi_at_0.1_kl88[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.12*Kmax_at_0.1)$y

}
obs.sp.vcurves <- do.call(rbind.data.frame, obs.predict.list)

save(bci.AB, file = file.path(results.folder, "bci.AB.Rdata"))
save(gap.models, file = file.path(results.folder, "gap.models.Rdata"))
save(data.model.AB, file = file.path(results.folder, "data.model.AB.Rdata"))
save(sp.soft.filled, file = file.path(results.folder, "sp.soft.filled.Rdata"))
save(sp.exp.param, file = file.path(results.folder, "sp.exp.param.Rdata"))

ggplot(subset(data.model.AB, !is.na(data.B)) %>% droplevels(),
       aes(y = data.B, x = SG100C_AVG)) + geom_point()
ggplot(subset(data.model.AB, !is.na(data.B)) %>% droplevels(),
       aes(y = data.B, x = psi_kl80, color = SG100C_AVG)) +
  geom_point(size = 3) +
  scale_color_viridis(direction = -1)
load(file = file.path("results/obs_vcurves_leaf.Rdata"))

pred.df.1 <- pred.df %>% left_join(data.model.AB, by = "sp") %>%
  subset(model_type == "Exponential")
obs.vcurves_fits <- ggplot(pred.df.1, aes(x = psi, y = predicted)) +
  geom_line(aes(group = sp, col = SG100C_AVG), size = 0.5) +
  xlab("Leaf Water Potential (-MPa)") +
  ylab(expression(atop(italic('K')['max, leaf'], ""^(mmol*~m^-2*~s^-1*~MPa^-1)))) +
  # scale_color_gradient(low = "blue", high = "Red")
  scale_color_viridis(name = expression("WSG ("*g*cm^-3*")"), direction = -1) +
  theme(legend.position = c(0.8, 0.65),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent"))
ggsave(plot = obs.vcurves_fits, file.path(figures.folder.kleaf, paste0("obs.range_vcurves_fits_ggplot.tiff")),
       device = "tiff", height = 3, width = 3.5, units='in')

obs.sp.vcurves.1 <- obs.sp.vcurves %>% left_join(data.model.AB, by = "sp") %>%
  group_by(sp) %>%
  mutate(k.predict.percent = 100 - range01(k.predict)*100)

obs.vcurves_full.pre <- ggplot(obs.sp.vcurves.1, aes(x = psi, y = k.predict)) +
  xlab("Leaf Water Potential (-MPa)") +
  ylab(expression(atop(italic('K')['max, leaf'], ""^(mmol*~m^-2*~s^-1*~MPa^-1)))) +
  # scale_color_gradient(low = "blue", high = "Red")
  theme(legend.position = c(0.8, 0.65),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent"))
obs.vcurves_full.wsg <- obs.vcurves_full.pre +
  geom_line(aes(group = sp, col = SG100C_AVG), size = 0.5) +
  scale_color_viridis(name = expression("WSG ("*g*cm^-3*")"), direction = -1, option = "B")
ggsave(plot = obs.vcurves_full.wsg, file.path(figures.folder.kleaf, paste0("obs_vcurves_fits_ggplot_wsg.tiff")),
       device = "tiff", height = 3, width = 3.5, units='in')
obs.vcurves_full.lma <- obs.vcurves_full.pre +
  geom_line(aes(group = sp, col = LMALAM_AVD.filled), size = 0.5) +
  scale_color_viridis(name = expression("LMA ("*g*cm^-2*")"), direction = -1, option = "D")
ggsave(plot = obs.vcurves_full.lma, file.path(figures.folder.kleaf, paste0("obs_vcurves_fits_ggplot_lma.tiff")),
       device = "tiff", height = 3, width = 3.5, units='in')

obs.plccurves.pre <- ggplot(obs.sp.vcurves.1, aes(x = psi, y = k.predict.percent)) +
  xlab("Leaf Water Potential (-MPa)") +
  ylab("Loss of Conductivity (%)") +
  # scale_color_gradient(low = "blue", high = "Red")
  theme(legend.position = c(0.8, 0.4),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent"))

obs.plccurves.wsg <- obs.plccurves.pre +
  geom_line(aes(group = sp, col = SG100C_AVG), size = 0.5) +
  scale_color_viridis(name = expression("WSG ("*g*cm^-3*")"), direction = -1, option = "B")
ggsave(plot = obs.plccurves.wsg, file.path(figures.folder.kleaf, paste0("obs_plccurves_fits_ggplot_wsg.tiff")),
       device = "tiff", height = 3, width = 3.5, units='in')
obs.plccurves.lma <- obs.plccurves.pre +
  geom_line(aes(group = sp, col = LMALAM_AVD.filled), size = 0.5) +
  scale_color_viridis(name = expression("LMA ("*g*cm^-2*")"), direction = -1, option = "D")
ggsave(plot = obs.plccurves.lma, file.path(figures.folder.kleaf, paste0("obs_plccurves_fits_ggplot_lma.tiff")),
       device = "tiff", height = 3, width = 3.5, units='in')

comm.sp.vcurves.1 <- comm.sp.vcurves %>%
  left_join(data.model.AB, by = "sp") %>%
  group_by(sp) %>%
  mutate(k.predict.percent = 100 - range01(k.predict)*100)

save(comm.sp.vcurves.1, file = file.path(results.folder, "comm.sp.vcurves.1.Rdata"))
save(obs.sp.vcurves.1, file = file.path(results.folder, "obs.sp.vcurves.1.Rdata"))

comm.vcurves.wsg <- obs.vcurves_full.wsg %+% comm.sp.vcurves.1 +
  ylim(c(0, 15))
ggsave(plot = comm.vcurves.wsg, file.path(figures.folder.kleaf, paste0("comm_vcurves_fits_ggplot.wsg.tiff")),
       device = "tiff", height = 3, width = 3.5, units='in')
comm.vcurves.lma <- obs.vcurves_full.lma %+% comm.sp.vcurves.1 +
  ylim(c(0, 15))
ggsave(plot = comm.vcurves.lma, file.path(figures.folder.kleaf, paste0("comm_vcurves_fits_ggplot.lma.tiff")),
       device = "tiff", height = 3, width = 3.5, units='in')

comm.plccurves.wsg <- obs.plccurves.wsg %+% comm.sp.vcurves.1
ggsave(plot = comm.plccurves.wsg, file.path(figures.folder.kleaf, paste0("comm_plccurves_fits_ggplot.wsg.tiff")),
       device = "tiff", height = 3, width = 3.5, units='in')
comm.plccurves.lma <- obs.plccurves.lma %+% comm.sp.vcurves.1
ggsave(plot = comm.plccurves.lma, file.path(figures.folder.kleaf, paste0("comm_plccurves_fits_ggplot.lma.tiff")),
       device = "tiff", height = 3, width = 3.5, units='in')

data.model.AB.onlyboth <- sp.exp.param %>% rename(data.A = A, data.B = B) #%>%
  # subset(!is.na(model.A) | !is.na(model.B))
formula = y ~ x
plot.B <- ggplot(data.model.AB.onlyboth, aes(x = model.B, y = data.B)) +
  geom_abline(intercept = 0, slope = 1, color = "dodgerblue") +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 0.8, size = 3) +
  geom_text(data = data.frame(x = 4.2, y = 4.4, label = "1:1"), aes(x = x, y = y, label = label), size =5) +
  xlim(range(c(data.model.AB.onlyboth$model.B, data.model.AB.onlyboth$data.B), na.rm = TRUE)) +
  ylim(range(c(data.model.AB.onlyboth$model.B, data.model.AB.onlyboth$data.B), na.rm = TRUE)) +
  ylab(expression(italic(B)['Fitted to Data'])) +
  xlab(expression(italic(B)['Predicted from WSG, LMA'])) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.05, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste0("P", ifelse(p.value < 0.001, " < 0.001",
                                                  paste0(" = ", round(..p.value.., digits = 4))))),
                  npcx = 0.05, npcy = 0.8, size = 4)
ggsave(plot = plot.B, file.path(figures.folder.kleaf, paste0("B_data_vs_model.tiff")),
      device = "tiff", height = 2.2, width = 2.2, units='in')

plot.A <- ggplot(data.model.AB.onlyboth, aes(x = model.A, y = data.A)) +
  geom_abline(intercept = 0, slope = 1, color = "dodgerblue") +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 0.8, size = 3) +
  geom_text(data = data.frame(x = 8, y = 8.5, label = "1:1"), aes(x = x, y = y, label = label), size = 5) +
  xlim(range(c(data.model.AB.onlyboth$model.A, data.model.AB.onlyboth$data.A), na.rm = TRUE)) +
  ylim(range(c(data.model.AB.onlyboth$model.A, data.model.AB.onlyboth$data.A), na.rm = TRUE)) +
  ylab(expression(italic(A)['Fitted to Data'])) +
  xlab(expression(italic(A)['Predicted from LMA & '*italic(B)])) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.05, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste0("P", ifelse(p.value < 0.001, " < 0.001",
                                                 paste0(" = ", round(..p.value.., digits = 4))))),
                  npcx = 0.05, npcy = 0.8, size = 4)
ggsave(plot = plot.A, file.path(figures.folder.kleaf, paste0("A_data_vs_model.tiff")),
       device = "tiff", height = 2.2, width = 2.2, units='in')

plot.AB <- cowplot::plot_grid(plot.A, plot.B, labels = c('A', 'B'),
                                  label_size = 14, ncol = 2, rel_widths = c(1.1, 1))
ggsave("A_B_data_vs_model.tiff", plot = plot.AB, path =
         file.path(figures.folder.kleaf), device = "tiff", height = 2.2, width = 4.6, units ='in')
ggsave("A_B_data_vs_model.jpeg", plot = plot.AB, path =
         file.path(figures.folder.kleaf), device = "jpeg", height = 2.2, width = 4.6, units ='in')

###*****
## Plot K by PSI using fitted exponential curves
###*****

# Define colour pallete
pal = colorRampPalette(c("blue", "red"))
# pal = colorRampPalette(c("red", "grey", "blue")) #"#D60F73", "#8119D1",
# pal = colorRampPalette(c("#7F0000", "#00007F")) #colorRampPalette(c("red", "blue"))


# Use the following line with RColorBrewer
# pal = colorRampPalette(cols)

sp.exp.param.plot <- sp.exp.param %>%
  mutate(order.chl = findInterval(Chl, sort(Chl)),
         order.tlp = findInterval(TLP, sort(TLP)),
         order.wsg = findInterval(SG100C_AVG, sort(SG100C_AVG)),
         order.decilvl = findInterval(DeciLvl, sort(DeciLvl)),
         order.lma = findInterval(LMA, sort(LMA)),
         order.B = findInterval(B, sort(B)))

k.raw.data <- read.table(paste0("data-raw/traits/HydraulicTraits_Kunert/Panama_raw_hydraulics_data.csv"), header = T, sep = ",") %>%
  left_join(sp.exp.param.plot, by = "sp")

## choose among alternatives of which data or subset to plot and which variable
df.plot <- sp.exp.param.plot
df.name <- "kmax_psi_data_fitted_AB"
col.var <- "LMA"
order.col.var <- "order.lma"
legend.col.var <- "LMA"

# col.var <- "B"
# order.col.var <- "order.B"
# legend.col.var <- "B"

df.plot <- bci.AB %>%
  subset(B >= 0 & A >= 0) %>%
  # subset(sp %in% sp.exp.param.plot$sp) %>% droplevels() %>%
  # subset(sp %in% unique(iso.1.3.join$sp)) %>% droplevels() %>%
  left_join(deci, by = "sp") %>%
  mutate(order.wsg = findInterval(SG100C_AVG, sort(SG100C_AVG)),
         order.lma = findInterval(LMALAM_AVD, sort(LMALAM_AVD)),
         order.decilvl = findInterval(DeciLvl, sort(DeciLvl))) %>%
  mutate(deciduousness.label.2 = recode_factor(as.factor(deciduous), `E` = "Evergreen", `DB` = "Brevi\nDeciduous",
                                               `DF` = "Facultative\nDeciduous", `DO` = "Obligate\nDeciduous")) %>%
  transform(deciduousness.label.2 = factor(deciduousness.label.2,
                                           levels = c("Evergreen", "Brevi\nDeciduous",
                                                      "Facultative\nDeciduous", "Obligate\nDeciduous"), ordered = TRUE)) %>%
  mutate(Kmax.predict = Exponential(A = A, B = B, psi = 0))

# jpeg(file.path(figures.folder.kleaf, paste0("Sp.Gravity_by_deci.jpeg")),
#      width = 5, height = 2.5, units = "in", pointsize = 10,
#      quality = 100, res = 300)
# par(mar = c(5, 5, 2, 2))
# boxplot(SG100C_AVG ~ deciduousness.label.2, data = df.plot, col = alpha(1:4, 0.8), alpha = 0.3,
#         ylab = expression('Wood Specific Gravity'['100C']), notch = TRUE, boxwex = 0.5)
# dev.off()
# jpeg(file.path(figures.folder.kleaf, paste0("LMALAM_by_deci.jpeg")),
#      width = 5, height = 3.5, units = "in", pointsize = 10,
#      quality = 100, res = 300)
# par(mar = c(5, 5, 2, 2))
# boxplot(LMALAM_AVD ~ deciduousness.label.2, data = df.plot, col = alpha(1:4, 0.8), alpha = 0.3,
#         ylab = expression('LMA'['Lamina']), notch = TRUE, boxwex = 0.5)
# dev.off()

var.list <- list(df.name = c("predicted_AB", "predicted_AB_for_data_sp"), #, "predicted_AB_for_iso_sp"
                 col.var = c("LMALAM_AVD", "SG100C_AVG", "DeciLvl", "NILL"),
                 order.col.var = c("order.lma", "order.wsg", "order.decilvl", "NILL"),
                 legend.col.var = c(expression("LMA ("*g*m^-2*")"),
                                    expression("WSG ("*g*cm^-3*")"), "Deciduoousness", ""),
                 std.k = c("", "std.k.sp", "std.k.comm"),
                 ylim.k = c(7.5, 100, 100),
                 legend.pos = c("topright", "bottomright", "bottomright"))

for (i in 1:length(var.list$df.name)) {
  df.name <- var.list$df.name[i]
  ### Data Prep
  df.plot <- bci.AB %>%
    subset(B >= 0 & A >= 0)

  if (df.name == "predicted_AB_for_data_sp") {
    df.plot <- df.plot %>%
      subset(sp %in% sp.exp.param.plot$sp) %>% droplevels()
  }
  # if (df.name == "predicted_AB_for_iso_sp") {
  #   df.plot <- df.plot %>%
  #     subset(sp %in% unique(iso.1.3.join$sp)) %>% droplevels()
  # }
  df.plot <- df.plot %>%
    left_join(deci, by = "sp") %>%
    mutate(order.wsg = findInterval(SG100C_AVG, sort(SG100C_AVG)),
           order.lma = findInterval(LMALAM_AVD, sort(LMALAM_AVD)),
           order.decilvl = findInterval(DeciLvl, sort(DeciLvl))) %>%
    mutate(deciduousness.label.2 = recode_factor(as.factor(deciduous), `E` = "Evergreen", `DB` = "Brevi\nDeciduous",
                                                 `DF` = "Facultative\nDeciduous", `DO` = "Obligate\nDeciduous")) %>%
    transform(deciduousness.label.2 = factor(deciduousness.label.2,
                                             levels = c("Evergreen", "Brevi\nDeciduous",
                                                        "Facultative\nDeciduous", "Obligate\nDeciduous"), ordered = TRUE)) %>%
    mutate(Kmax.predict = Exponential(A = A, B = B, psi = 0)) %>% droplevels()

  ###
  Kmax.predict.max <- max(df.plot$Kmax.predict, na.rm = TRUE)

  # if(df.name == "predicted_AB_for_iso_sp") {
  #   xlim.to.plot <- c(-1, 3)
  # } else {
    xlim.to.plot <- c(0, 3)
  # }

  for (j in 1:length(var.list$col.var)) {
    col.var <- var.list$col.var[j]
    order.col.var <- var.list$order.col.var[j]
    legend.col.var <- var.list$legend.col.var[j]
    for (k in 1:length(var.list$std.k)) {
      std.k <- var.list$std.k[k]
      ylim.k <- var.list$ylim.k[k]
      legend.pos <- var.list$legend.pos[k]

      jpeg(file.path(figures.folder.kleaf, paste0(std.k, "kmax_by_psi_color_by_", col.var, "_", df.name, ".jpeg")),
           width = 2.7, height = 2.7, units = "in", pointsize = 10,
           res = 300) #quality = 100,
      if (std.k == "") {
        par(mar = c(4, 4.5, 1.5, 1.5))
        plot(1, type = "n", xlab = "Leaf Water Potential (-MPa)",
             xlim = xlim.to.plot, ylim = c(0, ylim.k))
        mtext(side = 2, text = "Leaf Hydraulic Conductance", line = 3)
        mtext(side = 2, text = expression("(mmol "*m^-2*s^-1*MPa^-1*")"), line = 2)
      } else {
        par(mar = c(4.5, 4.5, 1.5, 1.5))
        plot(1, type = "n", xlab = "Leaf Water Potential (-MPa)", ylab =
               expression("Loss of Conductance (%)"),
             xlim = xlim.to.plot, ylim = c(0, ylim.k))
      }
      for (m in 1:nrow(df.plot)) {
        params <- df.plot[m, ]
        if(df.name == "kmax_psi_data_fitted_AB"){
          col = "gray85"; pch = 20
        }
      }
      for (n in 1:nrow(df.plot)) {
        params <- df.plot[n, ]
        df <- data.frame(psi = seq(0, 3, length.out = 100)) %>%
          mutate(k.predict = Exponential(A = params$A, B = params$B, psi = psi)) %>%
          cbind.data.frame(params, row.names = NULL)
        if(std.k == "std.k.sp") {
          df <- df %>% mutate(k.predict = 100 - range01(k.predict)*100)
        }
        if(std.k == "std.k.comm") {
          df <- df %>% mutate(k.predict = 100 - k.predict/Kmax.predict.max*100)
        }
        if(col.var == "NILL") {
          lines(k.predict ~ psi, data = df) # "darkorange"
        }
        # Rank variable for colour assignment
        if(col.var == "deci") {
          lines(k.predict ~ psi, data = df, col = deciduousness) # "darkorange"
        } else {
          lines(k.predict ~ psi, data = df, col = pal(nrow(df.plot))[df.plot[n, order.col.var]], size = 0.3) # "darkorange"
        }
        if(df.name == "predicted_AB_for_iso_sp") {
          text(labels = params$sp, x = -0.2, y = df$k.predict[df$psi == df$psi[1]],
               col = pal(nrow(df.plot))[df.plot[n, order.col.var]], cex = 0.3)
        }
      }

      if(col.var == "deci") {
        ## using levels in deciduousness used by color
        legend(x = legend.pos, legend = levels(df.plot$deciduousness),
               col = 1:4, pch=19, bty = "n")
      } else if(col.var != "NILL") {
        lgd_ <- rep(NA, 4)
        lgd_[c(1,4)] <- c(round(sort(range(df.plot[, col.var], na.rm = TRUE),
                            decreasing = TRUE), 1))
        fill.col <- sort(pal(4), decreasing = TRUE) #colorRampPalette(colors = c('red', 'blue'))(4)
         legend(x = legend.pos, title = legend.col.var, fill=fill.col,
               legend = lgd_, bty = "n", y.intersp = 0.5)
        # legend("topright",
        #        legend = lgd_,
        #        fill = colorRampPalette(colors = c('black','red3','grey96'))(11),
        #        border = NA,
        #        y.intersp = 0.5,
        #        cex = 2, text.font = 2)
      }
      dev.off()
    }
  }
}


#******************************************************
## Predictors of KmaxS-SWP relationship --------
#******************************************************

figures.folder.kstem <- paste0("figures/PhenoDemoTraitsPsi/kmax_by_psi/Stem")
if(!dir.exists(file.path(figures.folder.kstem))) {dir.create(file.path(figures.folder.kstem))}

sp.exp.param.stem <- hyd %>%
  rename(Kmax = KmaxS, A = vc_a, B = vc_b, A_se = vc_a_se, B_se = vc_b_se,
         LMA_sd = lma_gm2sd, LMA_n = lma_n,
         WD_sd = xylem_den_sd) %>%
  select(sp, Kmax, A, B, A_se, B_se, WD, LMA, LMA_sd, LMA_n, WD_sd) %>%
  left_join(deci %>% select(-sp4, -deciduousness.label), by = "sp") %>%
  left_join(bci.traits %>%
              select(sp, SG100C_AVG, LMALEAF_AVD, LMALAM_AVD), by = "sp") %>%
  left_join(data.model.AB %>% select(sp, model.B, model.A, data.A, data.B), by = "sp")
range(sp.exp.param.stem$Kmax, na.rm = TRUE)
# 0.87 8.83

## ******
### PCA to check which variables strongly correlate with Parameters A & B
## ******

df.pca.stem <- sp.exp.param.stem %>% remove_rownames %>% column_to_rownames(var = "sp") %>%
  select_if(is.numeric) ## more species without Kmax data %>% select(-SafetyMargin.p50, -p50_Kmax, -Kmax)
# if (diff(range(df.pca$root.95, na.rm = TRUE)) == 0) {
df.pca.stem <- df.pca.stem %>% subset(complete.cases(df.pca.stem))
result.pca.stem <- prcomp(df.pca.stem, center = TRUE, scale = TRUE)

pdf(file.path(figures.folder.kstem ,"kamx_by_psi_params_pca.pdf"), height = 6, width = 6)
biplot(result.pca.stem, choices = 1:2, pc.biplot = TRUE, main = "")
dev.off()

## single explanatory variable
var.y <- c("A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B","B", "B", "B")
var.x <- c("WD", "LMA", "SG100C_AVG", "LMALAM_AVD", "LMALEAF_AVD", "B","data.A", "model.A",
           "WD", "LMA", "SG100C_AVG", "LMALAM_AVD", "LMALEAF_AVD", "A", "data.B", "model.B")
for(i in 1:length(var.x)){
  jpeg(file.path(figures.folder.kstem, paste0("kmax_by_psi_", var.y[i], "_by_", var.x[i], ".jpeg")),
       width = 4, height = 4, units = "in", pointsize = 10,
       quality = 100, res = 300)
  df <- sp.exp.param.stem[!is.na(sp.exp.param.stem[, var.x[i]]),]
  df$x <- df[, var.x[i]]
  df$y <- df[, var.y[i]]
  model1 <- lm(y ~ poly(x, 2), data = df)
  plot(y ~ x, data = df, xlab = var.x[i], ylab = var.y[i])
  x0 <- seq(min(df$x, na.rm = TRUE), max(df$x, na.rm = TRUE), length = 100)  ## prediction grid
  y0 <- predict(model1, newdata = list(x = x0))  ## predicted values
  lines(x0, y0, col = 2)
  ## rounded coefficients for better output
  cf <- round(coef(model1), 2)
  ## sign check to avoid having plus followed by minus for negative coefficients
  eq <- paste0("y = ", cf[1],
               ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), "x",
               ifelse(sign(cf[3])==1, " + ", " - "), abs(cf[3]), "x2")
  r2 <- paste0("R-squared = ", round(summary(model1)$r.squared, 2),
               ", P-val = ", round(lmp(model1), 2))
  mtext(eq, 3, line = -2)
  mtext(r2, 3, line = -4)
  dev.off()
}
# A related with LMA, then WSG & WD
# B highly related with LMA, then WSG, then LMALeaf

## Best R2 have two variables
k_by_psi.models.stem <- vector(mode = "list", length = 2)
k_by_psi.models.stem[[1]] <- lm(B ~ polym(LMA, degree = 2, raw = TRUE), data = sp.exp.param.stem)
summary(k_by_psi.models.stem[[1]])
# With only LMA Adjusted R-squared:  0.173, 0.043
# With LMA & SG100C_AVG, Adjusted R-squared:  0.1701;  p-value: 0.147
## with only SG100C_AVG, R-squared:  0.1354 ; p-value: 0.09001
## but not related with LMALEAF_AVD, so this is not helpful to generalise

k_by_psi.models.stem[[2]]<- lm(A ~ polym(LMA, SG100C_AVG, degree = 2, raw = TRUE), data = sp.exp.param.stem)
summary(k_by_psi.models.stem[[2]])
# Adjusted R-squared:  0.3933,  p-value: 0.01622
## but not related with LMALEAF_AVD or LMALAM_AVD, so this is not helpful to generalise

## But LMA can be estimated from LMALEAF_AVD so filling up some gaps in them
k_by_psi.models.stem[[3]] <- lm(LMA ~ LMALAM_AVD, data = sp.exp.param.stem)
summary(k_by_psi.models.stem[[3]])
# Adjusted R-squared:  0.639; p-value: 0.0003513
plot(LMA ~ LMALEAF_AVD, data = sp.exp.param.stem)
save(k_by_psi.models.stem, file = file.path(results.folder, "k_by_psi.models.stem.Rdata"))

###**********
## Predict A & B from soft traits and the above models
###**********

bci.AB.stem <- bci.traits %>% select(sp, WSG_CHAVE, SG60C_AVG, SG100C_AVG, LMALEAF_AVD, LMALAM_AVD)
bci.AB.stem <- bci.AB.stem %>% mutate(LMALAM_AVD = ifelse(is.na(LMALAM_AVD),
                             predict(gap.models$LMA.LAM.LEAF, newdata = bci.AB.stem), LMALAM_AVD))
bci.AB.stem <- bci.AB.stem %>% mutate(LMA = predict(k_by_psi.models.stem[[3]], newdata = bci.AB.stem))
bci.AB.stem <- bci.AB.stem %>% mutate(B = predict(k_by_psi.models.stem[[1]], newdata = bci.AB.stem),
                                      A = predict(k_by_psi.models.stem[[2]], newdata = bci.AB.stem))

###**********
## Do params A & B predicted from soft traits match those fitted to data?
###**********
data.model.AB.stem <- bci.AB.stem %>%
  subset(B >= 0 & A >= 0) %>%
  rename(model.A = A, model.B = B) %>%
  left_join(sp.exp.param.stem %>% select(sp, A, B) %>%
              rename(data.A = A, data.B = B), by = "sp")

for (i in 1:nrow(data.model.AB.stem)) {
  params <- data.model.AB.stem[i,]
  df <- data.frame(psi = seq(0, 5, length.out = 1000)) %>%
    mutate(k.predict = Exponential(A = params$model.A, B = params$model.B, psi = psi))

  Kmax <- max(df$k.predict, na.rm = TRUE)
  Kmax_at_0.1 <- df$k.predict[which(round(df$psi, 3) == 0.1)]

  data.model.AB.stem$psi_kl20[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.8*Kmax)$y
  data.model.AB.stem$psi_kl50[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.5*Kmax)$y
  data.model.AB.stem$psi_kl80[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.2*Kmax)$y
  data.model.AB.stem$psi_kl95[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.05*Kmax)$y

  data.model.AB.stem$Kmax[i] <- Kmax
  data.model.AB.stem$Kmax_at_0.1[i] <- Kmax_at_0.1
  data.model.AB.stem$psi_at_0.1_kl50[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.5*Kmax_at_0.1)$y
  data.model.AB.stem$psi_at_0.1_kl80[i] <- -approx(x = df$k.predict, y = df$psi, xout=0.2*Kmax_at_0.1)$y

}

save(data.model.AB.stem, file = file.path(results.folder, "data.model.AB.stem.Rdata"))

ggplot(data.model.AB.stem, aes(x = model.B, y = data.B)) +
  geom_abline(intercept = 0, slope = 1, color = "dodgerblue") +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 0.8, size = 3) +
  geom_text(data = data.frame(x = 4.6, y = 4.0, label = "1:1"), aes(x = x, y = y, label = label), size =5) +
  xlim(range(c(data.model.AB.stem$model.B, data.model.AB.stem$data.B), na.rm = TRUE)) +
  ylim(range(c(data.model.AB.stem$model.B, data.model.AB.stem$data.B), na.rm = TRUE)) +
  ylab(expression(italic(B)['Fitted to Data'])) +
  xlab(expression(italic(B)['Predicted from LMA'])) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.05, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.05, npcy = 0.8, size = 4)
ggsave(file.path(figures.folder.kstem, paste0("B_data_vs_model.jpeg")),
       device = "jpeg", height = 2.2, width = 2.2, units='in')

ggplot(data.model.AB.stem, aes(x = model.A, y = data.A)) +
  geom_abline(intercept = 0, slope = 1, color = "dodgerblue") +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 0.8, size = 3) +
  geom_text(data = data.frame(x = 4.6, y = 4.0, label = "1:1"), aes(x = x, y = y, label = label), size = 5) +
  xlim(range(c(data.model.AB.stem$model.B, data.model.AB.stem$data.B), na.rm = TRUE)) +
  ylim(range(c(data.model.AB.stem$model.B, data.model.AB.stem$data.B), na.rm = TRUE)) +
  ylab(expression(italic(A)['Fitted to Data'])) +
  xlab(expression(italic(A)['Predicted from LMA & WSG'])) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.05, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.05, npcy = 0.8, size = 4)
ggsave(file.path(figures.folder.kstem, paste0("A_data_vs_model.jpeg")),
       device = "jpeg", height = 2.2, width = 2.2, units='in')
## So poor predictors

##********************
## Plot K by Psi Stem
##********************
# df.name <- "predicted_AB_for_data_sp"
df.name <- "predicted_AB"

if (df.name == "predicted_AB_for_data_sp") {
  df.plot.stem <- bci.AB.stem %>%
    left_join(deci %>% select(-sp4, -deciduousness.label), by = "sp")
} else if(df.name == "predicted_AB") {
  df.plot.stem <- sp.exp.param.stem
}
## Plotting DataAB
df.plot.stem <- df.plot.stem %>%
  subset(B >= 0 & A >= 0) %>%
  mutate(order.wsg = findInterval(SG100C_AVG, sort(SG100C_AVG)),
         order.lma = findInterval(LMALAM_AVD, sort(LMALAM_AVD)),
         order.decilvl = findInterval(DeciLvl, sort(DeciLvl))) %>%
  mutate(deciduousness.label.2 = recode_factor(as.factor(deciduous), `E` = "Evergreen", `DB` = "Brevi\nDeciduous",
                                               `DF` = "Facultative\nDeciduous", `DO` = "Obligate\nDeciduous")) %>%
  transform(deciduousness.label.2 = factor(deciduousness.label.2,
                                           levels = c("Evergreen", "Brevi\nDeciduous",
                                                      "Facultative\nDeciduous", "Obligate\nDeciduous"), ordered = TRUE)) %>%
  mutate(Kmax.predict = Exponential(A = A, B = B, psi = 0))

# col.var <- "NILL"
# order.col.var <- "NILL"
# legend.col.var <- "NILL"
col.var <- "LMALAM_AVD"
order.col.var <- "order.lma"
legend.col.var <- "LMALAM_AVD"
# col.var <- "SG100C_AVG"
# order.col.var <- "order.wsg"
# legend.col.var <- "Specific Gravity"
# col.var <- "DeciLvl"
# order.col.var <- "order.decilvl"
# legend.col.var <- "DeciLvl"
# col.var <- "deci"
# legend.col.var <- "Deciduoousness"
std.k <- ""; ylim.k = 20; if(df.name == "predicted_AB") {ylim.k = max(df.plot.stem$Kmax, na.rm = TRUE)}
# std.k <- "std.k.sp"; ylim.k = 1
# std.k <- "std.k.comm"; ylim.k = 1
Kmax.predict.max.stem <- max(df.plot.stem$Kmax.predict, na.rm = TRUE)

jpeg(file.path(figures.folder.kstem, paste0(std.k, "kmax_by_psi_color_by_", legend.col.var, "_", df.name, ".jpeg")),
     width = 2.7, height = 2.7, units = "in", pointsize = 10,
     quality = 100, res = 300)
if (std.k == "") {
  par(mar = c(4, 4.5, 1.5, 1.5))
  plot(1, type = "n", xlim = xlim.to.plot, ylim = c(0, ylim.k), ann = FALSE)
  mtext(side = 2, text = "Stem-Specific Hydraulic Conductivity", line = 3)
  mtext(side = 2, text = expression("(Kg "*m^-1*s^-1*MPa^-1*")"), line = 2)
  mtext(side = 1, text = "Stem Water Potential (-MPa)", line = 2)
} else {
  par(mar = c(4.5, 4.5, 1.5, 1.5))
  plot(1, type = "n", xlim = xlim.to.plot, ylim = c(0, ylim.k), ann = FALSE)
  mtext(side = 2, text = "Std. Stem-Specific", line = 3)
  mtext(side = 2, text = expression("Hydraulic Conductivity"), line = 2)
  mtext(side = 1, text = "Stem Water Potential (-MPa)", line = 2)
}
for (i in 1:nrow(df.plot)) {
  params <- df.plot[i, ]
  if(df.name == "kmax_psi_data_fitted_AB"){
    col = "gray85"; pch = 20
  }
}
for (i in 1:nrow(df.plot.stem)) {
  params <- df.plot.stem[i, ]
  df <- data.frame(psi = seq(0, 3, length.out = 100)) %>%
    mutate(k.predict = Exponential(A = params$A, B = params$B, psi = psi)) %>%
    cbind.data.frame(params, row.names = NULL)
  if(std.k == "std.k.sp") {
    df <- df %>% mutate(k.predict = range01(k.predict))
  }
  if(std.k == "std.k.comm") {
    df <- df %>% mutate(k.predict = k.predict/Kmax.predict.max.stem)
  }
  if(col.var == "NILL") {
    lines(k.predict ~ psi, data = df, alpha = 0.7) # "darkorange"
  }
  # Rank variable for colour assignment
  if(col.var == "deci") {
    lines(k.predict ~ psi, data = df, col = deciduousness) # "darkorange"
  } else {
    lines(k.predict ~ psi, data = df, col = pal(nrow(df.plot.stem))[df.plot.stem[i, order.col.var]]) # "darkorange"
  }
  if(df.name == "predicted_AB_for_iso_sp") {
    text(labels = params$sp, x = -0.2, y = df$k.predict[df$psi == df$psi[1]],
         col = pal(nrow(df.plot.stem))[df.plot.stem[i, order.col.var]], cex = 0.5)
  }
}
if(col.var == "deci") {
  ## using levels in deciduousness used by color
  legend("topright", legend = levels(df.plot.stem$deciduousness),
         col = 1:4, pch=19, bty = "n")
} else if(col.var != "NILL") {
  legend("topright", title = legend.col.var, col=pal(2), pch=19,
         legend=c(round(sort(range(df.plot.stem[, col.var], na.rm = TRUE),
                             decreasing = TRUE), 1)), bty = "n")
}
dev.off()

#******************************************************
### Load LAI data ------
#******************************************************
figures.folder.phen <- paste0("figures/PhenoDemoTraitsPsi/leaf_fall")
if(!dir.exists(file.path(figures.folder.phen))) {dir.create(file.path(figures.folder.phen))}

lai <- read.csv(file.path("data-raw/BCI_lai_matteo.csv")) %>%
  mutate(date = as.Date(date, format = "%d-%b-%y"),
         mon = format(date, "%m"),
         month = format(date, "%b"),
         month.plot = factor(month,
                             levels=unique(month[order(mon)]),
                             ordered=TRUE)) %>%
  group_by(month, month.plot, mon) %>%
  summarise(lai = mean(lai, na.rm = TRUE),
            se = mean(lai, na.rm = TRUE)) %>%
  mutate(group = "group") %>%
  mutate(date = as.Date(paste0("01-", month), format("%d-%b")),
         doy = as.numeric(format(date, "%j")))

lai.plot.1 <- ggplot(lai, aes(y = lai, x = month.plot)) +
  geom_point(size = 3) + geom_line(aes(group = group)) +
  ylab("LAI") + xlab("Month")
ggsave(("LAI_by_month_BCI.jpeg"),
       plot = lai.plot.1, file.path(figures.folder.phen), device = "jpeg", height = 3, width = 4.5, units='in')

lai.plot.2 <- ggplot(lai, aes(y = lai, x = doy)) +
  geom_point(size = 3) + geom_line(aes(group = group)) +
  ylab("LAI") + xlab("DOY") + theme(plot.margin = margin(1,1,1,1, "cm"))
ggsave(("LAI_by_DOY_BCI.jpeg"),
       plot = lai.plot.2, file.path(figures.folder.phen), device = "jpeg", height = 3, width = 4.5, units='in')

#******************************************************
### Load Weekly Leaf-fall data ------
#******************************************************
leaf.fall <- read.csv(file.path("data-raw/traits/Wright_Osvaldo_BCI_weekly_leaf-fall_data/Rutuja.csv")) %>%
  rename(date = mean_date,
         sp4 = sp) %>% mutate(date = as.Date(date)) %>%
  subset(site != "San_Lorenzo") %>%
  # subset(site == "BCI50-ha") %>%
  # mutate(site = "BCI") %>%
  left_join(deci %>% dplyr::select(sp, sp4), by = "sp4") %>%
  mutate(sp_site = paste(sp, site, sep = "_"))

f1 <- ggplot(leaf.fall, aes(x = date, y = leaf_gm)) +
  geom_line(aes(group = sp, color = sp), show.legend = FALSE)
ggsave(("leaf.fall.daily_BCI.jpeg"),
       plot = f1, file.path(figures.folder.phen), device = "jpeg", height = 4.5, width = 9, units='in')

## Need to fill gaps in dates and then get an interpolated estimate of leaf fall for those days
## So first converting the weekly sums to leaf_fall rate for the past week
## i.e. leaf_fall rate per day of the census interval = sum of leaf fall durign census interval/no. of days in the census interval
## to get leaf_fall on each day of the interval, the above rate is actually treated as
## the rate observed on the census date (even though it is the mean rate of the interval)

df.sp.site <- split(leaf.fall, f = list(leaf.fall$sp_site), drop = TRUE)

# leaf.fall.full <- vector(mode = "list", length = length(df.sp.site))
# names(leaf.fall.full) <- names(df.sp.site) # "psi.p50.g1", "psi.p50.g2"

fill.day.gaps <- function(df) {
  full.date.df <- data.frame(date = seq(from = min(df$date, na.rm = TRUE),
                                        to = max(df$date, na.rm = TRUE), by = 1)) %>%
    mutate(day_number = as.numeric(difftime(date, min(df$date, na.rm = TRUE) - 1)))
  df.1 <- df %>% arrange(date) %>%
    mutate(n.days = as.numeric(date - lag(date))) %>%
    #  leaf_fall rate per day, mean for the census interval
    mutate(leaf_gm_rate = leaf_gm/n.days) %>%# day
    full_join(full.date.df, by = "date") %>% arrange(date, site)
  return(df.1)
}
# leaf_fall on each day of the interval
leaf.fall.daygaps <- lapply(df.sp.site, fill.day.gaps)

## Interpolating from weely sums to daily leaf_gm
leaf.interp.approx <- function(df) {
  x <- df$day_number
  y <- df$leaf_gm_rate
  xout <- df$day_number[is.na(df$leaf_gm_rate)]
  # yout <- approx(x, y, xout, method = "linear")
  ## method = "constant" would be more parsimonious
  yout <- approx(x, y, xout, method = "constant")
  df.1 <- df %>%
    left_join(data.frame(day_number = yout$x, leaf_gm.int.raw = yout$y), by = "day_number") %>%
    ## filling interpolation gap on the day of the census
    mutate(leaf_gm.int.raw = ifelse(is.na(leaf_gm_rate), leaf_gm.int.raw, leaf_gm_rate),
           sp = df$sp[1],
           site = df$site[1])
  return(df.1)
}
## Absolute biomass among the two sites on BCI are not comparable
## So normalising them.
## Here we are interested in the seasonal pattern of leaffall
## So rescaling weekly leaf-fall as a fraction of the total leaf fall of the year--beginning DOY 120
set.k = 7
leaf.fall.int <- lapply(lapply(leaf.fall.daygaps, leaf.interp.approx),
                        as.data.frame) %>%
  bind_rows(.id = "sp.site") %>%
  group_by(sp, date) %>%
  summarise(leaf_gm.int.raw = sum(leaf_gm.int.raw)) %>%
  mutate(doy = as.numeric(format(date, "%j")),
         year = as.numeric(format(date, "%Y")),
         sp.leaf.fall.year = ifelse(doy < 150, year-1, year)) %>%
  group_by(sp, sp.leaf.fall.year) %>%
  mutate(annual.leaf.fall = sum(leaf_gm.int.raw, na.rm = TRUE),
         leaf_gm.int = leaf_gm.int.raw/annual.leaf.fall,
         # Assuming full leaf cover from June through October
         leaf_gm.int.mod = ifelse(doy >= 150 & doy < 300, 0, leaf_gm.int),
         leaf_gm.cum = cumsum(leaf_gm.int.mod),
         leaf_cover = 1 - leaf_gm.cum) %>%
  ungroup(sp, sp.leaf.fall.year) %>%
  mutate(leaf_gm.int.mov = rollmean(leaf_gm.int, k = set.k, fill = NA)) %>%
  group_by(sp, doy) %>%
  mutate(leaf_gm.int.mean = mean(leaf_gm.int, na.rm = TRUE),
         leaf_gm.int.sd = sd(leaf_gm.int, na.rm = TRUE),
         leaf_cover.mean = mean(leaf_cover, na.rm = TRUE),
         leaf_cover.sd = sd(leaf_cover, na.rm = TRUE)) %>%
  ungroup(sp, doy) %>%
  subset(sp != "na") %>%
  left_join(deci %>% dplyr::select(-deciduousness.label, -sp4), by = "sp") %>%
  mutate(sp.year = paste(sp, year, sep = "."))

sp.leaf_cover <- leaf.fall.int %>%
  select(sp, date, doy, year, leaf_cover, leaf_cover.mean, leaf_cover.sd)
sp.leaf_cover.mean <- leaf.fall.int %>%
  group_by(sp, doy) %>%
  summarise(leaf_cover.mean = mean(leaf_cover, na.rm = TRUE),
            leaf_cover.sd = sd(leaf_cover, na.rm = TRUE)) %>%
  ungroup(sp, doy)

save(sp.leaf_cover, file = file.path(results.folder, "sp.leaf_cover.Rdata"))
save(sp.leaf_cover.mean, file = file.path(results.folder, "sp.leaf_cover.mean.Rdata"))

f2 <- ggplot(leaf.fall.int, aes(x = date, y = leaf_gm.int)) +
  geom_line(aes(group = sp, color = sp), show.legend = FALSE)
ggsave(("leaf.fall.int_time_series_BCI.jpeg"),
       plot = f2, file.path(figures.folder.phen), device = "jpeg", height = 4.5, width = 9, units='in')

f2.1 <- ggplot(leaf.fall.int %>% subset(deciduousness != "Evergreen"), # sp %in% unique(iso.1.3.join$sp)
               aes(x = doy, y = leaf_gm.int)) +
  facet_wrap(sp ~ ., scales = "free_y") +
  geom_hline(yintercept = 0) +
  geom_line(aes(group = sp.year, color = deciduousness), size = 0.3) +
  # geom_ribbon(aes(ymin=leaf_gm.int.mean + leaf_gm.int.sd, ymax=leaf_gm.int.mean - leaf_gm.int.sd),
  #             fill='pink', alpha=0.8) +
  geom_line(aes(y = leaf_gm.int.mean), color = "blue") +
  ylab(expression('Leaf Fall/Annual Total')) + xlab("DOY") +
  geom_vline(aes(xintercept = 150), color = "red") +
  guides(color = guide_legend(order = 1, title = NULL, direction = "horizontal",
                              override.aes = list(size = 3))) +
  theme(legend.position = "top", legend.title = element_blank()) +
  scale_color_viridis_d(drop = FALSE) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 1, hjust = 1))
ggsave("leaf.fall.seasonality_BCI.jpeg",
       plot = f2.1, file.path(figures.folder.phen), device = "jpeg", height = 6, width = 10, units='in')

f2.2 <- f2.1 %+% subset(leaf.fall.int, deciduousness == "Evergreen")
ggsave("leaf.fall.seasonality_evergreens_BCI.jpeg",
       plot = f2.2, file.path(figures.folder.phen), device = "jpeg", height = 12, width = 18, units='in')

f2.3 <- ggplot(leaf.fall.int %>% subset(deciduousness != "Evergreen"), # sp %in% unique(iso.1.3.join$sp)
               aes(x = doy, y = leaf_cover)) +
  facet_wrap(sp ~ ., scales = "free_y") +
  geom_hline(yintercept = 0) +
  geom_line(aes(group = sp.year, color = deciduousness), size = 0.3) +
  # geom_ribbon(aes(ymin=leaf_gm.int.mean + leaf_gm.int.sd, ymax=leaf_gm.int.mean - leaf_gm.int.sd),
  #             fill='pink', alpha=0.8) +
  geom_line(aes(y = leaf_cover.mean), color = "blue") +
  ylab(expression('Leaf Cover Fraction')) + xlab("DOY") +
  # geom_vline(aes(xintercept = 120), color = "red") +
  guides(color = guide_legend(order = 1, title = NULL, direction = "horizontal",
                              override.aes = list(size = 3))) +
  theme(legend.position = "top", legend.title = element_blank()) +
  scale_color_viridis_d(drop = FALSE) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 1, hjust = 1))
ggsave(("leaf.cover_BCI.jpeg"),
       plot = f2.3, file.path(figures.folder.phen), device = "jpeg", height = 6, width = 10, units='in')

f2.4 <- f2.3 %+% subset(leaf.fall.int, deciduousness == "Evergreen")
ggsave("leaf.cover_evergreens_BCI.jpeg",
       plot = f2.4, file.path(figures.folder.phen), device = "jpeg", height = 12, width = 18, units='in')

f3 <- ggplot(leaf.fall.int,
             aes(x = doy, y = leaf_cover.mean)) +
  geom_line(aes(group = sp, color = deciduousness), size = 0.5) +
  theme(legend.position = c(0.7, 0.6), legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent")) +
  ylab("Leaf Cover Fraction") + xlab("DOY")
ggsave(("leaf.cover_BCI_single_panel.jpeg"),
       plot = f3, file.path(figures.folder.phen), device = "jpeg", height = 3, width = 4.5, units='in')

## Substituting 100% Leaf cover through the year for Evergreens
sp.leaf_cover.for.model <- sp.leaf_cover.mean %>%
  left_join(deci %>% select(sp, deciduousness), by = "sp") %>%
  mutate(leaf_cover = ifelse(deciduousness == "Evergreen", 1, leaf_cover.mean)) %>%
  select(-leaf_cover.sd) %>%
  as.data.frame()
f4 <- ggplot(sp.leaf_cover.for.model,
             aes(x = doy, y = leaf_cover)) +
  geom_line(aes(group = sp, color = deciduousness), size = 0.5) +
  theme(legend.position = c(0.7, 0.6), legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent")) +
  ylab("Leaf Cover Fraction") + xlab("DOY")
ggsave(("leaf.cover_BCI_single_panel_evergreensAt100.jpeg"),
       plot = f4, file.path(figures.folder.phen), device = "jpeg", height = 3, width = 4.5, units='in')

save(sp.leaf_cover.for.model, file = file.path(results.folder, "sp.leaf_cover.for.model.Rdata"))

#******************************************************
### Load Leaf Cohort tracking data from the crane sites------
#******************************************************
figures.folder.cohort <- paste0("figures/PhenoDemoTraitsPsi/leaf_cohort")
if(!dir.exists(file.path(figures.folder.cohort))) {dir.create(file.path(figures.folder.cohort))}

rama95 <- read_excel(file.path("data-raw/traits/Panama_cranes_leaf_phenology_data/RAMA95.xlsx")) %>%
  mutate(site = "PNM")
rama97 <- read_excel(file.path("data-raw/traits/Panama_cranes_leaf_phenology_data/RAMA97.xlsx")) %>%
  mutate(site = "PNM")
sherman <- read_excel(file.path("data-raw/traits/Panama_cranes_leaf_phenology_data/SHERMAN.xlsx")) %>%
  mutate(site = "SNL")

cohort <- rbind(rama95, rama97) %>%
  rbind(sherman) %>%
  rename_all(tolower) %>%
  rename(sp4 = sp) %>%
  mutate(dead = as.Date(dead), born = as.Date(born),
         born.doy = as.numeric(format(born, "%j")),
         dead.doy = as.numeric(format(dead, "%j"))) %>%
  left_join(deci %>% dplyr::select(sp, sp4, deciduous, deciduousness), by = "sp4") %>%
  unite("deci_sp", deciduous, sp, remove = FALSE) %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(deciduousness)]), ordered=TRUE),
         deci_sp.plot = factor(deci_sp, levels=unique(deci_sp[order(deciduousness)]), ordered=TRUE)) %>%
  subset(!is.na(deciduousness)) # only one species does not have a deciduousness label

coh.plot.base <- ggplot(cohort) +
  guides(color = guide_legend(title = "Deciduousness"))

coh.plot1 <- coh.plot.base +
  geom_point(aes(x = born, y = deci_sp.plot, color = deciduousness), alpha = 0.5) +
  theme(axis.text.y = element_text(size = 3))
ggsave(("born.doy_by_sp.jpeg"),
       plot = coh.plot1, file.path(figures.folder.cohort), device = "jpeg", height = 4.5, width = 9, units='in')

coh.plot2 <- coh.plot.base +
  geom_density(aes(x = born.doy, color = deciduousness, linetype = "Leaf Born"), size = 1) +
  geom_density(aes(x = dead.doy, color = deciduousness, linetype = "Leaf Dead"), size = 1) +
  scale_linetype_manual(name = "Event", values = c("solid", "twodash")) +
  xlab("DOY")
ggsave(("born.dead.doy_density_by_deciduousness.jpeg"),
       plot = coh.plot2, file.path(figures.folder.cohort), device = "jpeg", height = 4.5, width = 9, units='in')

rectangles.1 <- data.frame(
  xmin = 120,
  xmax = 335,
  ymin = 0,
  ymax = 0.05
)
coh.plot3.base <- coh.plot.base +
  scale_linetype_manual(name = "Event", values = c("solid", "twodash")) +
  facet_grid(site ~ deciduousness) +
  guides(color = NULL) +
  xlab("DOY") + ylab("Density") +
  geom_rect(data=rectangles.1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill='gray80', alpha=0.8)
coh.plot3.1 <- coh.plot3.base +
  geom_density(aes(x = born.doy, color = sp, linetype = "Leaf Born"), show.legend = FALSE) +
  ## limiting outliers
  scale_y_continuous(limits = c(0, 0.05)) + ggtitle("Timing of Leaf Births")
ggsave(("born.doy_density_by_sp_deciduousness.jpeg"),
       plot = coh.plot3.1, file.path(figures.folder.cohort), device = "jpeg", height = 4.5, width = 9, units='in')
coh.plot3.2 <- coh.plot3.base +
  geom_density(aes(x = dead.doy, color = sp, linetype = "Leaf Dead"), show.legend = FALSE) +
  ggtitle("Timing of Leaf Deaths")
ggsave(("dead.doy_density_by_sp_deciduousness.jpeg"),
       plot = coh.plot3.2, file.path(figures.folder.cohort), device = "jpeg", height = 4.5, width = 9, units='in')

coh.sp.summ <-  cohort %>%
  group_by(sp, site, deciduousness, deci_sp) %>%
  summarise(lifetime = mean(lifetime, na.rm = TRUE),
            born.mean = mean(born, na.rm = TRUE),
            dead.mean = mean(dead, na.rm = TRUE),
            born.sd = sd(born, na.rm = TRUE),
            dead.sd = sd(dead, na.rm = TRUE))
formula = y~x
coh.plot4.base <- ggplot(coh.sp.summ,
       aes(x = lifetime)) +
  guides(color = guide_legend(title = "Deciduousness")) +
  facet_wrap(site ~ .)
coh.plot4.1 <- ggplot(coh.sp.summ) +
  guides(color = guide_legend(title = "Deciduousness")) +
  facet_wrap(site ~ .) +
  geom_density(aes(x = born.mean, color = deciduousness)) +
  xlab("Mean of leaf born dates")
ggsave(("mean_leaf_born_dates.jpeg"),
         plot = coh.plot4.1, file.path(figures.folder.cohort), device = "jpeg", height = 4, width = 7, units='in')

coh.plot4.2 <- ggplot(coh.sp.summ) +
  guides(color = guide_legend(title = "Deciduousness")) +
  facet_wrap(site ~ .) +
  geom_density(aes(x = dead.mean, color = deciduousness)) +
  xlab("Mean of leaf death dates")
ggsave(("mean_leaf_death_dates.jpeg"),
       plot = coh.plot4.2, file.path(figures.folder.cohort), device = "jpeg", height = 4, width = 7, units='in')

coh.plot5 <- ggplot(coh.sp.summ,
                    aes(x = lifetime, y = born.sd)) +
  geom_point(aes(color = deciduousness), alpha = 0.8, size = 3) +
  guides(color = guide_legend(title = "Deciduousness")) +
  ylab("SD of Leaf Born Dates") + xlab("Leaf Longevity (Days)") +
  geom_smooth(data = subset(coh.sp.summ, site == "SNL"),
              method = "lm", se = FALSE) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.8, rr.digits = 2,
               formula = formula, parse = TRUE, size = 6) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.7, size = 6) +
  facet_wrap(site ~ .)
ggsave(("leaf_longevity_vs_leaf_born_days_concentration.jpeg"),
       plot = coh.plot5, file.path(figures.folder.cohort), device = "jpeg", height = 4, width = 7, units='in')

coh.plot6 <- ggplot(coh.sp.summ,
                    aes(x = lifetime, y = dead.sd)) +
  geom_point(aes(color = deciduousness), alpha = 0.8, size = 3) +
  guides(color = guide_legend(title = "Deciduousness")) +
  ylab("SD of Leaf Death Dates") + xlab("Leaf Longevity (Days)") +
  geom_smooth(data = subset(coh.sp.summ, site == "SNL"),
              method = "lm", se = FALSE) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.85, rr.digits = 2,
               formula = formula, parse = TRUE, size = 6) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.75, size = 6) +
  facet_wrap(site ~ .)
ggsave(("leaf_longevity_vs_leaf_death_days_concentration.jpeg"),
       plot = coh.plot6, file.path(figures.folder.cohort), device = "jpeg", height = 4, width = 7, units='in')

coh.summ <-  cohort %>%
  group_by(sp, site, deciduousness, deci_sp) %>%
  summarise(lifetime = mean(lifetime, na.rm = TRUE),
            born.doy.mean = mean(born.doy, na.rm = TRUE),
            dead.doy.mean = mean(dead.doy, na.rm = TRUE),
            born.doy.sd = sd(born.doy, na.rm = TRUE),
            dead.doy.sd = sd(dead.doy, na.rm = TRUE)) %>%
  ungroup(sp, site, deciduousness, deci_sp) %>%
  group_by(site, deciduousness) %>%
  summarise(lifetime = mean(lifetime, na.rm = TRUE),
            born.doy.sd = sd(born.doy.mean, na.rm = TRUE),
            dead.doy.sd = sd(dead.doy.mean, na.rm = TRUE),
            born.doy.mean = mean(born.doy.mean, na.rm = TRUE),
            dead.doy.mean = mean(dead.doy.mean, na.rm = TRUE))

save(cohort, file = file.path(results.folder, "cohort.Rdata"))
save(coh.sp.summ, file = file.path(results.folder, "coh.sp.summ.Rdata"))

