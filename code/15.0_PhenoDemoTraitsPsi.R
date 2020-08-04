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


rm(list=ls())
#*******************************************
####   Load Libraries, Prep for graphics, folders  ##
#*******************************************

if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, readxl, forcats, agricolae, gridExtra,
               scales, GGally, ggpmisc, Evapotranspiration,
               data.table, bci.elm.fates.hydro, mgcv)
# graphics info
theme_set(theme_bw())
theme_update(text = element_text(size = 14),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)

figures.folder <- paste0("figures/PhenoDemoTraitsPsi")
if(!dir.exists(file.path(figures.folder))) {dir.create(file.path(figures.folder))}
results.folder <- paste0("results/PhenoDemoTraitsPsi")
if(!dir.exists(file.path(results.folder))) {dir.create(file.path(results.folder))}
#****************************
###   Custom Functions   ####
#****************************

range01 <- function(x){(x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))}

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

psi.corr.fun.ls.2 <- list(
  "gr.Psi" =
    function(df, dflc) {
      result.df <-
        psi.study[, psi.mod := range01(Exponential(A = df$A, B = df$B, psi = -psi))][
            , keyby = .(depth, interval, par.sam), .(gfac = mean(psi.mod, na.rm = TRUE))][
              , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(list(result.df = result.df))
    },
  "gr.Psi.VPD" =
    function(df, dflc) {
      result.df <-
        psi.study[, psi.mod := range01(Exponential(A = df$A, B = df$B, psi = -psi))][
          , keyby = .(depth, interval, par.sam), .(gfac = mean(psi.mod*std.VPD, na.rm = TRUE))][
            , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(list(result.df = result.df))
    },
  "gr.Psi.leaf" =
    function(df, dflc) {
      # dflc.dt <- data.table(doy = dflc$doy, leaf_cover = dflc$leaf_cover)
      result.df <-
        as.data.table(psi.study)[data.table(dflc), on = 'doy'][,
          psi.mod := range01(Exponential(A = df$A, B = df$B, psi = -psi))][
            , keyby = .(depth, interval, par.sam), .(gfac = mean(psi.mod*leaf_cover, na.rm = TRUE))][
              , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(list(result.df = result.df))
    },
  "gr.Psi.VPD.leaf" =
    function(df, dflc) {
      dflc.dt <- data.table(doy = dflc$doy, leaf_cover = dflc$leaf_cover)
      result.df <-
        as.data.table(psi.study)[dflc.dt, on = 'doy'][,
        psi.mod := range01(Exponential(A = df$A, B = df$B, psi = -psi))][
          , keyby = .(depth, interval, par.sam), .(gfac = mean(c(psi.mod*std.VPD*leaf_cover), na.rm = TRUE))][
            , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(list(result.df = result.df))
  }
)

psi.corr.fun.ls <- list(
  "gr.Psi" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := Exponential(A = df$A, B = df$B, psi = -psi)][
            , keyby = .(depth, interval, par.sam), .(gfac = mean(psi.mod, na.rm = TRUE))][
              , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    },
  "gr.Psi.Rad" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := Exponential(A = df$A, B = df$B, psi = -psi)][
          , keyby = .(depth, interval, par.sam), .(gfac = mean(psi.mod + std.Rs, na.rm = TRUE))][
            , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    } ,
  "gr.Psi.VPD" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := Exponential(A = df$A, B = df$B, psi = -psi)][
          , keyby = .(depth, interval, par.sam), .(gfac = mean(psi.mod + std.VPD, na.rm = TRUE))][
            , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    } ,
  "gr.Psi.Rad.VPD" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := Exponential(A = df$A, B = df$B, psi = -psi)][
          , keyby = .(depth, interval, par.sam), .(gfac = mean(psi.mod + std.Rs.VPD, na.rm = TRUE))][
            , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    } ,
  "gr.Psi.Rad.PET" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := Exponential(A = df$A, B = df$B, psi = -psi)][
          , keyby = .(depth, interval, par.sam), .(gfac = mean(psi.mod + std.Rs.pet.PM, na.rm = TRUE))][
            , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    } ,
  "mr.Psi" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := psi - df$psi_kl80][
          psi.mod > 0, psi.mod := 0][
            , keyby = .(depth, interval, par.sam), .(gfac = sum(psi.mod, na.rm = TRUE))][
              , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    } ,
  "mr.Psi.PET" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := psi - df$psi_kl80][
          psi.mod > 0, psi.mod := 0][
            , keyby = .(depth, interval, par.sam), .(gfac = sum(psi.mod*std.pet.PM, na.rm = TRUE))][
              , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    } ,
  "mr.Psi.VPD" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := psi - df$psi_kl80][
          psi.mod > 0, psi.mod := 0][
            , keyby = .(depth, interval, par.sam), .(gfac = sum(psi.mod*std.VPD, na.rm = TRUE))][
              , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    } ,
  "mr.Psi.I" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := indicator(psi, df$psi_kl80, greater.than = FALSE)][
          , keyby = .(depth, interval, par.sam), .(gfac = sum(psi.mod, na.rm = TRUE))][
            , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    } ,
  "mr.Psi.PET.I" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := indicator(psi, df$psi_kl80, greater.than = FALSE)][
          , keyby = .(depth, interval, par.sam), .(gfac = sum(psi.mod*std.pet.PM, na.rm = TRUE))][
            , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    },
  "mr.Psi.VPD.I" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := indicator(psi, df$psi_kl80, greater.than = FALSE)][
          , keyby = .(depth, interval, par.sam), .(gfac = sum(psi.mod*std.VPD, na.rm = TRUE))][
            , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(result.df)
    }
)

get.mfac.ls <- list(
  "mr.kl80.I" =
  function(df) {
    result.df <-
      as.data.table(psi.study)[, psi.mod := indicator(psi, df$psi_kl80, greater.than = FALSE)][
        , keyby = .(depth, interval, par.sam), .(mfac = sum(psi.mod, na.rm = TRUE))][
          , keyby = .(depth, interval), .(mfac = mean(mfac, na.rm = TRUE))]
    result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "mfac")
    return(result.df)
  },
  "mr.kl50.I" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := indicator(psi, df$psi_kl50, greater.than = FALSE)][
          , keyby = .(depth, interval, par.sam), .(mfac = sum(psi.mod, na.rm = TRUE))][
            , keyby = .(depth, interval), .(mfac = mean(mfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "mfac")
      return(result.df)
    },
  "mr.kl80.I.VPD" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := indicator(psi, df$psi_kl80, greater.than = FALSE)][
          , keyby = .(depth, interval, par.sam), .(mfac = sum(psi.mod*std.VPD, na.rm = TRUE))][
            , keyby = .(depth, interval), .(mfac = mean(mfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "mfac")
      return(result.df)
    },
  "mr.kl50.I.VPD" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := indicator(psi, df$psi_kl50, greater.than = FALSE)][
          , keyby = .(depth, interval, par.sam), .(mfac = sum(psi.mod*std.VPD, na.rm = TRUE))][
            , keyby = .(depth, interval), .(mfac = mean(mfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "mfac")
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

get.ml.depth.rsq.top <- function(result) {
  r2.sorted <- sort(result$R2, decreasing = FALSE, index.return = TRUE)
  return(list(depth = result$depth[r2.sorted$ix[which(r2.sorted$x >= r2.sorted$x[1] - 0.1)]],
              R2 = r2.sorted$x[which(r2.sorted$x >= r2.sorted$x[1] - 0.1)]))
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
### Load data -------
#******************************************************
source("code/load.R")

#******************************************************
### Calculate Correlation of growth rates with psi by depth -------
#******************************************************
## 1. psi: interval mean psi by depth
## 2. psi.p50.g1 or g2: interval mean psi by depth: but only for the period when psi was above p50/p80
## 3. psi.p50.g1.n: number of days psi was above p50/p80 or group 1 p50 vs. group 2 p80

## Limiting to evergreen species for 2 & 3: For non-evergreen species periods in question will have to be a subset of days when leaves were on
# growth.selection <- "size_class_predefined_cc_scaled"
# intervals <- 5
# load(file = paste0("results/gro.long.cc.med_", intervals, "_", growth.selection, ".Rdata"))
#
# growth.sub <- gro.long.cc.med %>% subset(size == "large") %>%
#   rename(demo.rate = med.dbh.resid)

growth.sub <- lapply(growth[grep("large", names(growth))], as.data.frame) %>%
  bind_rows(.id = "sp_size") %>%
  rename(demo.rate = median) %>%
  # rename(demo.rate = dbh.residuals) %>%
  separate(sp_size, c("sp", "size", sep = "_"), remove = FALSE, extra = "drop", fill = "right") %>%
  dplyr::select(-sp_size, -"_") %>%
  # remove species taht are not canopy
  subset(sp %in% bci.traits$sp[bci.traits$form1 == "T"])
length(unique(growth.sub$sp))

grate.plot <- ggplot(growth.sub, aes(x = interval, y = demo.rate)) +
  geom_line(aes(group = sp, color = sp), show.legend = FALSE) +
  facet_wrap(. ~ sp) +
  ylab("Std. Growth") + xlab("Interval") +
  geom_errorbar(aes(ymax = upr, ymin = lwr), width = 0.1, size = 0.5)
ggsave(paste0("Std.Growth_", growth.type,".jpeg"),
       plot = grate.plot, file.path(figures.folder), device = "jpeg", height = 7, width = 7, units='in')

grate.plot.iso.sp <- grate.plot %+%
  subset(growth.sub, sp %in% iso.1.3.join$sp[iso.1.3.join$source == "Meinzer et al.1999 Fig. 4"])
ggsave(paste0("Std.Growth_iso.sp_", growth.type,".jpeg"),
       plot = grate.plot.iso.sp, file.path(figures.folder), device = "jpeg", height = 3, width = 4, units='in')

census.meds <- readr::read_rds("results/census.mediandates.rds")
census.beg <- census.meds[3: length(census.meds)]
cut.breaks <- census.beg
cut.breaks.2 <- as.Date(paste0(seq(1990, 2015, by = 5), "-01-01"))
cut.labels.2 <- paste0(seq(1990, 2010, by = 5), "-", seq(1995, 2015, by = 5))
cut.labels.interval <- 3: (length(census.meds)-1)

mort.sub <- mrate.long %>% subset(size == "large") %>%
  mutate(interval = as.numeric(recode(census, `1985` = "1",
                                      `1990` = "2", `1995` = "3",
                                      `2000` = "4", `2005` = "5",
                                      `2010` = "6", `2015` = "7"))) %>%
  subset(interval %in% cut.labels.interval) %>%
  dplyr::select(sp_size, interval, mrate, diff.mrate) %>%
  rename(demo.rate = mrate) %>%
  separate(sp_size, c("sp", "size", sep = "_"), remove = FALSE, extra = "drop", fill = "right") %>%
  dplyr::select(-sp_size, -"_")

## getting Kmax_by_PSI parameters for species with growth data or mortality
## Using k_leaf_by_LWP curve parameters directly fitted to data, rather than predicted from soft traits, when available

## combining psi, PET and VPD

clim.daily.effect <- clim.daily %>%
  mutate(Rs.pet.PM.effect = as.numeric(predict(gpp.models$eq.gpp.rad.pet.gam, newdata = clim.daily)),
         Rs.VPD.effect = as.numeric(predict(gpp.models$eq.gpp.rad.vpd.gam, newdata = clim.daily)),
         std.Rs.pet.PM = range01(Rs.pet.PM.effect),
         std.Rs.VPD = range01(Rs.VPD.effect),
         pet.PM.effect = as.numeric(predict(gpp.models$eq.gpp.pet, newdata = clim.daily)),
         VPD.effect = as.numeric(predict(gpp.models$eq.gpp.vpd, newdata = clim.daily)),
         Rs.effect = as.numeric(predict(gpp.models$eq.gpp.rad, newdata = clim.daily)),
         std.pet.PM = range01(pet.PM.effect),
         std.VPD = range01(VPD.effect),
         std.Rs = range01(Rs.effect))
#
# plot(pet.PM.effect ~ pet.PM, clim.daily.effect[1:1000,])
# plot(std.pet.PM ~ pet.PM, clim.daily.effect[1:1000,])
# jpeg(file.path(figures.folder, "pet_vs.std.pet.jpeg"), width = 4, height = 4, units = "in", pointsize = 10,
#      quality = 100, res = 300)
# par(mar = c(5,4,2,2))
# plot(std.pet.PM ~ pet.PM, clim.daily.effect[1:1000,])
# dev.off()
# plot(VPD.effect ~ VPD, clim.daily.effect[1:1000,])
# plot(std.VPD ~ VPD, clim.daily.effect[1:1000,])
# jpeg(file.path(figures.folder, "vpd_vs.std.vpd.jpeg"), width = 4, height = 4, units = "in", pointsize = 10,
#      quality = 100, res = 300)
# par(mar = c(5,4,2,2))
# plot(std.VPD ~ VPD, clim.daily.effect[1:1000,])
# dev.off()
# jpeg(file.path(figures.folder, "rad_vs.std.rad.jpeg"), width = 4, height = 4, units = "in", pointsize = 10,
#      quality = 100, res = 300)
# par(mar = c(5,4,2,2))
# plot(std.Rs ~ Rs, clim.daily.effect[1:1000,])
# dev.off()

# psi.m <- psi.depths %>%
psi.m <- psi %>%
  mutate(interval = cut(date, include.lowest = TRUE, breaks = cut.breaks,
                        labels = cut.labels.interval, right = TRUE)) %>%
  # full_join(clim.daily %>%
  #             mutate(std.pet.PM = 1 - range01(pet.PM),
  #                    std.VPD = 1 - range01(VPD),
  #                    std.Rs = range01(Rs)) %>%
  full_join(clim.daily.effect %>%
              dplyr::select(date, std.Rs.pet.PM, std.Rs.VPD,
                     std.Rs, std.pet.PM, std.VPD), by = "date")
depth.sub <- c(soil.depths[5:length(soil.depths)])
depth.breaks <- c(soil.depths[5], soil.depths[7:length(soil.depths)])
depth.labels <- c(0.5, soil.depths[8:length(soil.depths)])

psi.study <- as.data.table(psi.m)[!is.na(interval),][,
  ':='(doy = format(date, "%j"), year = format(date, "%Y"))][
  (doy < 368) & (!year %chin% c("1990")) &
    (!year %chin% c("1991") | depth == 2.9) &
    (!year %chin% c("1991", "1992") | depth < 2.9)][,
     c("year") := NULL][,
      doy := as.numeric(as.character(doy))][depth %chin% depth.sub][, ':='(depth = cut(depth, include.lowest = TRUE, breaks = depth.breaks,
                                                                                                   labels = depth.labels, right = TRUE))][, by = .(date, interval, interval.yrs, par.sam, depth),
                                                                                                                                          lapply(.SD, median, na.rm = TRUE)]

# psi.study <- as.data.table(psi.m)[!is.na(interval),][,
#   ':='(doy = format(date, "%j"),
#        year = format(date, "%Y"))][
#   (doy < 368) & (!year %chin% c("1990")) &
#     (!year %chin% c("1991") | depth == 3) &
#     (!year %chin% c("1991", "1992") | depth < 3)][, c("year") := NULL][,
#     doy := as.numeric(as.character(doy))]

# psi.study <- as.data.frame(psi.study)

tlp.sp <-  traits.long %>% subset(trait == "TLP") %>% select(sp, value) %>% rename(tlp = value)
n.tlp.sp <- length(tlp.sp$sp)
tlp.sp <- tlp.sp %>%
  mutate(psi_threshold = tlp*0.8)
tlp.sp.ls <- split(tlp.sp, f = list(tlp.sp$sp), drop = TRUE)

data.model.AB.sub <- data.model.AB %>%
  ## Using parameters fitted to data when available
  mutate(A = ifelse(is.na(data.A), model.A, data.A),
         B = ifelse(is.na(data.B), model.B, data.B)) %>%
  subset(sp %in% union(growth.sub$sp, mort.sub$sp)) %>%
  # ## And LMA from LMALAM from LEAF and DISC
  subset(sp.LMA.sub == "original") %>%
  ## removing A & B predicted from gap-filled WSG
  subset(sp.WSG.sub == "original")

setdiff(unique(growth.sub$sp), unique(data.model.AB.sub$sp))

# sp not present in sp.leaf_cover.mean but in growth data sets
# sp.growth.not.leaf_cover <- setdiff(unique(growth.sub$sp), unique(sp.leaf_cover.mean$sp))
# sp.leaf_cover.for.model <- sp.leaf_cover.mean %>%
#   bind_rows(data.frame(sp = rep(sp.growth.not.leaf_cover, each = max(sp.leaf_cover.mean$doy)),
#                        doy = rep(unique(sp.leaf_cover.mean$doy), times = length(sp.growth.not.leaf_cover)),
#                        leaf_cover.mean = NA, leaf_cover.sd = NA)) %>%
#   left_join(deci %>% select(sp, deciduous), by = "sp") %>%
#   mutate(leaf_cover = ifelse(deciduous == "E", 1, leaf_cover.mean)) %>%
#   select(-leaf_cover.sd) %>%
#   as.data.frame()
  ## But "ficutr", "pri2co", "termob" do not have leaf_cover data and also not Evergreen either but DB,
  ## so taking average DB species' leaf-cover,  "pri2co" 's deciduousness not known, so can't gap-fill
#   # First joinging average deci leaf_cover then substituing for missing data by species
# sp.leaf_cover.for.model <- sp.leaf_cover.for.model %>%
#   left_join(sp.leaf_cover.for.model %>% group_by(deciduous, doy) %>%
#               summarise(deci.leaf_cover = mean(leaf_cover, na.rm = TRUE), .groups = "drop"),
#             by = c("deciduous", "doy")) %>%
#   mutate(leaf_cover = ifelse(is.na(leaf_cover), deci.leaf_cover, leaf_cover))
#
# setdiff(unique(growth.sub$sp), unique(sp.leaf_cover.for.model$sp))
## unique(sp.leaf_cover.for.model[is.na()]$sp)

AB.sp.ls <- split(data.model.AB.sub, f = list(data.model.AB.sub$sp), drop = TRUE)
leaf_cover.sp.ls <- split(sp.leaf_cover.for.model, f = list(sp.leaf_cover.for.model$sp), drop = TRUE)

sp.ab.leaf_cover <- intersect(unique(data.model.AB.sub$sp), unique(sp.leaf_cover.for.model$sp))

## both needs to have the same species (in the same order)
AB.sp.ls <- AB.sp.ls[sp.ab.leaf_cover]
leaf_cover.sp.ls <- leaf_cover.sp.ls[sp.ab.leaf_cover]

# for (i in 1:length(names.gfac)) {
#   gfac.interval[[i]] <- lapply(lapply(AB.sp.ls, psi.corr.fun.ls[[i]]),
#                               as.data.frame) %>%
#     bind_rows(.id = "sp")
# }

## for species with traits data
# growth.sub <- growth[names(growth) %in% paste0(unique(c(hyd$sp, traits$sp)), "_large")]
names.gfac <- names(psi.corr.fun.ls.2) #c("g.Psi", "g.Psi.Rad.VPD", "g.Psi.Rad.PET")
## Preparing PSI matrices to compare against
gfac.interval <- vector(mode = "list", length = length(names.gfac))
names(gfac.interval) <- names.gfac  # "psi.p50.g1", "psi.p50.g2"
for (i in 1:length(names.gfac)) {
  gfac.interval[[i]] <-  lapply(mapply(FUN = psi.corr.fun.ls.2[[i]],
                                      AB.sp.ls, leaf_cover.sp.ls),
                                as.data.frame) %>%
    bind_rows(.id = "sp") %>% separate(sp, c("sp", NA, NA))
}

ml.ls <- vector("list", length(gfac.interval))
ml.dens <- vector("list", length(gfac.interval))
ml.rsq.ls <- vector("list", length(gfac.interval))
ml.corr.ls <- vector("list", length(gfac.interval))

for (i in names(gfac.interval)) {
  if(grepl("gr.", i)) {
    demo <- growth.sub
  } else {
    demo <- mort.sub
  }
  demo.psi <- demo %>%
    mutate(interval = as.factor(interval)) %>%
    inner_join(gfac.interval[[i]], by = c("interval", "sp")) %>%
    subset(!is.na(demo.rate) & is.finite(demo.rate)) %>% droplevels()
  demo.psi.ls <- split(demo.psi,
                       f = list(demo.psi$sp), drop = TRUE)
  ml.ls[[i]] <- lapply(demo.psi.ls, get.ts.lk)
  ml.dens[[i]] <- sapply(ml.ls[[i]], get.ml.max)
  ml.corr.ls[[i]] <- do.call(rbind.data.frame, ml.ls[[i]])
  ml.rsq.ls[[i]] <- do.call(rbind, lapply(ml.ls[[i]], get.ml.depth.rsq))
}

# for(n in 4:length(ml.dens)) {
#   if(n == 1) plot(density(ml.dens[[n]][1:length(ml.dens[[n]])]))
#   lines(density(ml.dens[[n]]), col = terrain.colors(length(ml.dens))[n])
#   abline(v = median(ml.dens[[n]]))
# }

ml.corr <- vector("list", length(gfac.interval))
ml.corr.best <- vector("list", length(gfac.interval))
for (i in names(gfac.interval)) {
  ml.corr[[i]] <- ml.corr.ls[[i]] %>%
    mutate(sp.depth = row.names(ml.corr.ls[[i]])) %>%
    separate(sp.depth, c("sp", "depth.fac", sep = "."), remove = FALSE, extra = "drop", fill = "right") %>%
    dplyr::select(-".", -depth.fac) %>%
    mutate(size = "large",
           depth = as.numeric(as.character(depth))) %>%
    arrange(sp, depth) %>%
    left_join(deci %>% dplyr::select(-sp4), by = "sp") %>%
    droplevels() %>%
    transform(deciduousness = factor(deciduousness,
                                     levels = c("Evergreen", "Brevideciduous",
                                                "Facultative Deciduous", "Obligate Deciduous"), ordered = TRUE)) %>%
    unite("deci_sp", deciduous, sp, remove = FALSE) %>%
    mutate(sp.plot = factor(sp, levels=unique(sp[order(deciduousness)]), ordered=TRUE),
           deci_sp.plot = factor(deci_sp, levels=unique(deci_sp[order(deciduousness)]), ordered=TRUE))

  # left_join(traits.long.hyd %>% select(deci_sp, deci_sp.plot, sp, sp.plot,  deciduousness), by = "sp") #%>%
  ml.corr.best[[i]] <- ml.corr[[i]] %>% group_by(sp) %>%
    mutate(corr.max = max(corr, na.rm = TRUE)) %>%
    subset(corr == corr.max) %>% dplyr::select(-corr.max) %>%
    ungroup(sp)
}
## species by depth by heat rsq
ml.rsq.combine <- dplyr::bind_rows(ml.corr, .id = "corr.func") %>%
  transform(corr.func = factor(corr.func, levels = names.gfac))
ml.rsq.combine.best <- dplyr::bind_rows(ml.corr.best, .id = "corr.func") %>%
  transform(corr.func = factor(corr.func, levels = names.gfac)) %>%
  unite(corr.func_sp_depth, corr.func,sp, depth, remove = FALSE)


## Plot chosen ERD

df.erd.to.plot <- ml.rsq.combine.best %>%
  subset(corr.func == "gr.Psi.VPD") %>%
  mutate(size = as.character(size)) %>%
  left_join(bci.traits %>% dplyr::select(sp, form1), by = "sp") %>%
  subset(form1 == "T" &
           corr >= 0 & R2 >= 0.1 & !duplicated(sp) & !is.na(depth)) %>% droplevels() %>%
  transform(sp = reorder(sp, depth))
length(unique(df.erd.to.plot$sp))
# 32
erd.sp.plot <- ggplot(df.erd.to.plot,
       aes(x = sp, y = depth)) +
  geom_point(aes(color = sp), show.legend = FALSE) +
  ylab("Effective Rooting Depth (m)") + xlab("Species") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(trans = reverselog_trans(10), breaks = ml.rsq.combine.best$depth)
ggsave("ERD_by_sp_large_canopy.jpeg",
       plot = erd.sp.plot, file.path(figures.folder), device = "jpeg", height = 3.5, width = 5, units='in')

grate.gfac <- gfac.interval.long <- vector(mode = "list", length = length(names.gfac))
names(gfac.interval.long) <- names(grate.gfac) <- names.gfac
for (i in 1:length(names.gfac)) {
  gfac.interval.long[[i]] <- gfac.interval[[i]] %>%
    pivot_longer(cols = c(-sp, -interval),
                 names_to = "depth", values_to = "gfac") %>%
    rename(interval.num =  interval) %>%
    mutate(depth = as.numeric(depth),
           interval.num =  as.numeric(as.character(interval.num)),
           size = "large")
  grate.gfac[[i]] <- growth.sub %>% rename(interval.num = interval) %>%
    left_join(gfac.interval.long[[i]], by = c("interval.num", "sp", "size")) %>%
    subset(!is.na(gfac))
}

grate.gfac.best <- dplyr::bind_rows(grate.gfac, .id = "corr.func") %>%
  mutate(censusint.m = recode(interval.num, `1` = "1982-85", `2` = "1985-90",
                              `3` = "1990-95", `4` = "1995-00", `5` = "2000-05", `6` = "2005-10", `7` = "2010-15")) %>%
  transform(corr.func = factor(corr.func, levels = names.gfac)) %>%
  unite(corr.func_sp_depth, corr.func,sp, depth, remove = FALSE) %>%
  left_join(bci.traits %>% dplyr::select(sp, form1), by = "sp") %>%
  subset(form1 = "T") %>%
  group_by(sp, size, corr.func, depth) %>%
  mutate(std.gfac = scale(gfac),
         std.growth = scale(demo.rate)) %>%
  left_join(ml.rsq.combine.best %>%
            dplyr::select(sp, size, corr.func, R2),
            by = c("sp", "size", "corr.func"))
grate.gfac.best.sub <- grate.gfac.best %>%
  subset(#corr.func %in% names.gfac &
           sp %in% iso.1.3.join$sp[iso.1.3.join$source == "Meinzer et al.1999 Fig. 4"] &
           ## those that are likely leafless in Mar-Apr
           !sp %in% as.character(leafless_mar.apr$sp[leafless_mar.apr$leafless_in_mar_apr_from_notes == "Yes"]) &
           corr.func %in% unique(grate.gfac.best$corr.func)[grep("gr.", unique(grate.gfac.best$corr.func))])
grate.gfac.best.plot <- ggplot(grate.gfac.best.sub %>%
           subset(corr.func_sp_depth %in% ml.rsq.combine.best$corr.func_sp_depth) %>% droplevels(),
       aes(x = censusint.m)) +
  geom_line(aes(y = std.growth, group = corr.func_sp_depth, color = sp, linetype = "Std.Growth"), size = 1) +
  geom_line(aes(y = std.gfac, group = corr.func_sp_depth, color = sp, linetype = "Std.Growth Factor"), size = 1) +
  facet_grid(corr.func ~ sp , scales = "free_y") +
  scale_color_brewer(palette = "Dark2") +
  xlab("Census Interval") + ylab("Std.Growth/Growth Factor") +
  guides(color = "none", linetype = guide_legend(order = 1, title = NULL, direction = "horizontal", label.position = "top", override.aes =
                                   list(linetype = c("Std.Growth" = "solid", "Std.Growth Factor" = "dotted")))) +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_text(aes(label = round(R2, 1), x = "2005-10", y = 1.2))
ggsave("Std.Growth Vs. Std.Growth Factor.jpeg",
       plot = grate.gfac.best.plot, file.path(figures.folder), device = "jpeg", height = 7, width = 14, units='in')

grate.gfac.best.plot.depths <- ggplot(grate.gfac.best.sub, #%>%
                                      # subset(depth %in% c(0.1, 1.0, 2.0)) %>% droplevels(),
                               aes(x = censusint.m)) +
  geom_line(aes(y = std.gfac, group = corr.func_sp_depth, color = depth), size = 1) +
  geom_line(aes(y = std.growth, group = corr.func_sp_depth,
                linetype = "Std.Growth"), size = 1, color = "red") +
  facet_grid(corr.func ~ sp , scales = "free_y") +
  # scale_color_discrete(name = "Depth (m)", drop = FALSE) +
  scale_color_continuous(name = "Depth (m)", trans = "reverse") +
  xlab("Census Interval") + ylab("Std.Growth/Growth Factor") +
  scale_linetype_manual(name = "", values = c("solid")) +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_text(aes(label = round(R2, 1), x = "2000-05", y = 1.2))
ggsave("Std.Growth Vs. Std.Growth Factor_depths.jpeg",
       plot = grate.gfac.best.plot.depths, file.path(figures.folder), device = "jpeg", height = 7, width = 14, units='in')

grate.gfac.best.plot.depths.best <- ggplot(grate.gfac.best.sub %>%
           subset(corr.func_sp_depth %in% ml.rsq.combine.best$corr.func_sp_depth),
           aes(x = censusint.m)) +
  geom_line(aes(y = std.gfac, group = corr.func_sp_depth, color = as.factor(depth)), size = 1) +
  geom_line(aes(y = std.growth, group = corr.func_sp_depth,
                linetype = "Std.Growth"), size = 1) +
  facet_grid(corr.func ~ sp , scales = "free_y") +
  scale_color_discrete(name = "Depth (m)", drop = FALSE) +
  # scale_color_continuous(name = "Depth (m)", trans = "reverse") +
  xlab("Census Interval") + ylab("Std.Growth/Growth Factor") +
  scale_linetype_manual(name = "", values = c("solid")) +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_text(aes(label = round(R2, 1), x = "2000-05", y = 1.2))
ggsave("Std.Growth Vs. Std.Growth Factor_depths.best.jpeg",
       plot = grate.gfac.best.plot.depths.best, file.path(figures.folder), device = "jpeg", height = 7, width = 14, units='in')

iso.2 <- iso.2.raw %>%
  # subset(source == "Meinzer et al.1999 Fig. 5A") %>%
  mutate(source = "Meinzer et al.1999 Fig. 5B", location = "BCI") %>%
  group_by(sp, source) %>%
  summarise(se = sd(Xylem_sap_deltaD_permil, na.rm = TRUE)/sqrt(n()),
            Xylem_sap_deltaD_permil = mean(Xylem_sap_deltaD_permil, na.rm = TRUE),
            n = n(),
            DBH = mean(DBH, na.rm = TRUE))

## those species that were likely leafless at the time of Xylem sap isotopes collection
## in Mar & April need to be removed
# pse1se zuelgu sponra huracr pla2el

ml.rsq.combine.best <- ml.rsq.combine.best %>%
  # left_join(iso.2 %>%
  left_join(iso.1.3.join %>%
              subset(!sp %in% as.character(leafless_mar.apr$sp[leafless_mar.apr$leafless_in_mar_apr_from_notes == "Yes"])) %>%
              dplyr::select(sp, Xylem_sap_deltaD_permil, se, source), by = "sp") %>%
  droplevels()
ml.rsq.combine.best <- ml.rsq.combine.best %>%
  # left_join(iso.2 %>%
  left_join(iso.1.3.join %>%
              subset(!sp %in% as.character(leafless_mar.apr$sp[leafless_mar.apr$leafless_in_mar_apr_from_notes == "Yes"])) %>%
              group_by(sp) %>%
              summarise(Xylem_sap_deltaD_permil.mean = mean(Xylem_sap_deltaD_permil, na.rm = TRUE),
                        se.mean = mean(se, na.rm = TRUE)) %>%
              dplyr::select(sp, Xylem_sap_deltaD_permil.mean, se.mean), by = "sp") %>%
  droplevels()

depth.rsq.isotopes <- ml.rsq.combine.best %>%
  group_by(corr.func, sp, size) %>%
  summarise_at(c("depth", "Xylem_sap_deltaD_permil", "se"), mean, na.rm = TRUE) %>%
  ungroup(corr.func, sp, size)
save(depth.rsq.isotopes, file = file.path(results.folder, "depth.rsq.isotopes.Rdata"))
save(ml.rsq.combine.best, file = file.path(results.folder, "ml.rsq.combine.best.Rdata"))
save(ml.rsq.combine, file = file.path(results.folder, "ml.rsq.combine.Rdata"))
load(file = file.path(results.folder, "depth.rsq.isotopes.Rdata"))
load(file = file.path(results.folder, "ml.rsq.combine.best.Rdata"))
load(file = file.path(results.folder, "ml.rsq.combine.Rdata"))

### Plot best correlated depth against isotopic data and traits-----

# ml.rsq.combine.best <- ml.rsq.combine.best %>%
#   left_join(traits.wide.hyd %>% select(sp, KmaxL, Panama.moist.pref:HSMLWP.TLP), by = "sp")

heat.rsq <- ggplot(ml.rsq.combine %>% subset(corr >= 0) %>% droplevels(),
                   aes(y = deci_sp, x = as.factor(depth))) +
  geom_tile(aes(fill = corr)) +
  ylab("Species") + xlab("Depth (m)") +
  facet_wrap(. ~ corr.func, nrow = 1) +
  scale_fill_viridis_c(expression("Pearson's "*italic(rho)), direction = -1, option = "plasma") #+
ggsave("psi.corr_all.depths_phenology_heat_by_corr.func.jpeg",
       plot = heat.rsq, file.path(figures.folder), device = "jpeg", height = 5, width = 12, units='in')

# theme(axis.text.y = element_text(angle = 90, vjust = 0.5)) +
# scale_x_continuous(breaks = soil.depths[c(1,8:13)])
heat.best.rsq <- heat.rsq %+% subset(ml.rsq.combine.best, R2 >= 0.2)
ggsave("psi.corr_best.depth_phenology_heat_by_corr.func.jpeg",
       plot = heat.best.rsq, file.path(figures.folder), device = "jpeg", height = 5, width = 12, units='in')

xylem.label <- expression('Xylem Sap '*delta~""^2*"H (\u2030)"*'')
ml.rsq.combine.best <- ml.rsq.combine.best %>% left_join(bci.traits %>% dplyr::select(sp, form1), by = "sp") %>%
  mutate(depth = as.numeric(depth))
ml.rsq.combine.sub <- ml.rsq.combine.best %>%
  subset(corr >= 0 & form1 == "T" &
           corr.func %in% names.gfac[1:5] &
           !sp %in% c("guapst") &
           # !sp %in% c("guapst", "alsebl") &
           R2 >= 0.1 & #corr.func == "gr.Psi.Rad.VPD" &
           !is.na(Xylem_sap_deltaD_permil.mean)) %>%
  droplevels()

formula = y~x
p0 <- ggplot(ml.rsq.combine.sub,
             aes(x = Xylem_sap_deltaD_permil.mean, y = depth)) + #HSMTLP.80L)) +
  geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil.mean + se.mean,
                     xmin = Xylem_sap_deltaD_permil.mean - se.mean, color = deciduousness),
                 size = 0.5, height = 0.05) +
  facet_wrap( ~ corr.func, nrow = 1) +
  geom_text(aes(x =  Xylem_sap_deltaD_permil.mean, y = depth, label = sp, color = deciduousness), nudge_y = 0.1, nudge_x = 0.2,
            size = 4, show.legend = FALSE) +
  ylab(expression("Best Correlated Depth (m)")) + xlab(xylem.label) +
  scale_y_continuous(trans = reverselog_trans(10), breaks = ml.rsq.combine.sub$depth) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.2, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.1, size = 4) +
  geom_point(size = 1, show.legend = TRUE, aes(color = deciduousness)) +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.direction = "horizontal") + #, strip.text = element_blank()) +
  scale_color_brewer(palette = "Dark2")
ggsave("psi.corr_best.depth_xylem_sap_deltaD_phenology_mean_isotope_source.jpeg",
       plot = p0, file.path(figures.folder), device = "jpeg", height = 5, width = 8, units = 'in')

p1 <- ggplot(ml.rsq.combine.sub,
             aes(x = Xylem_sap_deltaD_permil, y = depth)) + #HSMTLP.80L)) +
  geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + se,
                     xmin = Xylem_sap_deltaD_permil - se, color = deciduousness),
                 size = 0.5, height = 0.05) +
  facet_wrap( ~ corr.func, nrow = 1) +
  geom_text(aes(x =  Xylem_sap_deltaD_permil, y = depth, label = sp, color = deciduousness), nudge_y = 0.1, nudge_x = 0.2,
            size = 4, show.legend = FALSE) +
  ylab(expression("Best Correlated Depth (m)")) + xlab(xylem.label) +
  scale_y_continuous(trans=reverselog_trans(10), breaks = ml.rsq.combine.sub$depth) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.9, npcy = 0.2, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.9, npcy = 0.1, size = 4) +
  geom_point(size = 2, show.legend = TRUE, aes(color = deciduousness)) +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.direction = "horizontal") + #, strip.text = element_blank()) +
  scale_color_brewer(palette = "Dark2")
ggsave("psi.corr_best.depth_xylem_sap_deltaD_phenology_two_isotope_sources.jpeg",
       plot = p1, file.path(figures.folder), device = "jpeg", height = 5, width = 8, units = 'in')
p2 <- p1 %+% subset(ml.rsq.combine.sub, source == "Meinzer et al.1999 Fig. 4")
ggsave("psi.corr_best.depth_xylem_sap_deltaD_phenology_Meinzer.jpeg",
       plot = p2, file.path(figures.folder), device = "jpeg", height = 5, width = 8, units = 'in')

p3 <- ggplot(ml.rsq.combine.sub %>% subset(source == "Meinzer et al.1999 Fig. 4"),
       aes(x = Xylem_sap_deltaD_permil, y = depth)) + #HSMTLP.80L)) +
  coord_cartesian(ylim = c(13, 0.3)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", size = 0.5) +
  geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + se,
                     xmin = Xylem_sap_deltaD_permil - se, color = deciduousness),
                 size = 0.5, height = 0.3) +
  facet_wrap( ~ corr.func, nrow = 1) +
  geom_text(aes(x =  Xylem_sap_deltaD_permil, y = depth, label = sp, color = deciduousness), nudge_y = 0.1, nudge_x = 0.2,
            size = 4, show.legend = FALSE) +
  # position=position_jitter(width=ifelse(ml.rsq.combine.sub$sp=='cordal',1,0),
  #                        height=ifelse(ml.rsq.combine.sub$sp=='cordal',1,0))
  ylab(expression("Water Uptake Depth (m)")) + xlab(xylem.label) +
  scale_y_continuous(trans=reverselog_trans(10), breaks = ml.rsq.combine.sub$depth) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.9, npcy = 0.2, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.9, npcy = 0.1, size = 4) +
  geom_point(size = 2, show.legend = TRUE, aes(color = deciduousness)) +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.direction = "horizontal") +
  scale_color_brewer(palette = "Dark2")
ggsave("psi.corr_best.depth_xylem_sap_deltaD_phenology_Meinzer.jpeg",
       plot = p3, file.path(figures.folder), device = "jpeg", height = 4, width = 8, units = 'in')


ml.rsq.combine.sub <- ml.rsq.combine.sub %>% transform(models.plot1 = factor(corr.func,
                                                                            labels = c("A", "B", "C", "D")),
                                                       models.plot2 = factor(corr.func,
                                                                            labels = c(expression(italic(K[italic('Leaf')])),
                                                                                                  expression(italic(K[italic('Leaf')]*'VPD')),
                                                                                                             expression(italic(K[italic('Leaf')]*'LeafFrac')),
                                                                                                                        expression(italic(K[italic('Leaf')]*'VPD'*'LeafFrac')))))

p3.2 <- ggplot(ml.rsq.combine.sub %>% subset(source == "Meinzer et al.1999 Fig. 4"),
             aes(x = Xylem_sap_deltaD_permil, y = depth)) + #HSMTLP.80L)) +
  coord_cartesian(ylim = c(13, 0.3)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", size = 0.5) +
  geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + se,
                     xmin = Xylem_sap_deltaD_permil - se, color = sp),
                 size = 0.5, height = 0.05) +
  facet_wrap( ~ models.plot1, nrow = 1) +
  geom_text(aes(x =  Xylem_sap_deltaD_permil, y = depth, label = sp, color = sp), nudge_y = 0.15, nudge_x = 0.2,
            size = 3) +
  ylab(expression("Water Uptake Depth (m)")) + xlab(xylem.label) +
  scale_y_continuous(trans=reverselog_trans(10), breaks = ml.rsq.combine.sub$depth) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.98, npcy = 0.13, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.98, npcy = 0.05, size = 4) +
  geom_point(shape = 21, color = "white", aes(fill = sp), alpha = 1, size = 3) +
  theme(legend.position = "top",
        legend.direction = "horizontal") +
  guides(fill = "none", color = "none")
  # guides(fill = guide_legend(title = "Species"))
ggsave("psi.corr_best.depth_xylem_sap_deltaD_sp_color_Meinzer.jpeg",
       plot = p3.2, file.path(figures.folder), device = "jpeg", height = 3, width = 9, units = 'in')

p4 <- ggplot(ml.rsq.combine.sub %>% subset(source == "Meinzer et al.1999 Fig. 4" &
                                             corr.func == "gr.Psi.VPD"),
             aes(x = Xylem_sap_deltaD_permil, y = depth)) +
  geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + se,
                     xmin = Xylem_sap_deltaD_permil - se, color = sp),
                 size = 0.5, height = 0.05, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, color = "black", size = 0.5) +
  ylab(expression("Water Uptake Depth (m)")) + xlab(xylem.label) +
  scale_y_continuous(trans=reverselog_trans(10), breaks = ml.rsq.combine.sub$depth) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.9, npcy = 0.2, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.9, npcy = 0.1, size = 4) +
  geom_point(shape = 21, color = "white", aes(fill = sp), alpha = 1, size = 3.5) +
  guides(fill = guide_legend(title = "Species"))
ggsave("psi.corr_best.depth_xylem_sap_deltaD_phenology_Meinzer_gr.Psi.VPD.jpeg",
       plot = p4, file.path(figures.folder), device = "jpeg", height = 3, width = 4, units = 'in')

g3 <- ggplot(ml.rsq.combine.best %>% subset(R2 >= 0.1 & !duplicated(sp) & !is.na(depth)),
             aes(x = deci_sp.plot, y = depth)) +
  facet_wrap(. ~ corr.func) +
  geom_col(aes(fill = deciduousness)) +
  guides(fill = guide_legend(title = "")) +
  theme(legend.position = "top") +
  ylab("Depth (m)") + xlab("Species") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(trans = reverselog_trans(10), breaks = ml.rsq.combine.best$depth)
g4 <- g3 + coord_flip() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))
ggsave("psi.corr_best.depth_phenology.jpeg",
       plot = g4, file.path(figures.folder), device = "jpeg", height = 6, width = 9, units = 'in')


### Plot against hydraulic traits------

## Soil preference vs traits
hyd <- hyd %>% left_join(depth.rsq.isotopes %>% ungroup() %>%
                           subset(corr.func == "gr.Psi.VPD") %>%
                           select(sp, depth), by = "sp") %>%
  left_join(iso.1.3.join %>% subset(source == "Meinzer et al.1999 Fig. 4") %>%
              dplyr::select(sp, Xylem_sap_deltaD_permil, se), by = "sp") %>%
  droplevels()

traits.labels.table.1 <- data.frame(trait = factor(c("depth", "Xylem_sap_deltaD_permil",
                                                     "lwp.min_Predawn", "lwp.min_Diurnal", "TLP",
                                                     "KmaxS", "p50S", "p88S",
                                                     "HSMTLP", "HSM50S","HSM88S",
                                                     "HSMTLP.50S", "HSMTLP.88S",
                                                     "CWR_Total", "CWR_Xylem", "CWR_Bark",
                                                     "Felbow_Xylem", "Felbow_Bark", "HSMFelbow_Xylem", "HSMFelbow_Bark",
                                                     "Fcap_Xylem", "Fcap_Bark","WD",
                                                     "Panama.moist.pref", "Plot.swp.pref", "LMA"),
                                                   levels = c("depth", "Xylem_sap_deltaD_permil",
                                                              "lwp.min_Predawn", "lwp.min_Diurnal", "TLP",
                                                              "KmaxS", "p50S", "p88S",
                                                              "HSMTLP", "HSM50S","HSM88S",
                                                              "HSMTLP.50S", "HSMTLP.88S",
                                                              "CWR_Total", "CWR_Xylem", "CWR_Bark",
                                                              "Felbow_Xylem", "Felbow_Bark", "HSMFelbow_Xylem", "HSMFelbow_Bark",
                                                              "Fcap_Xylem", "Fcap_Bark","WD",
                                                              "Panama.moist.pref", "Plot.swp.pref", "LMA"), ordered = TRUE)) %>%
  transform(trait.plot = factor(trait, labels = c(expression(Depth[italic('Rsq')]), expression(italic(delta)^2*H[Xylem]),
                                                  expression(Psi[Predawn]), expression(Psi[min]), expression(Psi[TLP]),
                                                  expression(italic('K')['max, Stem']),  expression(Psi['50, Stem']),  expression(Psi['88, Stem']),
                                                  expression(Psi[min]*' - '*Psi[TLP]),
                                                  expression(Psi[min]*' - '*Psi['50, Stem']),
                                                  expression(Psi[min]*' - '*Psi['88, Stem']),
                                                  expression(Psi[TLP]*' - '*Psi['50, Stem']),
                                                  expression(Psi[TLP]*' - '*Psi['88, Stem']),
                                                  expression('CWR'['Total']), expression('CWR'['Xylem']), expression('CWR'['Bark']),
                                                  expression(italic('F')['Elbow, Xylem']), expression(italic('F')['Elbow, Bark']),
                                                  expression(Psi[min]*'-'*italic('F')['Elbow, Xylem']),
                                                  expression(Psi[min]*'-'*italic('F')['Elbow, Bark']),
                                                  expression(italic('F')['Cap, Xylem']), expression(italic('F')['Cap, Bark']),
                                                  expression('WD'[stem]),
                                                  expression('Panama'[wet]), expression('Plot'[wet]), "LMA")),
            trait.plot.chart = factor(trait, labels = c(expression(Depth[italic('Rsq')]), expression(italic(delta)^2*H[Xylem]),
                                                        expression(Psi[Predawn]), expression(Psi[min]), expression(Psi[TLP]),
                                                        expression(italic('K')['max, Stem']), expression(Psi['50,Stem']),  expression(Psi['88,Stem']),
                                                        expression(Psi[min]*'-'*Psi[TLP]),
                                                        expression(Psi[min]*'-'*Psi['50,Stem']),
                                                        expression(Psi[min]*'-'*Psi['88,Stem']),
                                                        expression(Psi[TLP]*'-'*Psi['50,Stem']),
                                                        expression(Psi[TLP]*'-'*Psi['88,Stem']),
                                                        expression('CWR'['Total']), expression('CWR'['Xylem']), expression('CWR'['Bark']),
                                                        expression(italic('F')['Elbow, Xylem']), expression(italic('F')['Elbow,Bark']),
                                                        expression(Psi[min]*'-'*italic('F')['Elbow,Xylem']),
                                                        expression(Psi[min]*'-'*italic('F')['Elbow,Bark']),
                                                        expression(italic('F')['Cap,Xylem']), expression(italic('F')['Cap,Bark']),
                                                        expression('WD'[stem]),
                                                        expression('Panama'[wet]), expression('Plot'[wet]), "LMA"))) %>% droplevels()

hyd.long <- hyd %>% select(-DeciLvl) %>%
  select(sp, deciduousness, deciduous, location, TLP, KmaxS, p50S, p88S,
         CWR_Total, Fcap_Xylem, CWR_Xylem, Felbow_Xylem, Fcap_Bark, CWR_Bark, Felbow_Bark, WD, LMA,
         HSM50S, HSM88S, HSMTLP, HSMFelbow_Xylem, HSMFelbow_Bark, HSMTLP.50S, HSMTLP.88S,
         Panama.moist.pref, Plot.swp.pref, lwp.min_Diurnal, lwp.min_Predawn, depth, Xylem_sap_deltaD_permil) %>% # , , se
  gather(trait, value, -sp, -deciduousness, -deciduous,  -location) %>% #, -se
  unite("deci_sp", deciduous, sp, remove = FALSE) %>%
  subset(deciduousness != "NA") %>%
  droplevels()
hyd.long.for.kruskal <- hyd.long %>% subset(trait != "Xylem_sap_deltaD_permil")
kruskal.list <- list()
for(i in unique(hyd.long.for.kruskal$trait)) {
  xx <- hyd.long.for.kruskal %>% subset(trait == i)
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

##
traits <- traits %>% left_join(depth.rsq.isotopes %>% ungroup() %>%
                                 subset(corr.func == "gr.Psi.VPD") %>%
                                 select(sp, depth), by = "sp") %>%
  left_join(iso.1.3.join %>% subset(source == "Meinzer et al.1999 Fig. 4") %>%
              dplyr::select(sp, Xylem_sap_deltaD_permil, se), by = "sp") %>%
  droplevels()

traits.labels.table.2 <- data.frame(trait = factor(c("depth", "Xylem_sap_deltaD_permil",
                                                     "KmaxL", "lwp.min_Predawn", "lwp.min_Diurnal", "TLP", "p50L", "p80L",
                                                     "HSMLWP.TLP", "HSMLWP.50L", "HSMTLP.50L",
                                                     "HSMLWP.80L", "HSMTLP.80L",
                                                     "Panama.moist.pref", "Plot.swp.pref", "SG100C_AVG", "Chl"),
                                                   levels = c("depth", "Xylem_sap_deltaD_permil",
                                                              "KmaxL", "lwp.min_Predawn", "lwp.min_Diurnal", "TLP", "p50L", "p80L",
                                                              "HSMLWP.TLP", "HSMLWP.50L", "HSMTLP.50L",
                                                              "HSMLWP.80L", "HSMTLP.80L",
                                                              "Panama.moist.pref", "Plot.swp.pref", "SG100C_AVG", "Chl"), ordered = TRUE)) %>%
  transform(trait.plot = factor(trait, labels = c(expression(Depth[italic('Rsq')]), expression(italic(delta)^2*H[Xylem]),
                                                  expression(italic(K)['max, Leaf']), expression(Psi[predawn]), expression(Psi[min]),
                                                  expression(Psi[TLP]), expression(Psi['50, Leaf']), expression(Psi['80, Leaf']),
                                                  expression(Psi[min]*' - '*Psi[TLP]),
                                                  expression(Psi[min]*' - '*Psi['50, Leaf']),
                                                  expression(Psi[TLP]*' - '*Psi['50, Leaf']),
                                                  expression(Psi[min]*' - '*Psi['80, Leaf']),
                                                  expression(Psi[TLP]*' - '*Psi['80, Leaf']),
                                                  expression('Panama'[wet]), expression('Plot'[wet]), expression('SG'[100*~degree*C]), "LMA")),
            trait.plot.chart = factor(trait, labels = c(expression(Depth[italic('Rsq')]), expression(italic(delta)^2*H[Xylem]),
                                                        expression(italic(K)['max, Leaf']), expression(Psi[predawn]), expression(Psi[min]),
                                                        expression(Psi[TLP]), expression(Psi['50,Leaf']), expression(Psi['80,Leaf']),
                                                        expression(Psi[min]*'-'*Psi[TLP]),
                                                        expression(Psi[min]*'-'*Psi['50,Leaf']),
                                                        expression(Psi[TLP]*'-'*Psi['50,Leaf']),
                                                        expression(Psi[min]*'-'*Psi['80,Leaf']),
                                                        expression(Psi[TLP]*'-'*Psi['80,Leaf']),
                                                        expression('Panama'[wet]), expression('Plot'[wet]), expression('SG'[100*~degree*C]), "LMA")))




## Kunert traits
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
select.traits <- c("lwp.min_Predawn", "lwp.min_Diurnal", "TLP", "KmaxS", "p50S", "p88S",
                   "HSMTLP", "HSM50S", "HSM88S", "HSMTLP.50S", "HSMTLP.88S")

traits.long <- traits.long %>%
  unite("deci_sp", deciduous, sp, remove = FALSE) %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(deciduousness)]), ordered=TRUE),
         deci_sp.plot = factor(deci_sp, levels=unique(deci_sp[order(deciduousness)]), ordered=TRUE))
# just for sp with hyd.traits, but traits.long does not have all those sp, and hab preference and WSG traits will be missed
## so beginning with those other traits
traits.wide <- traits.long %>% select(-trait.plot, -trait.plot.chart) %>%
  pivot_wider(names_from = trait, values_from = value)
traits.long.hyd <- data.frame(sp = unique(c(hyd$sp, traits$sp))) %>%
  full_join(bci.traits %>% select(sp, form1, SG100C_AVG), by = "sp") %>%
  full_join(deci %>% select(-sp4), by = "sp") %>%
  subset(sp %in% unique(c(depth.rsq.isotopes$sp, hyd$sp, traits$sp))) %>%
  left_join(traits.wide %>% select(-form1, -deciduous, -deciduousness,
                                   -SG100C_AVG, -Panama.moist.pref, -Plot.swp.pref), by = "sp") %>%
  pivot_longer(cols = c(-sp, -form1, -deciduous, -deciduousness, -DeciLvl,
                        -deciduousness.label, -sp.plot, -deci_sp, -deci_sp.plot),
               names_to = "trait", values_to = "value") %>%
  unite("deci_sp", deciduous, sp, remove = FALSE) %>%
  left_join(traits.labels.table.2, by = "trait") %>%
  mutate(sp.plot = factor(sp, levels = unique(sp[order(deciduousness)]), ordered=TRUE),
         deci_sp.plot = factor(deci_sp, levels=unique(deci_sp[order(deciduousness)]), ordered=TRUE)) %>%
  droplevels()

save(traits.long, file = file.path(results.folder, "kunert.traits.key.long_depth_isotopes.RData"))
save(traits.long.hyd, file = file.path(results.folder, "kunert.traits.key.long_in_Wolfe_traits_species_list_depth_isotopes.RData"))


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
select.traits.1 <- c("TLP", "KmaxS","p50S", "p88S", "HSMTLP", "HSM50S",
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
                     "HSMLWP.TLP", "Panama.moist.pref", "Plot.swp.pref", "SG100C_AVG", "Chl")

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

depth.traits.hyd <- hyd.long %>%
  subset(trait %in% c(select.traits.1, select.traits.2)) %>%
  subset(trait != "depth") %>% left_join(hyd.pairs.2 %>% select(sp, `Depth[italic("Rsq")]`), by = "sp")

depth.traits.hyd.plot <- ggplot(depth.traits.hyd,
                                aes(y = `Depth[italic("Rsq")]`, x = value)) +
  geom_smooth(method = "lm") +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 0.8, size = 2.5) +
  scale_y_reverse() +
  coord_cartesian(ylim = c(10, 0)) +
  ylab("Effective Rooting Depth (m)") + xlab("") +
  facet_wrap(. ~ trait.plot, scales = "free_x", labeller = label_parsed) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.10, npcy = 0.25, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.10, npcy = 0.1, size = 4) +
  theme(panel.spacing = unit(1, "lines"))
ggsave(file.path(figures.folder, paste0("Wolfe_traits_depth.jpeg")),
       plot = depth.traits.hyd.plot, height = 6, width = 6, units ='in')

depth.traits.kunert <- traits.long %>%
  subset(trait %in% c(select.traits.3, select.traits.4)) %>%
  subset(trait != "depth") %>% left_join(traits.pairs.2 %>% select(sp, `Depth[italic("Rsq")]`), by = "sp")

depth.traits.kunert.plot <- ggplot(depth.traits.kunert %>%
                                     subset(trait != "Xylem_sap_deltaD_permil"),
                                   aes(y = `Depth[italic("Rsq")]`, x = value)) +
  geom_smooth(method = "lm") +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 0.8, size = 2.5) +
  scale_y_reverse() +
  coord_cartesian(ylim = c(10, 0)) +
  ylab("Effective Rooting Depth (m)") + xlab("") +
  facet_wrap(. ~ trait.plot, scales = "free_x", labeller = label_parsed) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.10, npcy = 0.25, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.10, npcy = 0.1, size = 4) +
  theme(panel.spacing = unit(1.5, "lines"))
ggsave(file.path(figures.folder, paste0("Kunert_traits_depth.jpeg")),
       plot = depth.traits.kunert.plot, height = 6.5, width = 7, units ='in')

traits.labels.select <- data.frame(trait = factor(c("KmaxS", "TLP", "p88S", "HSM88S"),
                                        levels = c("KmaxS", "TLP", "p88S", "HSM88S"), ordered = TRUE)) %>%
  transform(trait.plot = factor(trait, labels = c(expression(italic('K')['max, Stem']), expression(Psi[TLP]),
                                                  expression(Psi['88, Stem']),
                                                  expression(Psi[min]*' - '*Psi['88, Stem']))))

hyd.error <- hyd %>% select(sp, KmaxS_se, vc_b_se, vc_a_se, tlp_sd) %>%
  rename(KmaxS = KmaxS_se, vc_b = vc_b_se, vc_a = vc_a_se, TLP = tlp_sd) %>%
  gather(trait, se, -sp) ## But note that for TLP it's not se but sd

erd.stem.traits <- depth.traits.hyd %>% subset(trait %in% c("KmaxS", "TLP", "p88S", "HSM88S")) %>%
  select(deci_sp, sp, trait, `Depth[italic("Rsq")]`, value) %>%
  # bind_rows(depth.traits.kunert %>% subset(trait == "KmaxL") %>%
  #             select(deci_sp, sp, trait, `Depth[italic("Rsq")]`, value)) %>%
  left_join(traits.labels.select %>% select(trait, trait.plot), by = "trait") %>%
  droplevels() %>%
  left_join(hyd.error, by = c("sp", "trait"))

tlp.hyd.kunert <- hyd.pairs.1 %>% select(sp, `Psi[TLP]`) %>%
  rename(`Wolfe et al. Psi[TLP]` = `Psi[TLP]`) %>%
  left_join(traits.pairs.2 %>% select(sp, `Psi[TLP]`) %>%
              rename(`Kunert et al. Psi[TLP]` = `Psi[TLP]`), by = "sp")
tlp.min.plot = min(c(tlp.hyd.kunert$`Kunert et al. Psi[TLP]`, tlp.hyd.kunert$`Wolfe et al. Psi[TLP]`), na.rm = TRUE)
tlp.max.plot = max(c(tlp.hyd.kunert$`Kunert et al. Psi[TLP]`, tlp.hyd.kunert$`Wolfe et al. Psi[TLP]`), na.rm = TRUE)
tlp.hyd.kunert.plot <- ggplot(tlp.hyd.kunert, aes(x = `Wolfe et al. Psi[TLP]`, y = `Kunert et al. Psi[TLP]`)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.05, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 5) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.05, npcy = 0.82, size = 5) +
  ylim(c(tlp.min.plot, tlp.max.plot)) + xlim(c(tlp.min.plot, tlp.max.plot))+
  geom_text(data = data.frame(x = -1.2, y = -1), aes(x = x, y = y, label = "1:1"))
ggsave(file.path(figures.folder, paste0("TLP_Kunert_by_Wolfe.jpeg")),
       plot = tlp.hyd.kunert.plot, height = 3, width = 3, units ='in')

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
  scale_x_discrete(name = "",
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
n.threshold = 50
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
m2.1 <- ggplot(mrate.long %>%
                 subset(size == "large" & avg.abund >= n.threshold & !is.na(deciduousness)),
               aes(x = censusint.m, y = diff.mrate, color = avg.abund)) +
  scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
                       low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
  facet_grid(size.num ~ ., scales = "free_y") +
  geom_hline(aes(yintercept = 0), color = "blue", size = 0.5) +
  geom_boxplot(stat = "boxplot", notch = TRUE) +
  scale_fill_brewer(name = "", palette = "Greens", direction = -1) +
  theme(legend.position = "top") +
  ylab(y.label.1) + xlab("Census Interval") +
  ggtitle("Mortality Anomaly by Census Interval") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(width = 0.05, shape = 21, fill = "darkgray", color = "black", show.legend = FALSE, alpha = 0.7)
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_by_size_large_points.jpeg")), plot = m2.1, height = 3, width = 6, units='in')

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
  ggtitle(expression('Species leafless in early wet season, with increasing HSM '*Psi[TLP]*' - '*Psi['88,Stem']))
  #  does not work: fill = "Facultative Deciduous"
  # guides(fill = guide_legend(title = "Deciduousness",
  #                            override.aes = list(fill = c("Facultative Deciduous" = "#35B779FF"))))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_deci_HSMTLP.88S_increasing_spp_leafless in early wet season.jpeg")),
       plot = m3.1, height = 3, width = 9, units='in')

m4.1 <- m3.base %+% subset(mrate.long.hyd, sp %in% hyd$sp & size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Mortality for Evergreen Species with increasing HSM '*Psi[TLP]*' - '*Psi['88,Stem']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_HSMTLP.88S_increasing_evergreens.jpeg")),
       plot = m4.1, height = 3, width = 9, units='in')

mrate.long.hyd <- mrate.long.hyd %>% mutate(sp.plot = factor(sp, levels=unique(sp[order(-p88S)]), ordered=TRUE))
m3.2 <- m3.base + geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Species leafless in early wet season, with increasingly more negative '*Psi['88,Stem']))
ggsave(file.path(paste0(figures.folder, "/sp_Mortality_rate_by_period_deci_p88S_increasing_spp_leafless in early wet season.jpeg")),
       plot = m3.2, height = 3, width = 9, units='in')

m4.2 <- m3.base %+% subset(mrate.long.hyd, sp %in% hyd$sp & size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Mortality for Evergreen Species with increasingly more negative '*Psi['88,Stem']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_p88S_increasing_evergreens.jpeg")),
       plot = m4.2, height = 3, width = 9, units='in')
mrate.long.hyd <- mrate.long.hyd %>% mutate(sp.plot = factor(sp, levels=unique(sp[order(-HSM50S)]), ordered=TRUE))
m4.3 <- m3.base %+% subset(mrate.long.hyd, sp %in% hyd$sp & size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Mortality for Evergreen Species with increasing HSM '*Psi[min]*' - '*Psi['50,Stem']))
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
  ggtitle(expression('Mortality for Evergreen Species with increasing HSM '*Psi[min]*' - '*Psi['50, Leaf']))
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

## Ordered along Rooting Depth Index
mrate.long.depth <- mrate.long %>%
  left_join(subset(depth.rsq.isotopes, corr.func == "gr.Psi.Rad.VPD"), by = "sp") %>%
  left_join(bci.traits %>% select(form1, sp), by = "sp")
  mutate(sp.plot = factor(sp, levels = unique(sp[order(depth)]), ordered = TRUE))
m4.6 <- m3.base %+% subset(mrate.long.depth, size == "large" &
                             form1 == "T" &!is.na(deciduousness) &!is.na(depth)) +
    geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
    # facet_grid(deciduousness ~ censusint.m) +
    theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 0.5, hjust = 1, size = 4)) +
    ggtitle(expression('Mortality for Canopy Species with increasing Rooting Depth Index'))
ggsave(file.path(paste0(figures.folder,
                          "/sp_Mortality_rate_by_period_best_corr_depth_increasing_canopy_sp.jpeg")),
         plot = m4.6, height = 3, width = 10, units='in')

## For each species get % days spent below the kl80 cutoff in the depth of best-correlation (and the depths below)?
## Test whether that explains mortality in that census (in reality this would be cumulative)


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

g2.1 <- ggplot(grate.long %>%
                 subset(size == "large" & !is.na(deciduousness)),
               aes(x = censusint.m, y = median)) +
  scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
                       low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
  facet_grid(size.num ~ ., scales = "free_y") +
  geom_hline(aes(yintercept = 0), color = "blue", size = 0.5) +
  geom_boxplot(stat = "boxplot", notch = TRUE) +
  scale_fill_brewer(name = "", palette = "Greens", direction = -1) +
  theme(legend.position = "top") +
  ylab(y.label.2) + xlab("Census Interval") +
  ggtitle("Growth Anomaly by Census Interval") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(width = 0.05, shape = 21, fill = "darkgray", color = "black", show.legend = FALSE, alpha = 0.7)
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_by_size_large_points.jpeg")), plot = g2.1, height = 3, width = 6, units='in')

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
  ggtitle(expression('Species leafless in early wet season, with increasing HSM '*Psi[TLP]*' - '*Psi['88,Stem']))
#  does not work: fill = "Facultative Deciduous"
# guides(fill = guide_legend(title = "Deciduousness",
#                            override.aes = list(fill = c("Facultative Deciduous" = "#35B779FF"))))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_deci_HSMTLP.88S_increasing_spp_leafless in early wet season.jpeg")),
       plot = g3.1, height = 4, width = 9, units='in')

g4.1 <- g3.base %+% subset(grate.long.hyd, sp %in% hyd$sp & size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = median, fill = deciduousness)) +
  ggtitle(expression('Growth for Evergreen Species with increasing HSM '*Psi[TLP]*' - '*Psi['88,Stem']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_HSMTLP.88S_increasing_evergreens.jpeg")),
       plot = g4.1, height = 4, width = 9, units='in')

grate.long.hyd <- grate.long.hyd %>% mutate(sp.plot = factor(sp, levels=unique(sp[order(-p88S)]), ordered=TRUE))
g3.2 <- g3.base + geom_col(aes(x = sp.plot, y = median, fill = deciduousness)) +
  ggtitle(expression('Species leafless in early wet season, with increasingly more negative '*Psi['88,Stem']))
ggsave(file.path(paste0(figures.folder, "/sp_Growth_rate_by_period_deci_p88S_increasing_spp_leafless in early wet season.jpeg")),
       plot = g3.2, height = 4, width = 9, units='in')

g4.2 <- g3.base %+% subset(grate.long.hyd, sp %in% hyd$sp & size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = median, fill = deciduousness)) +
  ggtitle(expression('Growth for Evergreen Species with increasingly more negative '*Psi['88,Stem']))
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


## Plot mortality by time spent below a threshold in the preferred depth-------

names.mfac <- names(get.mfac.ls)
mfac.interval <- vector(mode = "list", length = length(names.mfac))
names(mfac.interval) <- names.mfac  # "psi.p50.g1", "psi.p50.g2"

for (i in 1:length(names.mfac)) { #
  mfac.interval[[i]] <- lapply(lapply(AB.sp.ls, get.mfac.ls[[i]]),
                               as.data.frame) %>%
    bind_rows(.id = "sp")
}

mrate.psi <- mfac.interval.long <- vector(mode = "list", length = length(names.mfac))
names(mfac.interval.long) <- names(mrate.psi) <- names.mfac
for (i in 1:length(names.mfac)) {
  mfac.interval.long[[i]] <- mfac.interval[[i]] %>%
    pivot_longer(cols = c(-sp, -interval),
                 names_to = "depth", values_to = "mfac") %>%
    rename(interval.num =  interval) %>%
    mutate(depth = as.numeric(depth),
           interval.num =  as.numeric(as.character(interval.num)),
           size = "large")
  mrate.psi[[i]] <- mrate.long %>%
    left_join(mfac.interval.long[[i]], by = c("interval.num", "sp", "size"))
}

## Ordered along Rooting Depth Index
mfac.on <- "mr.kl50.I"
mrate.depth <-
  adult.mrate.long %>% mutate(size = "large") %>%
  # mrate.long %>%
  left_join(subset(depth.rsq.isotopes, corr.func == "gr.Psi.VPD") %>%
              rename(rdi.gr = depth) %>%
              dplyr::select(sp, size, rdi.gr), by = c("sp", "size")) %>%
  left_join(subset(depth.rsq.isotopes, corr.func == "mr.Psi.VPD.I") %>%
              rename(rdi.mr = depth) %>%
              dplyr::select(sp, size, rdi.mr), by = c("sp", "size")) %>%
  left_join(bci.traits %>% dplyr::select(form1, sp), by = "sp") %>%
  # mutate(sp.plot = factor(sp, levels = unique(sp[order(rdi.gr)]), ordered = TRUE)) %>%
  mutate(size = as.character(size)) %>%
  subset(size == "large" & form1 == "T") %>% droplevels()
mrate.mfac.depth <- mrate.depth %>%
  right_join(mfac.interval.long[[mfac.on]] %>%
               mutate(censusint.m = recode(interval.num, `1` = "1982-85", `2` = "1985-90",
                                           `3` = "1990-95", `4` = "1995-00", `5` = "2000-05", `6` = "2005-10", `7` = "2010-15")),
             by = c("censusint.m", "sp", "size")) %>%
  mutate(sp_size = paste(sp, size, sep = "_")) %>%
  group_by(sp, size, censusint.m) %>%
  mutate(mfac.soil.column = sum(mfac, na.rm = TRUE)) %>%
  ungroup(sp, censusint.m)
mrate.mfac.depth.to.rdi.gr <- mrate.mfac.depth %>%
  group_by(sp, size) %>%
  subset(!depth > rdi.gr) %>%
  ungroup(sp, size) %>%
  group_by(sp, size, censusint.m) %>%
  mutate(mfac.soil.column.gr = sum(mfac, na.rm = TRUE)) %>%
  ungroup(sp, size, censusint.m)
mrate.mfac.depth.to.rdi.gr.total.int <- mrate.mfac.depth.to.rdi.gr %>%
  dplyr::select(sp, size, mfac.soil.column.gr, censusint.m, mrate, diff.mrate, depth) %>%
  group_by(sp, size, censusint.m) %>%
  summarise(mfac.soil.column.gr = mean(mfac.soil.column.gr, na.rm = TRUE),
            mrate = mean(mrate, na.rm = TRUE),
            diff.mrate = mean(diff.mrate, na.rm = TRUE),
            depth = mean(depth, na.rm = TRUE)) %>%
  ungroup(sp, size, censusint.m)
mrate.mfac.depth.to.rdi.gr.total <- mrate.mfac.depth.to.rdi.gr.total.int %>%
  group_by(sp, size) %>%
  summarise(mfac.soil.column.total.gr = sum(mfac.soil.column.gr, na.rm = TRUE)) %>%
  ungroup(sp, size)
mrate.mfac.depth.to.rdi.gr.study <- mrate.mfac.depth.to.rdi.gr %>%
  subset(depth == rdi.gr) %>%
  group_by(sp, size) %>%
  summarise(mfac.total.gr = sum(mfac, na.rm = TRUE),
         mrate.sum = sum(mrate, na.rm = TRUE)) %>%
  ungroup(sp, size)

save(mrate.depth, file = file.path(results.folder, "mrate.depth.Rdata"))
save(mrate.mfac.depth, file = file.path(results.folder, "mrate.mfac.depth.Rdata"))
## for rdi.mr
mrate.mfac.depth.to.rdi.mr <- mrate.mfac.depth %>%
  group_by(sp, size) %>%
  subset(!depth > rdi.mr) %>%
  ungroup(sp, size) %>%
  group_by(sp, size, censusint.m) %>%
  mutate(mfac.soil.column.mr = sum(mfac, na.rm = TRUE)) %>%
  ungroup(sp, size, censusint.m)
mrate.mfac.depth.to.rdi.mr.total.int <- mrate.mfac.depth.to.rdi.mr %>%
  dplyr::select(sp, size, mfac.soil.column.mr, censusint.m, mrate, diff.mrate, depth) %>%
  group_by(sp, size, censusint.m) %>%
  summarise(mfac.soil.column.mr = mean(mfac.soil.column.mr, na.rm = TRUE),
            mrate = mean(mrate, na.rm = TRUE),
            diff.mrate = mean(diff.mrate, na.rm = TRUE),
            depth = mean(depth, na.rm = TRUE)) %>%
  ungroup(sp, size, censusint.m)
mrate.mfac.depth.to.rdi.mr.total <- mrate.mfac.depth.to.rdi.mr.total.int %>%
  group_by(sp, size) %>%
  summarise(mfac.soil.column.total.mr = sum(mfac.soil.column.mr, na.rm = TRUE)) %>%
  ungroup(sp, size)
mrate.mfac.depth.to.rdi.mr.study <- mrate.mfac.depth.to.rdi.mr %>%
  subset(depth == rdi.mr) %>%
  group_by(sp, size) %>%
  summarise(mfac.total.mr = sum(mfac, na.rm = TRUE),
            mrate.sum = sum(mrate, na.rm = TRUE)) %>%
  ungroup(sp, size)

mrate.mfac.column <- mrate.mfac.depth %>%
  group_by(sp, size, censusint.m) %>%
  mutate(mfac.soil.column = sum(mfac, na.rm = TRUE)) %>%
  ungroup(sp, size, censusint.m)
mrate.mfac.column.total.int <- mrate.mfac.column %>%
  dplyr::select(sp, size, mfac.soil.column, censusint.m, mrate, diff.mrate, depth) %>%
  group_by(sp, size, censusint.m) %>%
  summarise(mfac.soil.column = mean(mfac.soil.column, na.rm = TRUE),
            mrate = mean(mrate, na.rm = TRUE),
            diff.mrate = mean(diff.mrate, na.rm = TRUE),
            depth = mean(depth, na.rm = TRUE)) %>%
  ungroup(sp, size, censusint.m)

mfac.plot.7 <- ggplot(mrate.mfac.depth.to.rdi.gr.study,
       aes(x = mfac.total.gr, y = mrate.sum)) +
  # facet_wrap(censusint.m ~ ., nrow = 1) +
  geom_point() +
  geom_smooth(method = "lm") +
  ylab(expression('Total Mortality Rate (% '*'year'^1*')')) +
  xlab(expression('Days '*Psi['Soil,z = ERD']*'<'*Psi['P80,Leaf'])) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4)
ggsave(file.path(paste0(figures.folder,
                        "/sp_pc_during_the_study_period_total days below kl80_in_z through_rdi.mr.jpeg")),
       plot = mfac.plot.7, height = 3, width = 3, units='in')

## Mean mortality rate vs. RDI
# mean.mrate.rdi <- mrate.mfac.depth %>%
#   subset(depth == rdi) %>%
#   group_by(sp, size, rdi) %>%
#   summarise(mrate.mean = mean(mrate, na.rm = TRUE),
#             mrate.se = sd(mrate, na.rm = TRUE)/sqrt(n()))
y.label.1 <- expression(atop(Mortality[Interval], '-'~Mean[Mortality]~('%'*yr^{-1})))
mfac.plot.8 <- ggplot(mrate.mfac.depth %>% subset(depth == rdi.mr),
                      aes(y = diff.mrate, x = mfac)) +
  geom_point() + ylab(y.label.1) +
  xlab(expression('Days '*Psi['Soil, z = RDI.mr']*' < '*Psi['crit'])) +
  geom_smooth(method = "lm") +
  facet_grid(. ~ censusint.m ) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/diff.mortality_rate_mfac in rdi.mr.jpeg")),
       plot = mfac.plot.8, height = 3, width = 9, units='in')

mfac.plot.9 <- mfac.plot.8 %+% subset(mrate.mfac.depth, depth == rdi.gr) +
  xlab(expression('Days '*Psi['Soil, z = ERD']*' < '*Psi['P80, Leaf'])) +
  geom_point(aes(color = depth)) +
  scale_color_continuous(trans = "reverse", guide = guide_colorbar(title = "Depth\n(m)"))
ggsave(file.path(paste0(figures.folder, "/diff.mortality_rate_mfac in rdi.gr.jpeg")),
       plot = mfac.plot.9, height = 3, width = 9, units='in')

mfac.plot.10 <- ggplot(mrate.mfac.depth.to.rdi.mr.total.int,
                      aes(y = diff.mrate, x = mfac.soil.column.mr)) +
  geom_point(aes(color = depth)) +
  scale_color_continuous(trans = "reverse", guide = guide_colorbar(title = "Depth\n(m)")) +
  ylab(y.label.1) +
  xlab(expression('Days '*Psi['Soil, z <= RDI.mr']*' < '*Psi['P80, Leaf'])) +
  geom_smooth(method = "lm") +
  facet_grid(. ~ censusint.m ) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/diff.mortality_rate_mfac upto rdi.mr.jpeg")),
       plot = mfac.plot.10, height = 3, width = 9, units='in')

mfac.plot.11 <- ggplot(mrate.mfac.depth.to.rdi.gr.total.int,
                       aes(y = diff.mrate, x = mfac.soil.column.gr)) +
  geom_point(aes(color = depth)) +
  scale_color_continuous(trans = "reverse", guide = guide_colorbar(title = "Depth\n(m)")) +
  ylab(y.label.1) +
  xlab(expression('Mean Days '*Psi['Soil, z <= ERD']*' < '*Psi['P80, Leaf'])) +
  geom_smooth(method = "lm") +
  facet_grid(. ~ censusint.m ) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/diff.mortality_rate_mfac upto rdi.gr.jpeg")),
       plot = mfac.plot.11, height = 3, width = 9, units='in')

mfac.plot.12 <- ggplot(mrate.mfac.column.total.int,
                       aes(y = mrate, x = mfac.soil.column)) +
  geom_point() + ylab(expression('Mortality Rate (% '*'year'^1*')')) +
  xlab(expression('Mean Days '*Psi['Soil, all z']*' < '*Psi['P80, Leaf'])) +
  geom_smooth(method = "lm") +
  facet_grid(. ~ censusint.m ) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_mfac all z.jpeg")),
       plot = mfac.plot.12, height = 3, width = 9, units='in')

mfac.plot.14 <- ggplot(mrate.mfac.depth %>% subset(depth == rdi.mr),
                      aes(y = mrate, x = rdi.mr)) +
  geom_point() + ylab(expression('Mortality Rate (% '*'year'^1*')')) +
  xlab("Depth best-correlated with Mortality Rate (m)") +
  geom_smooth(method = "lm") +
  facet_grid(. ~ censusint.m ) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by rdi.mr.jpeg")),
       plot = mfac.plot.14, height = 3, width = 9, units='in')


# mfac.plot.16 <- mfac.plot.15 %+% subset(mrate.mfac.depth, depth == rdi.gr) +
#   facet_grid(. ~ censusint.m) +
#   ylab(expression('Mortality Rate (% '*'year'^1*')'))
# ggsave(file.path(paste0(figures.folder, "/mortality_rate_by rdi.gr_interval.jpeg")),
#        plot = mfac.plot.16, height = 3, width = 9, units='in')

adult.mrate.mean <- adult.mrate.long %>%
  group_by(sp, deciduousness) %>%
  summarize_at(vars(mean.mrate, avg.abund, mean.grate), mean, na.rm = TRUE) %>%
  mutate(mean.mrate = ifelse(!is.finite(mean.mrate),
                        rep(NA, length(mean.mrate)), mean.mrate)) %>%
  mutate(size = "large") %>%
  left_join(subset(depth.rsq.isotopes, corr.func == "gr.Psi.VPD") %>%
              dplyr::select(sp, size, depth), by = c("sp", "size")) %>%
  left_join(bci.traits %>% dplyr::select(form1, sp), by = "sp") %>%
  subset(size == "large" & form1 == "T")

pm.2 <- ggplot(adult.mrate.mean %>%
                 subset(avg.abund >= n.threshold & deciduousness %in%
                          c("Facultative Deciduous", "Evergreen")),
               aes(x = depth, y = mean.mrate)) +
  geom_point() +
  scale_x_continuous(trans="sqrt", breaks = soil.depths[-c(2,3, 4, 6, 7, 9)]) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.9, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 5) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                  npcx = 0.9, npcy = 0.8, size = 5) +
  ylab(expression('Mean Mortality Rate (% '*'year'^1*')')) +
  xlab("Effective Rooting Depth (m)") +
  facet_grid(. ~ deciduousness)  +
  geom_smooth(method = "lm", se = TRUE) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  ggtitle("Adult tree mortality (>= 10 cm DBH)") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(paste0(figures.folder, "/adult_Mortality_vs_udi_with_outliers_avg.abund_above",
                        n.threshold, "_sp_deci.jpeg")), plot = pm.2, height = 4, width = 6, units='in')

pg.2 <- ggplot(adult.mrate.mean,
               aes(x = depth, y = mean.grate)) +
  geom_point() +
  scale_x_continuous(trans="sqrt", breaks = soil.depths[-c(2,3, 4, 6, 7, 9)]) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.9, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 5) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                  npcx = 0.9, npcy = 0.8, size = 5) +
  ylab(expression('Mean Growth Rate (cm.'*'year'^1*')')) +
  xlab("Effective Rooting Depth (m)") +
  geom_smooth(method = "lm", se = TRUE) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  ggtitle("Adult Trees (>=10cm DBH)") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(paste0(figures.folder, "/adult_Growth_vs_udi_with_outliers_sp.jpeg")), plot = pg.2, height = 4, width = 4, units='in')

pg.2.deci <- pg.2 %+% subset(adult.mrate.mean,!is.na(deciduousness)) +
  facet_grid(. ~ deciduousness)
ggsave(file.path(paste0(figures.folder, "/adult_Growth_vs_udi_with_outliers_sp_deci.jpeg")), plot = pg.2.deci, height = 4, width = 8, units='in')

#******************************************************
### Yearly psi dynamics versus climatology-------
#******************************************************

## psi does not fit a normal/exp/gamma distribution # so treating non-parametrically
psi.stat.1 <- psi %>%
  group_by(interval.yrs, date, depth) %>%
  summarise(median = -median(psi, na.rm = TRUE),
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
  geom_line(aes(x = date, y = median, group = depth, color = depth), size = 0.3) +
  scale_color_continuous(trans = "reverse", guide = guide_colorbar(title = "Depth\n(cm)")) +
  geom_vline(xintercept = census.beg) + # as.Date(paste0(c(1990:2015), "-01-01")) +
  scale_y_reverse(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 10, 15), limits = c(2.5, 0)) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%Y")) +
  coord_cartesian(xlim = c(as.Date("1990-01-01"), as.Date("2018-12-31"))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  # facet_grid(interval.yrs ~ .) +
  ylab("-Soil Water Potential [MPa]") +
  theme(text = element_text(size = 12)) + xlab("")
  # ggtitle("PSI:Violins of interval medians of best-fit ensembles")
ggsave("psi_model_daily_bestfit_params.top.few_CI_full.jpeg", plot = plot.psi.stat.1,
       file.path(figures.folder), device = "jpeg", height = 3, width = 20, units='in')

psi.2 <- psi %>%
  mutate(interval.yrs.to.plot = forcats::fct_explicit_na(cut(date, include.lowest = TRUE, breaks = cut.breaks.2,
                                    labels = cut.labels.2, right = TRUE)))

psi.stat.2 <- psi.2 %>%
  group_by(interval.yrs.to.plot, date, depth) %>%
  summarise(median = median(psi, na.rm = TRUE),
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
  geom_line(aes(x = days, y = median, group = depth, color = depth), size = 0.3) +
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
  summarise(median = median(psi, na.rm = TRUE),
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
  summarise(median = median(psi, na.rm = TRUE),
            q97.5 = quantile(psi, probs = 0.975),
            q2.5 = quantile(psi, probs = 0.025),
            q5 = quantile(psi, probs = 0.05)) %>%
  ungroup(doy, year, depth) %>%
  mutate(doy = as.numeric(doy))
psi.stat.5 <- psi.stat.4 %>%
  group_by(doy, depth) %>%
  summarise(median.clim = median(median, na.rm = TRUE),
            q97.5.clim = quantile(median, probs = 0.975),
            q2.5.clim = quantile(median, probs = 0.025),
            q10.clim = quantile(median, probs = 0.1),
            q5.clim = quantile(median, probs = 0.05)) %>%
  ungroup(doy, depth) %>%
  mutate(doy = as.numeric(doy))
rectangles.3 <- data.frame(
  xmin = 120,
  xmax = 335,
  ymin = 0,
  ymax = -2.5
)
## For selected depths used in inverse modeling
# depth.sub, depth.breaks, depth.labels are defined above in section:
### Calculate Correlation of growth rates with psi by depth

psi.stat.4.select <- psi %>%
  mutate(doy = format(date, "%j"),
         year = format(date, "%Y")) %>%
  ## due to time nneded for model initialisation soil layers deeper than 2.9 m
  ## are not fully recharged until the end of 1992, layer 2.9 m recharges by the end of 1992
  subset(!year %in% c("1990")) %>%
  subset(!year %in% c("1990", "1991") | depth == 2.9) %>%
  subset(!year %in% c("1990", "1991", "1992") | depth < 2.9) %>%
  subset(depth %in% depth.sub) %>%
  mutate(depth = cut(depth, include.lowest = TRUE, breaks = depth.breaks,
      labels = depth.labels, right = TRUE)) %>%
  group_by(date, doy, year, interval.yrs, depth) %>%
  summarise(median = median(psi, na.rm = TRUE),
            q97.5 = quantile(psi, probs = 0.975),
            q2.5 = quantile(psi, probs = 0.025),
            q5 = quantile(psi, probs = 0.05)) %>%
  ungroup(doy, year, depth) %>%
  mutate(doy = as.numeric(doy))

psi.stat.5.select <- psi.stat.4.select %>%
  group_by(doy, depth) %>%
  summarise(median.clim = median(median, na.rm = TRUE),
            q97.5.clim = quantile(median, probs = 0.975),
            q2.5.clim = quantile(median, probs = 0.025),
            q10.clim = quantile(median, probs = 0.1),
            q5.clim = quantile(median, probs = 0.05)) %>%
  ungroup(doy, depth) %>%
  mutate(doy = as.numeric(doy))

source("code/Utilities/plot.ticks.R")
plot.psi.stat.5.base <- ggplot(psi.stat.5 %>% droplevels()) +
  # geom_rect(data=rectangles.3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
  #           fill='gray80', alpha=0.8) +
  scale_x_continuous(breaks = c(seq(0, 360, by = 60))) +
  coord_cartesian(ylim = c(-2.5, 0)) +
  theme(panel.grid.major.y = element_line()) +
  ylab(expression(Psi[soil]*~"(MPa)")) + xlab("Day of the Year")
plot.psi.stat.5 <- plot.psi.stat.5.base +
  geom_ribbon(aes(x = doy, ymin = q2.5.clim, ymax = q97.5.clim), alpha = 0.3, fill = "grey20") +
  geom_line(aes(x = doy, y = median.clim, group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
  guides(color = guide_legend(title = "Depth(m)", order = 2, override.aes = list(size = 3)))
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_censuspanels_climatology.jpeg",
       plot = plot.ticks(plot.psi.stat.5),
       file.path(figures.folder), device = "jpeg", height = 3, width = 5, units='in')

plot.psi.stat.5.over <- plot.psi.stat.5.base %+%
  subset(psi.stat.5, depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)) +
  geom_line(aes(x = doy, y = median.clim, linetype = "climatology", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
  geom_line(data = subset(psi.stat.4, year == "2016" & depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)),
            aes(x = doy, y = median, linetype = "2016", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
  guides(color = guide_legend(title = "Depth(m)", legend.position = "right", order = 1, override.aes = list(size = 3)),
         linetype = guide_legend(order = 2, title = NULL, legend.position = "top", override.aes =
                                   list(linetype = c("climatology" = "dashed", "2016" = "solid")))) +
  coord_cartesian(ylim = c(-3, 0)) + ggtitle("2016")
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_censuspanels_climatology_over.jpeg",
       plot = plot.ticks(plot.psi.stat.5.over),
       file.path(figures.folder), device = "jpeg", height = 3, width = 5, units='in')

pdf(paste0(figures.folder, "/psi_model_daily_bestfit_params.top.few_CI_full_censuspanels_climatology_over_by_year.pdf"), height = 4, width = 7)
for (i in unique(psi.stat.4$year)) {
  plot.psi.stat.5.yr <- plot.psi.stat.5.base %+% subset(psi.stat.5, depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)) +
    geom_line(aes(x = doy, y = median.clim, linetype = "climatology", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
    geom_line(data = subset(psi.stat.4, year == i & depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)),
              aes(x = doy, y = median, linetype = "Year", group = as.factor(depth), color = as.factor(depth)), size = 0.5) +
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
  mutate(below.q10 = ifelse(median < q10.clim, median, NA),
         below.q5 = ifelse(median < q5.clim, median, NA),
         below.q2.5 = ifelse(median < q2.5.clim, median, NA),
         depth_year = paste(depth, year, sep = "_"))

psi.stat.4.select <- psi.stat.4.select %>%
  mutate(interval.yrs.2 = forcats::fct_explicit_na(cut(date, include.lowest = TRUE, breaks = c(cut.breaks, max(psi.stat.4$date, na.rm = TRUE)),
                                                       labels = c(cut.labels.2, "2015-2018"), right = TRUE))) %>%
  left_join(psi.stat.5.select, by = c("doy", "depth")) %>%
  mutate(below.q10 = ifelse(median < q10.clim, median, NA),
         below.q5 = ifelse(median < q5.clim, median, NA),
         below.q2.5 = ifelse(median < q2.5.clim, median, NA),
         depth_year = paste(depth, year, sep = "_"))

save(psi.stat.4, file = file.path(results.folder, "psi.stat.4.select.Rdata"))
save(psi.stat.4.select, file = file.path(results.folder, "psi.stat.4.select.Rdata"))

pct.drought.days <- psi.stat.4 %>%
  mutate(season = ifelse(doy < 120, "Dry Season", "Wet Season")) %>%
  group_by(depth, interval.yrs, interval.yrs.2, season) %>%
  summarise(pct.days.below.q10 = 100*round(sum(!is.na(below.q10))/n(), 3),
            pct.days.below.q5 = 100*round(sum(!is.na(below.q5))/n(), 3),
            pct.days.below.q2.5 = 100*round(sum(!is.na(below.q2.5))/n(), 3),
            pct.days.above.q5 = 100 - pct.days.below.q5,
            pct.days.above.q2.5 = 100 - pct.days.below.q2.5,
            pct.days.below.q5.0.5 = 100*round(sum(!is.na(below.q5) & below.q5 < -0.5)/n(), 3),
            pct.days.above.q5.0.5 = 100 - pct.days.below.q5.0.5,
            pct.days.below.q2.5.0.5 = 100*round(sum(!is.na(below.q2.5) & below.q2.5 < -0.5)/n(), 3),
            pct.days.above.q2.5.0.5 = 100 - pct.days.below.q2.5.0.5) %>%
  ungroup(depth, interval.yrs) %>% #subset(depth %in% c(0.06, 0.62, 1))
  mutate(depth.fac = factor(depth, levels = sort(unique(psi.stat.4$depth), decreasing = TRUE))) %>%
  mutate(int.ssn = paste(interval.yrs, season, sep = "_"))

heat.fr.drought.days.base <- ggplot(pct.drought.days %>% subset(interval.yrs.2 != "(Missing)"),
                   aes(x = interval.yrs.2, y = depth.fac)) +
  facet_wrap(~season, nrow = 2) +
  xlab("Census Interval") + ylab("Depth (m)") +
  scale_fill_viridis_c(expression('% Days '*Psi['Soil, DOY']*' < '*italic(Q)[italic(p)*'=0.025, '*Psi['Soil, DOY']]),
                       direction = -1, option = "plasma") +
  theme(legend.position = "top", legend.direction = "horizontal") +
  geom_tile(aes(fill = pct.days.below.q2.5)) +
  theme(axis.text.x = element_text(size = 10))
ggsave("pct.days.below.q2.5.clim_by depth & intervalyrs&ssn_full.jpeg",
       plot = heat.fr.drought.days.base, file.path(figures.folder), device = "jpeg", height = 6, width = 5, units='in')

heat.fr.drought.days.study <- heat.fr.drought.days.base %+%
  subset(pct.drought.days, interval.yrs != "(Missing)")
ggsave("pct.days.below.q2.5.clim_by depth & intervalyrs&ssn_study_period.jpeg",
       plot = heat.fr.drought.days.study, file.path(figures.folder), device = "jpeg", height = 6.5, width = 5, units='in')

heat.fr.drought.days.0.5 <- heat.fr.drought.days.base +
  geom_tile(aes(fill = pct.days.below.q2.5.0.5))
ggsave("pct.days.below.q2.5.clim.2.5_by depth & intervalyrs&ssn_full.jpeg",
       plot = heat.fr.drought.days.2.5, file.path(figures.folder), device = "jpeg", height = 6.5, width = 5, units='in')

## alpha = 0.05
heat.fr.drought.days.base.q5 <-  heat.fr.drought.days.base +
  geom_tile(aes(fill = pct.days.below.q5)) +
  scale_fill_viridis_c(expression('% Days '*Psi['Soil, DOY']*' < '*italic(Q)[italic(p)*'=0.05, '*Psi['Soil, DOY']]),
                     direction = -1, option = "plasma")
ggsave("pct.days.below.q5.clim_by depth & intervalyrs&ssn_full.jpeg",
       plot = heat.fr.drought.days.base.q5, file.path(figures.folder), device = "jpeg", height = 6, width = 5, units='in')
heat.fr.drought.days.study.q5 <- heat.fr.drought.days.base.q5 %+%
  subset(pct.drought.days, interval.yrs != "(Missing)")
ggsave("pct.days.below.q5.clim_by depth & intervalyrs&ssn_study_period.jpeg",
       plot = heat.fr.drought.days.study.q5, file.path(figures.folder), device = "jpeg", height = 6.5, width = 5, units='in')

heat.fr.drought.days.q5.0.5 <- heat.fr.drought.days.base.q5 +
  geom_tile(aes(fill = pct.days.below.q5.0.5))
ggsave("pct.days.below.q5.clim.2.5_by depth & intervalyrs&ssn_full.jpeg",
       plot = heat.fr.drought.days.q5.0.5, file.path(figures.folder), device = "jpeg", height = 6.5, width = 5, units='in')

rectangles.4 <- data.frame(
  xmin = 120,
  xmax = 335,
  ymin = 0,
  ymax = -3.0
)

## add legends for ribbon areas for drought rarity, lines for climatology,
# Could add a line for sd2
plot.psi.stat.6.interval.base <- plot.psi.stat.5.base %+%
  subset(psi.stat.4, depth %in% c(0.12, 0.62, 1)) +
  facet_wrap(. ~ interval.yrs.2) +
  geom_rect(data=rectangles.4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill='gray90', alpha=0.8) +
  geom_line(aes(x = doy, y = median.clim, group = as.factor(depth), color = as.factor(depth)), size = 0.3, linetype = "solid") +
  geom_ribbon(aes(x = doy, ymin = below.q5, ymax = median.clim, group = as.factor(depth_year),
                  fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE) +
    theme(panel.grid.major.y = element_line(size = 0.1)) +
  scale_color_discrete(name = "Depth (cm)", labels = c("10","60","100")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  coord_cartesian(ylim = c(-2, 0), xlim = c(0, 200)) +
  theme(legend.position = "top")
plot.psi.stat.6.interval.q10 <- plot.psi.stat.6.interval.base  +
  geom_ribbon(aes(x = doy, ymin = below.q10, ymax = median.clim, group = as.factor(depth_year),
                fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE)
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_q10.clim.jpeg",
       plot = plot.psi.stat.6.interval.q10, file.path(figures.folder), device = "jpeg", height = 3.25, width = 5, units='in')

plot.psi.stat.7.interval.q10 <- plot.psi.stat.6.interval.base %+%
  subset(psi.stat.4, depth %in% c(0.12, 0.62, 1) & interval.yrs != "(Missing)") +
  facet_wrap(. ~ interval.yrs, nrow = 1) +
  geom_ribbon(aes(x = doy, ymin = below.q10, ymax = median.clim, group = as.factor(depth_year),
                  fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE)
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_study_period_q10.jpeg",
       plot = plot.psi.stat.7.interval.q10, file.path(figures.folder), device = "jpeg", height = 2.5, width = 7, units='in')
# q5
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_q5.jpeg",
       plot = plot.psi.stat.6.interval.base, file.path(figures.folder), device = "jpeg", height = 3.25, width = 5, units='in')
plot.psi.stat.7.interval.q5 <- plot.psi.stat.6.interval.base %+%
  subset(psi.stat.4, depth %in% c(0.12, 0.62, 1) & interval.yrs != "(Missing)") +
  facet_wrap(. ~ interval.yrs, nrow = 1)
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_study_period_q5.jpeg",
       plot = plot.psi.stat.7.interval.q5, file.path(figures.folder), device = "jpeg", height = 2.5, width = 6, units='in')
# q2.5
plot.psi.stat.6.interval.q2.5 <- plot.psi.stat.6.interval.base  +
  geom_ribbon(aes(x = doy, ymin = below.q2.5, ymax = median.clim, group = as.factor(depth_year),
                  fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE)
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_q2.5.jpeg",
       plot = plot.psi.stat.6.interval.q2.5, file.path(figures.folder), device = "jpeg", height = 3.25, width = 5, units='in')
plot.psi.stat.7.interval.q2.5 <- plot.psi.stat.6.interval.base %+%
  subset(psi.stat.4, depth %in% c(0.12, 0.62, 1) & interval.yrs != "(Missing)") +
  facet_wrap(. ~ interval.yrs, nrow = 1) +
  geom_ribbon(aes(x = doy, ymin = below.q2.5, ymax = median.clim, group = as.factor(depth_year),
                  fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE)
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_study_period_q2.5.jpeg",
       plot = plot.psi.stat.7.interval.q2.5, file.path(figures.folder), device = "jpeg", height = 2.5, width = 6, units='in')

## plot over observed
plot.psi.stat.7.interval.median <- plot.psi.stat.6.interval.base %+%
  subset(psi.stat.4, depth %in% c(0.12, 0.62, 1) & interval.yrs.2 != "2015-2020") +
  facet_wrap(. ~ interval.yrs, nrow = 1) +
  geom_ribbon(aes(x = doy, ymin = median, ymax = median.clim, group = as.factor(depth_year),
                  fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE)
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_study_period_all_observed.jpeg",
       plot = plot.psi.stat.7.interval.median, file.path(figures.folder), device = "jpeg", height = 2.5, width = 7, units='in')

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
    geom_line(aes(x = doy, y = median.clim, linetype = "climatology", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
    geom_line(data = subset(psi.stat.4, year == year.on & depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)),
              aes(x = doy, y = median, linetype = "Year", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
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

### Plot pct.drought.days vs. mortality -------
###***************************************
pct.drought.days <- pct.drought.days %>%
  mutate(interval.num = as.numeric(as.character(recode(interval.yrs.2, `1982-1985` = "1",
                                 `1985-1990` = "2", `1990-1995` = "3",
                                 `1995-2000` = "4", `2000-2005` = "5",
                                 `2005-2010` = "6", `2010-2015` = "7", `2015-2018` = "8"))))
pct.drought.days.mean <- pct.drought.days %>%
  dplyr::select(-season, -int.ssn) %>%
  group_by(interval.num, interval.yrs, interval.yrs.2, depth, depth.fac) %>%
    summarise_all(list(~sum(., na.rm = TRUE)))

mrate.drought <- mrate.long %>% left_join(pct.drought.days.mean, by = "interval.num")  %>%
  left_join(bci.traits %>% dplyr::select(sp, form1), by = "sp")

mrate.drought.ssn <- mrate.long %>% left_join(pct.drought.days, by = "interval.num")  %>%
  left_join(bci.traits %>% dplyr::select(sp, form1), by = "sp")

mrate.drought.plot <- ggplot(mrate.drought %>%
                               subset(size == "large" & form1 == "T" &
                                        depth %in% c(0.12, 0.37, 1.00)),
                      aes(y = mrate, x = pct.days.below.q2.5)) +
  geom_point(aes(color = sp), show.legend = FALSE) + ylab(expression('Mortality Rate (% '*'year'^1*')')) +
  # xlab('% Drought Days (below 2.5% Quantile)') +
  xlab(expression('% Drought Days ('*Psi['z, DOY']*' < '*italic(Q)[italic(p)*'=0.025, '*Psi['Soil, DOY']]*')')) +
  geom_smooth(method = "loess") +
  facet_grid(depth  ~ .) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'loess',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by_pct.drought.days_by_depth.jpeg")),
       plot = mrate.drought.plot, height = 3, width = 9, units='in')

mrate.drought.plot.2 <- ggplot(mrate.drought %>%
                               subset(size == "large" & form1 == "T" &
                                        depth %in% c(0.12, 0.37, 1.00)),
                             aes(y = mrate, x = censusint.m)) +
  geom_point(aes(color = pct.days.below.q2.5)) +
  ylab(expression('Mortality Rate (% '*'year'^1*')')) + xlab('Census Interval') +
  guides(color = guide_legend(title = '% Drought Days (below 2.5% Quantile)')) +
  scale_color_continuous(trans = "reverse") +
  theme(legend.position = "top") +
  geom_smooth(method = "loess") +
  facet_grid(depth  ~ .) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'loess',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by_interval_colored_by_pct.drought.days_by_depth.jpeg")),
       plot = mrate.drought.plot.2, height = 3, width = 9, units='in')

###
mrate.drought.plot.3 <- ggplot(mrate.drought %>%
                               subset(size == "large" & form1 == "T" &
                                        depth %in% c(0.12, 0.37, 1.00)),
                             aes(y = mrate, x = pct.days.below.q2.5.0.5)) +
  geom_point(aes(color = sp), show.legend = FALSE) + ylab(expression('Mortality Rate (% '*'year'^1*')')) +
  # xlab('% Drought Days (below 2.5% Quantile)') +
  xlab(expression('% Drought Days ('*Psi['z, DOY']*' < '*italic(Q)[italic(p)*'=0.025, '*Psi['Soil, DOY']]*' & < 0.5 MPa)')) +
  geom_smooth(method = "loess") +
  facet_grid(depth  ~ .) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'loess',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by_pct.drought.days_below0.5_by_depth.jpeg")),
       plot = mrate.drought.plot.3, height = 3, width = 9, units='in')

mrate.drought.plot.4 <- ggplot(mrate.drought %>%
                                 subset(size == "large" & form1 == "T" &
                                          depth %in% c(0.12, 0.37, 1.00)),
                               aes(y = diff.mrate, x = censusint.m)) +
  geom_point(aes(color = pct.days.below.q2.5.0.5)) +
  geom_line(aes(group = sp, color = pct.days.below.q2.5.0.5)) +
  geom_hline(aes(yintercept = 0)) +
  ylab(y.label.1) + xlab('Census Interval') +
  guides(color = guide_legend(title = '% Drought Days (<2.5% Quantile & <0.5MPa)')) +
  scale_color_continuous(trans = "reverse") +
  theme(legend.position = "top") +
  geom_smooth(method = "loess") +
  facet_grid(depth  ~ .) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'loess',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by_interval_colored_by_pct.drought.days_below0.5_by_depth.jpeg")),
       plot = mrate.drought.plot.4, height = 4, width = 6.5, units='in')

###
mrate.drought.plot.5 <- ggplot(mrate.drought.ssn %>%
                               subset(size == "large" & form1 == "T" &
                                        depth %in% c(0.12, 0.37, 1.00)),
                             aes(y = mrate, x = pct.days.below.q2.5)) +
  geom_point(aes(color = sp), show.legend = FALSE) +
  ylab(expression('Mortality Rate (% '*'year'^1*')')) +
  # xlab('% Drought Days (below 2.5% Quantile)') +
  xlab(expression('% Drought Days ('*Psi['z, DOY']*' < '*italic(Q)[italic(p)*'=0.025, '*Psi['Soil, DOY']]*')')) +
  geom_smooth(method = "lm") +
  facet_grid(depth  ~ season) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.6, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by_pct.drought.days_by_depth_ssn.jpeg")),
       plot = mrate.drought.plot.5, height = 3, width = 9, units='in')

mrate.drought.plot.6 <- ggplot(mrate.drought.ssn %>%
                                 subset(size == "large" & form1 == "T" &
                                          depth %in% c(0.12, 0.37, 1.00)),
                               aes(y = diff.mrate, x = censusint.m)) +
  geom_point(aes(color = pct.days.below.q2.5)) +
  geom_line(aes(group = sp, color = pct.days.below.q2.5)) +
  geom_hline(aes(yintercept = 0)) +
  ylab(y.label.1) + xlab('Census Interval') +
  guides(color = guide_legend(title = '% Drought Days (below 2.5% Quantile)')) +
  scale_color_continuous(trans = "reverse") +
  theme(legend.position = "top") +
  geom_smooth(method = "loess") +
  facet_grid(depth  ~ season) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.6, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by_interval_colored_by_pct.drought.days_by_depth_ssn.jpeg")),
       plot = mrate.drought.plot.6, height = 4, width = 9, units='in')

