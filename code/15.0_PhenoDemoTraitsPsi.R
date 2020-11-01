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
#******************************************************
### Load data -------
#******************************************************
source("code/load.R")

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

## Function used to estimate K_leaf from Psi and species specific parameters
Exponential <- function (A, B, psi) {
  A * exp(-B * psi)
}

psi.corr.fun.ls.2 <- list(
  "gr.Psi" =
    function(df, dflc) {
      result.df <-
        psi.study[, psi.mod := range01(Exponential(A = df$A, B = df$B, psi = -psi))][
            , keyby = .(depth, interval, par.sam), .(gfac = mean(psi.mod, na.rm = TRUE))]#[
              # , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(list(result.df = result.df))
    },
  "gr.Psi.VPD.add" =
    function(df, dflc) {
      result.df <-
        psi.study[, psi.mod := Exponential(A = df$A, B = df$B, psi = -psi)][
          , keyby = .(depth, interval, par.sam), .(gfac = mean(c(psi.mod + VPD.effect), na.rm = TRUE))]#[
      # , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(list(result.df = result.df))
    },
  "gr.Psi.VPD.multi" =
    function(df, dflc) {
      result.df <-
        psi.study[, psi.mod := range01(Exponential(A = df$A, B = df$B, psi = -psi))][
          , keyby = .(depth, interval, par.sam), .(gfac = mean(psi.mod*std.VPD, na.rm = TRUE))]#[
            # , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(list(result.df = result.df))
    },
  "gr.Psi.leaf" =
    function(df, dflc) {
      # dflc.dt <- data.table(doy = dflc$doy, leaf_cover = dflc$leaf_cover)
      result.df <-
        as.data.table(psi.study)[data.table(dflc), on = 'doy'][,
          psi.mod := range01(Exponential(A = df$A, B = df$B, psi = -psi))][
            , keyby = .(depth, interval, par.sam), .(gfac = mean(psi.mod*leaf_cover, na.rm = TRUE))]#[
              # , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(list(result.df = result.df))
    },
  "gr.Psi.VPD.leaf.add" =
    function(df, dflc) {
      dflc.dt <- data.table(doy = dflc$doy, leaf_cover = dflc$leaf_cover)
      result.df <-
        as.data.table(psi.study)[dflc.dt, on = 'doy'][,
                                                      psi.mod := range01(Exponential(A = df$A, B = df$B, psi = -psi))][
                                                        , keyby = .(depth, interval, par.sam), .(gfac = mean(c(psi.mod + std.VPD)*leaf_cover, na.rm = TRUE))]#[
      # , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(list(result.df = result.df))
    },
  "gr.Psi.VPD.leaf.multi" =
    function(df, dflc) {
      dflc.dt <- data.table(doy = dflc$doy, leaf_cover = dflc$leaf_cover)
      result.df <-
        as.data.table(psi.study)[dflc.dt, on = 'doy'][,
        psi.mod := range01(Exponential(A = df$A, B = df$B, psi = -psi))][
          , keyby = .(depth, interval, par.sam), .(gfac = mean(c(psi.mod*std.VPD*leaf_cover), na.rm = TRUE))]#[
            # , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
      return(list(result.df = result.df))
  }#,
  ###---not used----
  # "gr.Psi.Rad.leaf.multi" =
  #   function(df, dflc) {
  #     dflc.dt <- data.table(doy = dflc$doy, leaf_cover = dflc$leaf_cover)
  #     result.df <-
  #       as.data.table(psi.study)[dflc.dt, on = 'doy'][,
  #                                                     psi.mod := range01(Exponential(A = df$A, B = df$B, psi = -psi))][
  #                                                       , keyby = .(depth, interval, par.sam), .(gfac = mean(c(psi.mod*std.Rs*leaf_cover), na.rm = TRUE))]#[
  #     # , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
  #     result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
  #     return(list(result.df = result.df))
  #   },
  # "gr.Psi.VPD.Rad.leaf.multi" =
  #   function(df, dflc) {
  #     dflc.dt <- data.table(doy = dflc$doy, leaf_cover = dflc$leaf_cover)
  #     result.df <-
  #       as.data.table(psi.study)[dflc.dt, on = 'doy'][,
  #                                                     psi.mod := range01(Exponential(A = df$A, B = df$B, psi = -psi))][
  #                                                       , keyby = .(depth, interval, par.sam), .(gfac = mean(c(psi.mod*std.VPD*std.Rs*leaf_cover), na.rm = TRUE))]#[
  #     # , keyby = .(depth, interval), .(gfac = mean(gfac, na.rm = TRUE))]
  #     result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "gfac")
  #     return(list(result.df = result.df))
  #   }
  ###---not used----
)

###************************************************************
### Functions explored earlier but not used anymore : psi.corr.fun.ls
###************************************************************

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
  "mr.kl20.I" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := indicator(psi, df$psi_kl20, greater.than = FALSE)][
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
    },
  "mr.kl20.I.VPD" =
    function(df) {
      result.df <-
        as.data.table(psi.study)[, psi.mod := indicator(psi, df$psi_kl20, greater.than = FALSE)][
          , keyby = .(depth, interval, par.sam), .(mfac = sum(psi.mod*std.VPD, na.rm = TRUE))][
            , keyby = .(depth, interval), .(mfac = mean(mfac, na.rm = TRUE))]
      result.df <- data.frame(result.df) %>% pivot_wider(names_from = "depth", values_from = "mfac")
      return(result.df)
    }
)
# Originally written by Sean M McMahon, then modified

get.ts.lk <- function(df) {
  demo.ts <- df$demo.rate
  k <- length(demo.ts)
  psi.dat <- df[, grep("[0-9]", names(df))]
  psi.dat <- psi.dat[, colSums(!is.na(psi.dat)) > 0]
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
  dplyr::select(-sp_size, -"_", -mean, -sd, -trees) %>%
  # remove species taht are not canopy
  subset(sp %in% bci.traits$sp[bci.traits$form1 == "T"])
length(unique(growth.sub$sp))

grate.plot <- ggplot(growth.sub, aes(x = interval, y = demo.rate)) +
  geom_line(aes(group = sp, color = sp), show.legend = FALSE) +
  facet_wrap(. ~ sp) +
  ylab("Std. Growth") + xlab("Interval") +
  geom_errorbar(aes(ymax = demo.rate + se, ymin = demo.rate - se), width = 0.1, size = 0.5)
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
         # temporarily renaming VPD variable name for the model input to work
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
                     std.Rs, std.pet.PM, std.VPD, VPD.effect), by = "date")
depth.sub <- c(soil.depths[5:length(soil.depths)])
depth.breaks <- c(soil.depths[5], soil.depths[7:length(soil.depths)])
depth.labels <- c(0.4, soil.depths[8:length(soil.depths)]) # mean(c(0.21, 0.37, 0.62))

#cut(depth.sub, include.lowest = TRUE, breaks = depth.breaks, labels = depth.labels, right = TRUE)
psi.study <- as.data.table(psi.m)[!is.na(interval),][,
  ':='(doy = format(date, "%j"), year = format(date, "%Y"))][
  (doy < 368) & (!year %chin% c("1990")) &
    (!year %chin% c("1991") | depth == 2.9) &
    (!year %chin% c("1991", "1992") | depth < 2.9)][,
     c("year") := NULL][,
      doy := as.numeric(as.character(doy))][
        depth %chin% depth.sub][,':='(depth = cut(depth, include.lowest = TRUE,
                                                  breaks = depth.breaks, labels = depth.labels, right = TRUE))][,
                                                                    by = .(date, interval, interval.yrs, par.sam, depth), lapply(.SD, mean, na.rm = TRUE)]

##_______________________________________________________
## Note above: There are lines above at the far right end
##_______________________________________________________

##*******************************************
## Save soil Psi for plotting-----
## Not ERD model dependent, so can be turned off for ERD tests
##*******************************************

# psi.stat.4 <- psi %>%
#   mutate(doy = format(date, "%j"),
#          year = format(date, "%Y")) %>%
#   ## due to time nneded for model initialisation soil layers deeper than 2.9 m
#   ## are not fully recharged until the end of 1992, layer 2.9 m recharges by the end of 1992
#   subset(!year %in% c("1990")) %>%
#   subset(!year %in% c("1990", "1991") | depth == 2.9) %>%
#   subset(!year %in% c("1990", "1991", "1992") | depth < 2.9) %>%
#   group_by(date, doy, year, interval.yrs, depth) %>%
#   summarise(median = median(psi, na.rm = TRUE),
#             q100 = quantile(psi, probs = 1),
#             q97.5 = quantile(psi, probs = 0.975),
#             q2.5 = quantile(psi, probs = 0.025),
#             q5 = quantile(psi, probs = 0.05), .groups = "drop_last") %>%
#   ungroup(doy, year, depth) %>%
#   mutate(doy = as.numeric(doy))
# psi.stat.5 <- psi.stat.4 %>%
#   group_by(doy, depth) %>%
#   summarise(median.clim = median(median, na.rm = TRUE),
#             q100.clim = quantile(median, probs = 1),
#             q97.5.clim = quantile(median, probs = 0.975),
#             q2.5.clim = quantile(median, probs = 0.025),
#             q10.clim = quantile(median, probs = 0.1),
#             q5.clim = quantile(median, probs = 0.05), .groups = "drop_last") %>%
#   ungroup(doy, depth) %>%
#   mutate(doy = as.numeric(doy))
# rectangles.3 <- data.frame(
#   xmin = 120,
#   xmax = 335,
#   ymin = 0,
#   ymax = -2.5
# )
# ## For selected depths used in inverse modeling
# # depth.sub, depth.breaks, depth.labels are defined above in section:
# psi.stat.4.select <- psi %>%
#   mutate(doy = format(date, "%j"),
#          year = format(date, "%Y")) %>%
#   ## due to time nneded for model initialisation soil layers deeper than 2.9 m
#   ## are not fully recharged until the end of 1992, layer 2.9 m recharges by the end of 1992
#   subset(!year %in% c("1990")) %>%
#   subset(!year %in% c("1990", "1991") | depth == 2.9) %>%
#   subset(!year %in% c("1990", "1991", "1992") | depth < 2.9) %>%
#   subset(depth %in% depth.sub) %>%
#   mutate(depth = cut(depth, include.lowest = TRUE, breaks = depth.breaks,
#                      labels = depth.labels, right = TRUE)) %>%
#   group_by(date, doy, year, interval.yrs, depth) %>%
#   summarise(median = median(psi, na.rm = TRUE),
#             q100 = quantile(psi, probs = 1),
#             q97.5 = quantile(psi, probs = 0.975),
#             q2.5 = quantile(psi, probs = 0.025),
#             q5 = quantile(psi, probs = 0.05), .groups = "drop_last") %>%
#   ungroup(doy, year, depth) %>%
#   mutate(doy = as.numeric(doy))
#
# psi.stat.5.select <- psi.stat.4.select %>%
#   group_by(doy, depth) %>%
#   summarise(median.clim = median(median, na.rm = TRUE),
#             q100.clim = quantile(median, probs = 1),
#             q97.5.clim = quantile(median, probs = 0.975),
#             q2.5.clim = quantile(median, probs = 0.025),
#             q10.clim = quantile(median, probs = 0.1),
#             q5.clim = quantile(median, probs = 0.05), .groups = "drop_last") %>%
#   ungroup(doy, depth) %>%
#   mutate(doy = as.numeric(doy))
#
# psi.stat.4 <- psi.stat.4 %>%
#   mutate(interval.yrs.2 = forcats::fct_explicit_na(cut(date, include.lowest = TRUE, breaks = c(cut.breaks, max(psi.stat.4$date, na.rm = TRUE)),
#                                                        labels = c(cut.labels.2, "2015-2018"), right = TRUE))) %>%
#   left_join(psi.stat.5, by = c("doy", "depth")) %>%
#   mutate(below.q10 = ifelse(median < q10.clim, median, NA),
#          below.q5 = ifelse(median < q5.clim, median, NA),
#          below.q2.5 = ifelse(median < q2.5.clim, median, NA),
#          depth_year = paste(depth, year, sep = "_"))
#
# psi.stat.4.select <- psi.stat.4.select %>%
#   mutate(interval.yrs.2 = forcats::fct_explicit_na(cut(date, include.lowest = TRUE, breaks = c(cut.breaks, max(psi.stat.4$date, na.rm = TRUE)),
#                                                        labels = c(cut.labels.2, "2015-2018"), right = TRUE))) %>%
#   left_join(psi.stat.5.select, by = c("doy", "depth")) %>%
#   mutate(below.q10 = ifelse(median < q10.clim, median, NA),
#          below.q5 = ifelse(median < q5.clim, median, NA),
#          below.q2.5 = ifelse(median < q2.5.clim, median, NA),
#          depth_year = paste(depth, year, sep = "_")) %>%
#   mutate(freq.below.q10 = ifelse(median < q10.clim, 1, NA),
#             freq.below.q5 = ifelse(median < q5.clim, 1, NA),
#             freq.below.q2.5 = ifelse(median < q2.5.clim, 1, NA))
#
# psi.stat.4.select.freq <- psi.stat.4.select %>%
#   group_by(interval.yrs, doy, depth) %>%
#   summarise(freq.below.q10 = sum(freq.below.q10, na.rm = TRUE),
#             freq.below.q5 = sum(freq.below.q5, na.rm = TRUE),
#             freq.below.q2.5 = sum(freq.below.q2.5, na.rm = TRUE), .groups = "drop") %>%
#   mutate(depth = as.numeric(as.character(depth)))
#
# ## this is the date that's shown in the graph, so only interested in shwoing an extreme year
# ## defined for the shown doy range
# xlim.in.wet.season <- 200
# psi.stat.4.extreme.yr  <- psi.stat.4 %>%
#   subset(doy < xlim.in.wet.season) %>%
#   group_by(depth_year) %>%
#   dplyr::summarise(extreme.yr.q2.5 = if_else(any(!is.na(below.q2.5)), TRUE, FALSE), .groups = "drop")
#
# psi.stat.4 <- psi.stat.4 %>%
#   left_join(psi.stat.4.extreme.yr, by = "depth_year")
#
# psi.stat.4 <- psi.stat.4 %>% mutate(plot.depth = paste0(round(depth, 1), " m"))
#
# save(psi.stat.4, file = file.path(results.folder, "psi.stat.4.Rdata"))
# save(psi.stat.4.select, file = file.path(results.folder, "psi.stat.4.select.Rdata"))
# save(psi.stat.4.select.freq, file = file.path(results.folder, "psi.stat.4.select.freq.Rdata"))
# save(psi.stat.5, file = file.path(results.folder, "psi.stat.5.Rdata"))
# save(psi.stat.5.select, file = file.path(results.folder, "psi.stat.5.select.Rdata"))
#
# # psi.study <- as.data.table(psi.m)[!is.na(interval),][,
# #   ':='(doy = format(date, "%j"),
# #        year = format(date, "%Y"))][
# #   (doy < 368) & (!year %chin% c("1990")) &
# #     (!year %chin% c("1991") | depth == 3) &
# #     (!year %chin% c("1991", "1992") | depth < 3)][, c("year") := NULL][,
# #     doy := as.numeric(as.character(doy))]
#
# # psi.study <- as.data.frame(psi.study)
##*******************************************

data.model.AB.sub <- data.model.AB %>%
  ## Using parameters fitted to data when available
  mutate(A = ifelse(is.na(data.A), model.A, data.A),
         B = ifelse(is.na(data.B), model.B, data.B)) %>%
  subset(sp %in% union(growth.sub$sp, mort.sub$sp)) %>%
  # ## And LMA from LMALAM from LEAF and DISC
  subset(sp.LMA.sub == "original") %>%
  ## removing A & B predicted from gap-filled WSG
  subset(sp.WSG.sub == "original")

save(data.model.AB.sub, file = file.path(results.folder, "data.model.AB.sub.Rdata"))

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
#   ## But "ficutr", "pri2co", "termob" do not have leaf_cover data and also not Evergreen either but DB,
#   ## so taking average DB species' leaf-cover,  "pri2co" 's deciduousness not known, so can't gap-fill
#   # First joinging average deci leaf_cover then substituing for missing data by species
# sp.leaf_cover.for.model <- sp.leaf_cover.for.model %>%
#   left_join(sp.leaf_cover.for.model %>% group_by(deciduous, doy) %>%
#               summarise(deci.leaf_cover = mean(leaf_cover, na.rm = TRUE), .groups = "drop"),
#             by = c("deciduous", "doy")) %>%
#   mutate(leaf_cover = ifelse(is.na(leaf_cover), deci.leaf_cover, leaf_cover))

setdiff(unique(growth.sub$sp), unique(sp.leaf_cover.for.model$sp))
# unique(sp.leaf_cover.for.model[is.na()]$sp)

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
names(ml.dens) <- names(ml.rsq.ls) <- names(ml.corr.ls) <- names(ml.ls) <- names.gfac
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
                       f = list(demo.psi$sp, demo.psi$par.sam), drop = TRUE)
  ## Get for each sp-par.sam get AIC scores <= min(AIC) + 2 and corresponding depths, ml, R2 and corr
  ml.ls[[i]] <- lapply(demo.psi.ls, get.ts.lk)
  ## Retain the one with maximum likelihood
  ml.dens[[i]] <- sapply(ml.ls[[i]], get.ml.max)
  ## Get corresponding R, corr and depth
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
ml.corr.best.parsam <- vector("list", length(gfac.interval))
names(ml.corr) <- names(ml.corr.best) <- names(ml.corr.best.parsam) <- names.gfac
for (i in names(gfac.interval)) {
  ml.corr[[i]] <- ml.corr.ls[[i]] %>%
    mutate(sp.parsam.row = row.names(ml.corr.ls[[i]])) %>%
    separate(sp.parsam.row, c("sp", "par.sam", "row", sep = "."), remove = FALSE, extra = "drop", fill = "right") %>%
    dplyr::select(-".", -row) %>%
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
  # Select best corr within each sp and par.sam
  ml.corr.best.parsam[[i]] <- ml.corr[[i]] %>%
    group_by(sp, par.sam) %>%
    mutate(corr.max = max(corr, na.rm = TRUE)) %>%
    subset(corr == corr.max) %>% dplyr::select(-corr.max) %>%
    ungroup(sp, par.sam)

  # Get sumamy corr across the best corr within each sp and par.sam
  ml.corr.best[[i]] <- ml.corr.best.parsam[[i]] %>% group_by(sp) %>%
    summarise(corr.upper.CI = -quantile(corr, probs = 0.975),
              corr.lower.CI = -quantile(corr, probs = 0.025),
              corr = -median(corr, na.rm = TRUE), .groups = "drop_last")
}
## species by depth by heat rsq
ml.rsq.combine <- dplyr::bind_rows(ml.corr, .id = "corr.func") %>%
  transform(corr.func = factor(corr.func, levels = names.gfac))
ml.rsq.combine.best.parsam <- dplyr::bind_rows(ml.corr.best.parsam, .id = "corr.func") %>%
  transform(corr.func = factor(corr.func, levels = names.gfac)) %>%
  unite(corr.func_sp_depth, corr.func, sp, depth, remove = FALSE)
ml.rsq.combine.best <- ml.rsq.combine.best.parsam %>%
  group_by(sp, corr.func, size, deciduous, deciduousness, deciduousness.label, DeciLvl, sp.plot, deci_sp, deci_sp.plot) %>%
  # subset(corr >= 0.7071068) %>% # that is R2 >= 0.5 and corr >= 0
  subset(corr >= 0 & R2 >= 0.5) %>%
  summarise(depth.se = sd(depth, na.rm = TRUE)/sqrt(n()),
            depth = median(depth, na.rm = TRUE),
            corr.se = sd(corr, na.rm = TRUE)/sqrt(n()),
            corr = mean(corr, na.rm = TRUE),
            R2.se = sd(R2, na.rm = TRUE)/sqrt(n()),
            R2 = mean(R2, na.rm = TRUE), .groups = "drop_last") %>%
  dplyr::select(sp, depth, corr, R2, everything()) %>%
  unite(corr.func_sp_depth, corr.func, sp, depth, remove = FALSE)

###____________________________
## Compare with isotope data----
###____________________________
iso.2 <- iso.2.raw %>%
  # subset(source == "Meinzer et al.1999 Fig. 5A") %>%
  mutate(source = "Meinzer et al.1999 Fig. 5B", location = "BCI") %>%
  group_by(sp, source) %>%
  summarise(se = sd(Xylem_sap_deltaD_permil, na.rm = TRUE)/sqrt(n()),
            Xylem_sap_deltaD_permil = mean(Xylem_sap_deltaD_permil, na.rm = TRUE),
            n = n(),
            DBH = mean(DBH, na.rm = TRUE), .groups = "drop_last")

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
                        se.mean = mean(se, na.rm = TRUE), .groups = "drop_last") %>%
              dplyr::select(sp, Xylem_sap_deltaD_permil.mean, se.mean), by = "sp") %>%
  droplevels()

ml.rsq.combine.sub.iso <- ml.rsq.combine.best %>%
  mutate(depth = as.numeric(depth)) %>%
  subset(!sp %in% c("guapst") & !is.na(Xylem_sap_deltaD_permil.mean)) %>%
  left_join(bci.traits %>%
              dplyr::rename(Code = sp, Genus = GENUS., Species = SPECIES., Family = FAMILY.) %>%
              mutate(s.names = paste0(substr(Genus, start = 1, stop = 1), ". ", tolower(Species)),
                     sp = tolower(Code)) %>%
              dplyr::select(sp, s.names), by = "sp") %>%
  subset(source == "Meinzer et al.1999 Fig. 4") %>%
  droplevels()

depth.rsq.isotopes <- ml.rsq.combine.best %>%
  group_by(corr.func, sp, size) %>%
  summarise_at(c("depth", "depth.se","Xylem_sap_deltaD_permil", "se"), mean, na.rm = TRUE) %>%
  ungroup(corr.func, sp, size)

save(depth.rsq.isotopes, file = file.path(results.folder, "depth.rsq.isotopes.Rdata"))

erd.model.n.sp <- erd.model.iso.n.sp <- erd.model.p <- erd.model.r2 <- vector("numeric", length(names.gfac))
names(erd.model.n.sp) <- names(erd.model.iso.n.sp) <- names(erd.model.p) <- names(erd.model.r2) <- names.gfac
for(i in 1: length(names.gfac)) {
  erd.sp.data <- ml.rsq.combine.best %>%
    subset(corr.func == names.gfac[i]) %>%
    subset(!duplicated(sp) & !is.na(depth)) %>% droplevels()
  erd.model.n.sp[[i]] <- length(erd.sp.data$sp[!is.na(erd.sp.data$depth)])

  m.data <- ml.rsq.combine.sub.iso %>%
    # subset(sp != "guapst") %>%
    subset(corr.func == names.gfac[i])
  if(nrow(m.data) > 0) {
    lm.model <- lm(depth ~ Xylem_sap_deltaD_permil, data = m.data)
    if(is.na(summary(lm.model)$coefficients[2])) {
      erd.model.p[[i]] <- NA; erd.model.r2[[i]] <- NA
    } else {
      erd.model.p[[i]] <- lmp(lm.model)
      erd.model.r2[[i]] <- summary(lm.model)$adj.r.squared
      erd.model.iso.n.sp[[i]] <- length(m.data$sp[!is.na(m.data$depth)])
    }
  }
}
erd.model.n.sp
erd.model.iso.n.sp
erd.model.p
erd.model.r2
## among the models for which p-value =< 0.05, chose the one that has the highest R2
# erd.model.sig <- which(erd.model.p == 0.05 | erd.model.p < 0.05)
# chosen.model <- names(erd.model.p)[which(erd.model.r2 == max(erd.model.r2[erd.model.sig], na.rm = TRUE))]

chosen.model <- "gr.Psi.VPD.multi"
save(chosen.model, file = file.path(results.folder, "chosen.model.Rdata"))
save(erd.model.n.sp, file = file.path(results.folder, "erd.model.n.sp.Rdata"))
save(erd.model.p, file = file.path(results.folder, "erd.model.p.Rdata"))
## Plot chosen ERD
df.erd.to.plot <- ml.rsq.combine.best %>%
  subset(corr.func == chosen.model) %>%
  mutate(size = as.character(size)) %>%
  subset(!duplicated(sp) & !is.na(depth)) %>% droplevels() %>%
  transform(sp = reorder(sp, depth))
length(unique(df.erd.to.plot$sp))
# 36

save(ml.rsq.combine.best, file = file.path(results.folder, "ml.rsq.combine.best.Rdata"))
save(ml.rsq.combine, file = file.path(results.folder, "ml.rsq.combine.Rdata"))
save(df.erd.to.plot, file = file.path(results.folder, "df.erd.to.plot.Rdata"))

###____________________________
## Compare with stem hydraulics data----
###____________________________

hyd.mod <- hyd %>% left_join(depth.rsq.isotopes %>% ungroup() %>%
                               subset(corr.func == chosen.model) %>%
                           dplyr::select(sp, depth, depth.se), by = "sp") %>%
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
                                                  expression(Psi[Predawn]), expression(Psi[min]), expression(Psi[tlp]),
                                                  expression(italic('K')['max, stem']),  expression(Psi['50, stem']),  expression(Psi['88, stem']),
                                                  expression(Psi[min]*' - '*Psi[tlp]),
                                                  expression(Psi[min]*' - '*Psi['50, stem']),
                                                  expression(Psi[min]*' - '*Psi['88, stem']),
                                                  expression(Psi[tlp]*' - '*Psi['50, stem']),
                                                  expression(Psi[tlp]*' - '*Psi['88, stem']),
                                                  expression('CWR'['total']), expression('CWR'['xylem']), expression('CWR'['bark']),
                                                  expression(italic('F')['Elbow, xylem']), expression(italic('F')['Elbow, bark']),
                                                  expression(Psi[min]*'-'*italic('F')['elbow, xylem']),
                                                  expression(Psi[min]*'-'*italic('F')['elbow, bark']),
                                                  expression(italic('F')['cap, xylem']), expression(italic('F')['cap, bark']),
                                                  expression('WD'[stem]),
                                                  expression('Panama'[wet]), expression('Plot'[wet]), "LMA")),
            trait.plot.chart = factor(trait, labels = c(expression(Depth[italic('Rsq')]), expression(italic(delta)^2*H[Xylem]),
                                                        expression(Psi[Predawn]), expression(Psi[min]), expression(Psi[tlp]),
                                                        expression(italic('K')['max, stem']), expression(Psi['50,stem']),  expression(Psi['88,stem']),
                                                        expression(Psi[min]*'-'*Psi[tlp]),
                                                        expression(Psi[min]*'-'*Psi['50,stem']),
                                                        expression(Psi[min]*'-'*Psi['88,stem']),
                                                        expression(Psi[tlp]*'-'*Psi['50,stem']),
                                                        expression(Psi[tlp]*'-'*Psi['88,stem']),
                                                        expression('CWR'['total']), expression('CWR'['Xylem']), expression('CWR'['bark']),
                                                        expression(italic('F')['elbow, xylem']), expression(italic('F')['Elbow,bark']),
                                                        expression(Psi[min]*'-'*italic('F')['Elbow,xylem']),
                                                        expression(Psi[min]*'-'*italic('F')['Elbow,bark']),
                                                        expression(italic('F')['cap,xylem']), expression(italic('F')['cap,bark']),
                                                        expression('WD'[stem]),
                                                        expression('Panama'[wet]), expression('Plot'[wet]), "LMA"))) %>% droplevels()

hyd.long <- hyd.mod %>% dplyr::select(-DeciLvl) %>%
  dplyr::select(sp, deciduousness, deciduous, location, TLP, KmaxS, p50S, p88S,
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
              summarise(value = max(value, na.rm = TRUE)), by = c("trait"), .groups = "drop_last") %>%
  subset(deciduousness != "NA") %>%
  droplevels() %>%
  left_join(traits.labels.table.1 %>% dplyr::select(trait, trait.plot), by = "trait")

hyd.long <- hyd.long %>%
  left_join(traits.labels.table.1, by = "trait")

hyd.error <- hyd %>% dplyr::select(sp, KmaxS_se, vc_b_se, vc_a_se, tlp_sd) %>%
  rename(KmaxS = KmaxS_se, vc_b = vc_b_se, vc_a = vc_a_se, TLP = tlp_sd) %>%
  gather(trait, se, -sp) ## But note that for TLP it's not se but sd

hyd.depth <- hyd.long %>%
  subset(trait == "depth") %>%
  dplyr::select(sp, deciduousness, trait.plot.chart, value) %>%
  pivot_wider(names_from = trait.plot.chart, values_from = value)

erd.stem.traits <- hyd.long %>%
  subset(trait %in% c("KmaxS", "TLP", "p88S", "HSM88S", "lwp.min_Diurnal", "lwp.min_Predawn")) %>%
  subset(trait != "depth") %>%
  left_join(hyd.depth %>% dplyr::select(sp, `Depth[italic("Rsq")]`), by = "sp") %>%
  dplyr::select(deci_sp, sp, trait, `Depth[italic("Rsq")]`, value) %>%
  # bind_rows(depth.traits.kunert %>% subset(trait == "KmaxL") %>%
  #             select(deci_sp, sp, trait, `Depth[italic("Rsq")]`, value)) %>%
  left_join(hyd.error, by = c("sp", "trait"))

save(erd.stem.traits, file = file.path(results.folder, "erd.stem.traits.Rdata"))
save(hyd.long, file = file.path(results.folder, "hyd.long.prepped.Rdata"))


##______________________________________________________________________
## Plot mortality by time spent below a threshold in the preferred depth-------
# For each sp-depth calculate number of days below a threshold with an indicator function
##______________________________________________________________________

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

save(mfac.interval.long, file = file.path(results.folder, "mfac.interval.long.Rdata"))
## Ordered along Rooting Depth Index
# load(file = file.path(results.folder, "mfac.interval.long.Rdata"))
mfac.on <- "mr.kl50.I"
mrate.depth <-
  adult.mrate.long %>% mutate(size = "large") %>%
  # mrate.long %>%
  left_join(subset(depth.rsq.isotopes, corr.func == chosen.model) %>%
              rename(rdi.gr = depth) %>%
              dplyr::select(sp, rdi.gr, depth.se) %>%
              mutate(size = "large"), by = c("sp", "size")) %>%
  left_join(bci.traits %>% dplyr::select(form1, sp), by = "sp") %>%
  # mutate(sp.plot = factor(sp, levels = unique(sp[order(rdi.gr)]), ordered = TRUE)) %>%
  mutate(size = as.character(size)) %>%
  subset(size == "large" & form1 == "T") %>% droplevels()
mrate.mfac.depth <- mrate.depth %>%
  right_join(mfac.interval.long[[mfac.on]], #%>%
             # mutate(censusint.m = recode(interval.num, `1` = "1982-85", `2` = "1985-90",
             #                             `3` = "1990-95", `4` = "1995-00", `5` = "2000-05", `6` = "2005-10", `7` = "2010-15")),
             by = c("interval.num", "sp", "size")) %>%
  left_join(psi.study %>%
              subset(depth == unique(psi.study$depth)[1] & par.sam == unique(psi.study$par.sam)[1]) %>%
              mutate(censusint.m = recode(interval, `1` = "1981-85", `2` = "1985-90",
                                          `3` = "1990-95", `4` = "1995-00", `5` = "2000-05", `6` = "2005-10", `7` = "2010-15")) %>%
              group_by(censusint.m) %>%
              summarise(days = length(date), .groups = "drop"), #%>%
            # mutate(interval = as.numeric(interval)),
            by = c("censusint.m")) %>%
  mutate(mfac.rate = mfac*100/days) %>% # expressed in days per year
  mutate(sp_size = paste(sp, size, sep = "_")) %>%
  group_by(sp, size, censusint.m) %>%
  mutate(mfac.soil.column = sum(mfac, na.rm = TRUE)) %>%
  ungroup(sp, censusint.m)

save(mrate.depth, file = file.path(results.folder, "mrate.depth.Rdata"))
save(mrate.mfac.depth, file = file.path(results.folder, "mrate.mfac.depth.Rdata"))

