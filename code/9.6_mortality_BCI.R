rm(list=ls())
gc()
# load("/Library/Frameworks/R.framework/Versions/3.4/Resources/library/CTFSRPackage/CTFSRPackage.Rdata")
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, hms, ggpmisc)

# graphics info
theme_set(theme_bw())
theme_update(text = element_text(size=14),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)

census.years <- c(1982, 1985, 1990, 1995, 2000, 2005, 2010, 2015)

require(scales)
rev_sqrt_trans <- function() {
  scales::trans_new(
    name = "rev_sqrt",
    transform = function(x) -sqrt(abs(x)),
    inverse = function(x) x^2);
}

####********************************************************************
##### Load demographic data ----------------------------------------
##### If you dont need to update the above, start directly here:
####********************************************************************

demo.graphs <- function(level.folder = level.folder, n.threshold = n.threshold,
                        dbh.residuals = dbh.residuals,
                        dryseason = dryseason,
                        root.selection = root.selection,
                        iso.subset = iso.subset) {

  load("results/demo.sp.RData")
  load("results/demo.sp_size.RData")
  load("results/mrate.long.RData")
  load("results/sp.mrate.long.RData")
  load("results/adult.mrate.long.RData")

  ##plotting mortality rates by sp
  census.years.short <- format(strptime(census.years, "%Y"), "%y")
  # mrate.long$censusint <- mapvalues(mrate.long$census, from = census.years[-1], to =  paste(census.years[-length(census.years)],census.years[-1], sep = "-"))
  # mrate.long$censusint.s <- mapvalues(mrate.long$census, from = census.years[-1], to =  paste(census.years.short[-length(census.years)],census.years.short[-1], sep = "-"))
  mrate.long$censusint.m <- recode(mrate.long$census, `1985` = "1982-85", `1990` = "1985-90", `1995` = "1990-95", `2000` = "1995-00", `2005` = "2000-05", `2010` = "2005-10", `2015` = "2010-15")
  sp.mrate.long$censusint.m <- recode(sp.mrate.long$census, `1985` = "1982-85", `1990` = "1985-90", `1995` = "1990-95", `2000` = "1995-00", `2005` = "2000-05", `2010` = "2005-10", `2015` = "2010-15")
  adult.mrate.long$censusint.m <- recode(adult.mrate.long$census, `1985` = "1982-85", `1990` = "1985-90", `1995` = "1990-95", `2000` = "1995-00", `2005` = "2000-05", `2010` = "2005-10", `2015` = "2010-15")

  ## adding uptake depth index
  load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
  load(file.path("results/4.1GLUEsetup_part2_BCI.RData")) # has n.ensembles and growth and si matrix

  intervals <- info$intervals
  n.ensembles <- growth_by_si.info$n.ensembles
  growth.type <- growth_by_si.info$growth.type
  growth.selection <- growth_by_si.info$growth.selection
  si.type <- growth_by_si.info$si.type
  goodness.fit <- 0.3
  soil.depths <- unique(info$root.param.long$depth)
  ##
  ##
  file.extension.base4 <- paste0(goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection,
                                 "_", dbh.residuals, "_", intervals, "_dryseason_", dryseason, "_iso.subset_",
                                 iso.subset, "_root.selection_", root.selection)
  load(file = paste0("results/", level.folder, "/ds.bestfit_cor", file.extension.base4, ".Rdata"))

  figures.folder <- paste0("figures/mortality/", growth.type, "/", level.folder, "/", paste0("dryseason_", dryseason, "_dbh.residuals_", dbh.residuals, "_iso.subset_",
                                                                                             iso.subset, "_root.selection_", root.selection, "_avg.abund_above_", n.threshold))
  if(!dir.exists(figures.folder)) {dir.create(figures.folder)}

  ds <- ds.bestfit
  # load(file = paste("results/commlevel/ds.bestfit_cor", file.extension.base4, ".Rdata", sep = ""))
  # ds <- rbind(ds, ds.bestfit)
  ds <- ds %>% mutate(tlplevel = as.factor(tlplevel)) %>% subset(!is.na(udi)) %>% droplevels()

  head(ds)
  udi <- subset(ds, select = c("sp", "sp_size", "size", "udi.best", "tlplevel")) # best.type.rsq > 0.7

  head(udi)
  nrow(udi)

  deci <- read.csv("data-raw/traits/HydraulicTraits_Kunert/deciduous_species_Meakem.csv")
  deci.level_key <- c("Evg" = 1, "DF" = 2, "DB" = 3, "DO" = 4, "D" = "4") #c(a = "apple", b = "banana", c = "carrot")
  deci <- deci %>% mutate(sp = as.character(Species.code), Deciduousness = as.character(Deciduousness)) %>%
    select(sp, Deciduousness)

  head(sp.mrate.long)
  sp.mrate.long <- left_join(sp.mrate.long, subset(udi, size == "large") %>% select(sp, udi.best), by = "sp") %>%
    mutate(sp = as.character(sp)) %>%
    transform(sp = reorder(sp, -udi.best)) %>%
    left_join(deci, by = "sp") %>%
    mutate(Deciduousness = replace_na(Deciduousness, "Evg")) %>%
    mutate(DeciLvl = as.numeric(recode_factor(Deciduousness, !!!deci.level_key)))

  adult.mrate.long <- left_join(adult.mrate.long, subset(udi, size == "large") %>% select(sp, udi.best), by = "sp") %>%
    mutate(sp = as.character(sp)) %>%
    transform(sp = reorder(sp, -udi.best)) %>%
    left_join(deci, by = "sp") %>%
    mutate(Deciduousness = replace_na(Deciduousness, "Evg")) %>%
    mutate(DeciLvl = as.numeric(recode_factor(Deciduousness, !!!deci.level_key)))

  mrate.long <- left_join(mrate.long, select(udi, sp_size, udi.best), by = "sp_size") %>%
    mutate(sp_size = as.character(sp_size)) %>%
    transform(sp_size = reorder(sp_size, -udi.best)) %>%
    separate(sp_size, c("sp", "size", sep = "_"), remove = FALSE, extra = "drop", fill = "right") %>%
    select(-"_") %>%
    left_join(deci, by = "sp") %>%
    mutate(Deciduousness = replace_na(Deciduousness, "Evg")) %>%
    mutate(DeciLvl = as.numeric(recode_factor(Deciduousness, !!!deci.level_key)))


  # ggplot(mrate.long, aes(x = censusint.m, y = mrate)) +
  #   geom_line(aes(x = censusint.m, y = mrate, group = sp, color = sp), size = 1, show.legend = F) +
  #   my.theme + my.bg + my.adjust + theme(panel.grid.major.x = element_line(colour = "grey", size = 0.01)) +
  #   ylab(expression("Mortality Rate")) + xlab("Census Interval")
  # ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_Mortality_rate_above", dbhthreshold/10, "cmDBH_aboveN50.jpeg")), height = 5, width = 9, units='in')
  #
  # ggplot(mrate.long, aes(x = censusint.m, y = mrate)) +
  #   geom_line(aes(x = censusint.m, y = mrate, group = sp, color = sp), size = 1, show.legend = F) +
  #   geom_point(aes(x = censusint.m, y = mrate, color = sp), size = 3, show.legend = F) +
  #   facet_wrap( ~ sp, scales = "free_y") +
  #   theme(axis.text.y = element_text(size = 8),
  #         axis.text.x = element_text(size = 6, face = "plain", angle = 45, vjust = 1, hjust = 1)) +
  #   ylab(expression("Mortality Rate")) + xlab("Census Interval")
  # ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_Mortality_rate_by_period_above", dbhthreshold/10, "cmDBH_aboveN50.jpeg")), height = 12, width = 15, units='in')

  ## selecting sp with avg abundance to be greater than a threshold, say 50
  n.threshold <- n.threshold
  # ggplot(mrate.long %>% subset(tlplevel == "comm" &
  #                              avg.abund >= n.threshold),
  #        aes(x = udi.best, y = mrate, color = avg.abund)) +
  #   scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
  #                        low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
  #   facet_grid(size ~ censusint.m, scales = "free_y") +
  #   geom_point() +
  #   geom_smooth(method = "loess", color = "black") +
  #   ylab(expression("Mean Mortality Rate (per year)")) + xlab("Water Uptake Depth Index (m)")
  # ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_Mortality_rate_by_period_udi_by_size_aboveN", n.threshold, "_tlpcommn.jpeg")), height = 7, width = 12, units='in')

  ## Restricting analysis to only canopy and large canopy species---------
  # Metadata Panama Traits Definitions.doc:
  # GRWFRM2$ - Values are as in GRWFRM1 except free-standing species can have multiple values if maximum size varies widely within Panama
  # GRWFRM1$ - Values are Climber, HERB, S, U, M and T. S, U, M and T are free-standing species with maximum heights of 5, 10, 20 and > 30 m, respectively.

  bci.traits <- read.csv("data-raw/traits/BCITRAITS_20101220.csv") %>%
    rename(form1 = GRWFRM1., form2 = GRWFRM2., sp = SP.) %>% mutate(sp = tolower(sp))

  # load(file = "results/sp.mode.canopy.understorey")
  canopy.bci.traits <- bci.traits %>% subset(form1 %in% c("T")) %>% droplevels()
  canopy.sp <- unique(canopy.bci.traits$sp)

  ggplot(mrate.long %>% subset(!is.na(size) & !is.na(udi.best) & sp %in% canopy.sp),
         aes(x = udi.best, y = mrate, color = avg.abund)) +
    scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
                         low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
    facet_grid(size ~ censusint.m, scales = "free_y") +
    geom_point() +
    geom_smooth(method = "loess", color = "black") +
    scale_x_continuous(trans="sqrt", breaks = c(soil.depths[1], signif(soil.depths[c(8, 10, 12)], 1))) +
    ylab(expression("Mean Mortality Rate (% per year)")) + xlab("Water Uptake Depth Index (m)") +
    ggtitle("Species-level mortality (all sizes included)") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(paste0(figures.folder,
                          "/sp_Mortality_rate_by_period_udi_by_size_aboveN",
                          n.threshold, ".jpeg")), height = 7, width = 13, units='in')
  ggplot(sp.mrate.long %>% subset(!is.na(udi.best) & sp %in% canopy.sp),
         aes(x = udi.best, y = mrate, color = udi.best)) +
    scale_color_gradient(name = "UDI", trans = "rev_sqrt",
                         low = "red", high = "blue") +
    # scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
    #                      low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
    facet_grid(. ~ censusint.m) +
    geom_point() +
    geom_smooth(method = "loess", color = "black") +
    scale_x_continuous(trans="sqrt", breaks = c(soil.depths[1], signif(soil.depths[c(8, 10, 12)], 1))) +
    ylab(expression("Mean Mortality Rate (% per year)")) + xlab("Water Uptake Depth Index (m)") +
    ggtitle("Species-level mortality (all sizes included)") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(paste0(figures.folder,
                          "/sp_sp.Mortality_rate_by_period_udi_aboveN",
                          n.threshold, ".jpeg")), height = 3, width = 13, units='in')

  ggplot(sp.mrate.long %>% subset(!is.na(udi.best) & avg.abund >= n.threshold & sp %in% canopy.sp),
         aes(y = udi.best, x = mrate, color = udi.best)) +
    scale_color_gradient(name = "UDI", trans = "rev_sqrt",
                         low = "red", high = "blue", breaks = c(0, 0.5, 2.5, 5, 7.5)) +
    facet_grid(censusint.m ~ .) +
    geom_point() +
    geom_smooth(method = "loess", color = "black") +
    scale_y_continuous(trans="rev_sqrt",
                       breaks = c(round(soil.depths[c(8, 10, 12)], 0))) +
    xlab(expression("Mean Mortality Rate (% per year)")) + ylab("Water Uptake Depth Index (m)") +
    ggtitle("Species-level mortality (all sizes included)") +
    theme(plot.title = element_text(hjust = 0.3))
  ggsave(file.path(paste0(figures.folder,
                          "/sp_sp.Mortality_rate_by_period_udi_aboveN",
                          n.threshold, "_vertical.jpeg")), height = 13, width = 5, units='in')
  formula = y ~ x
  p0 <- ggplot(sp.mrate.long %>% subset(!is.na(udi.best) & sp %in% canopy.sp),
               aes(y = udi.best, x = mrate, color = udi.best)) +
    scale_color_gradient(name = "UDI", trans = "rev_sqrt",
                         low = "red", high = "blue", breaks = c(0.5, 2.5, 5, 7.5, 10)) +
    facet_grid(censusint.m ~ .) +
    geom_point() +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_poly_eq(aes(label = paste(..rr.label..)),
                 npcx = 0.85, npcy = 0.15, rr.digits = 2,
                 formula = formula, parse = TRUE, size = 3) +
    stat_fit_glance(method = 'lm',
                    method.args = list(formula = formula),
                    geom = 'text_npc',
                    aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                    npcx = 0.85, npcy = 0.05, size = 3) +
    scale_y_continuous(trans="rev_sqrt", breaks = c(soil.depths[1], signif(soil.depths[c(8, 10, 12)], 1))) +
    xlab(expression("Mean Mortality Rate (% per year)")) + ylab("Water Uptake Depth Index (m)") +
    ggtitle("Species-level mortality (all sizes included)") +
    theme(plot.title = element_text(hjust = 0.3))
  ggsave(file.path(paste0(figures.folder, "/sp.Mortality_rate_by_period_udi_aboveN",
                          n.threshold, "_vertical_lm.jpeg")),
         plot = p0, height = 9, width = 4, units='in')
  p1 <- p0 +
    facet_grid(censusint.m ~ Deciduousness)
  ggsave(file.path(paste0(figures.folder, "/sp.Mortality_rate_by_period_udi_aboveN",
                          n.threshold, "_vertical_lm_Deci.jpeg")),
         plot = p1, height = 9, width = 10, units='in')
  sp.mrate.long.class <- sp.mrate.long %>% subset(avg.abund >= n.threshold & sp %in% canopy.sp) %>%
    mutate(udi.best.class = forcats::fct_explicit_na(cut(udi.best, breaks = c(0, 2.5, 5, 11)))) %>%
    mutate(udi.best.class = factor(udi.best.class, levels = c("(Missing)", "(5,11]" , "(2.5,5]", "(0,2.5]"))) %>%
    group_by(udi.best.class, censusint.m) %>%
    summarise(mean.mrate = mean(mrate, na.rm = TRUE),
              se = sd(mrate, na.rm = TRUE)/sqrt(n()))

  ggplot(sp.mrate.long.class %>% subset(!is.na(udi.best.class)),
         aes(y = udi.best.class, x = mean.mrate)) +
    # scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
    #                      low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
    # scale_color_gradient(name = "UDI", trans = "rev_sqrt",
    #                      low = "red", high = "blue") +
    facet_grid(censusint.m ~ .) +
    geom_point() +
    geom_errorbarh(aes(xmax = mean.mrate + se, xmin = mean.mrate - se), size = 0.5) +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    xlab(expression("Mean Mortality Rate (% per year)")) + ylab("Water Uptake Depth Index (m)") +
    ggtitle("Species-level mortality (all sizes included)") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(paste0(figures.folder,
                          "/sp.Mortality_rate_by_period_udi_aboveN",
                          n.threshold, "_UDI.class.jpeg")), height = 13, width = 5, units='in')
  ###---------large trees' udi vs. sp-level tree mortality----------
  ## plotting diff drought-induced - background:
  sp.mrate.long.bck <- sp.mrate.long %>%
    subset(censusint.m %in% c("1982-85") & sp %in% canopy.sp & !is.na(udi.best)) %>% # c("1990-95", "1995-00")
    mutate(sp = as.character(sp)) %>%
    group_by(avg.abund, sp, udi.best, Deciduousness) %>%
    summarise(mrate = mean(mrate, na.rm = TRUE)) %>%
    mutate(period = "1982-85")
  sp.mrate.long.drought <- sp.mrate.long %>%
    subset(censusint.m %in% c("2000-05") & sp %in% canopy.sp & !is.na(udi.best)) %>% # c("2000-05")
    group_by(avg.abund, sp, udi.best, Deciduousness) %>%
    summarise(mrate = mean(mrate, na.rm = TRUE)) %>%
    mutate(period = "2000-2005")
  sp.mrate.long.diff <- bind_rows(sp.mrate.long.drought, sp.mrate.long.bck) %>% arrange(sp, period) %>%
    group_by(sp) %>%
    mutate(mrate.diff = mrate - lag(mrate)) %>% ungroup(sp)
  sp.mrate.long.diff
  formula <- y ~ x
  g0 <- ggplot(sp.mrate.long.diff, #  %>% subset(udi.best < 0.8)
               aes(y = udi.best, x = mrate.diff, color = udi.best)) +
    scale_color_gradient(name = "UDI", trans = "rev_sqrt",
                         low = "red", high = "blue", breaks = c(0, 0.5, 2.5, 5, 7.5)) +
    geom_point() +
    geom_vline(xintercept = 0) +
    geom_smooth(method = "lm", color = "black") +
    scale_y_continuous(trans="rev_sqrt",
                       breaks = c(0.01, signif(soil.depths[c(8, 10, 12)], 1))) +
    stat_poly_eq(aes(label = paste(..rr.label..)),
                 npcx = 0.9, npcy = 0.2, rr.digits = 2,
                 formula = formula, parse = TRUE, size = 4) +
    stat_fit_glance(method = 'lm',
                    method.args = list(formula = formula),
                    geom = 'text_npc',
                    aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                    npcx = 0.9, npcy = 0.1, size = 4) +
    xlab("2000-2005 - 1982-1985\nDifference Mean Mortality Rate (% per year)") +
    ylab("Water Uptake Depth Index (m)") +
    theme(axis.title.x = element_text(vjust = 0, hjust = 0.5)) +
    ggtitle("Species-level mortality (all sizes included)")
  g0.1 <- g0 + theme(plot.margin = margin(2, 2, 2, 2, "cm"))
  ggsave(file.path(paste0(figures.folder,
                          "/sp.Mortality_rate_by_period_udi_aboveN",
                          n.threshold, "_diff_mean_2000-2005 - 1982-1985.jpeg")),
         plot = g0.1, height = 5, width = 5.5, units='in')
  g1 <- g0 %+% subset(sp.mrate.long.diff, Deciduousness != "D") +
    facet_wrap(. ~ Deciduousness)
  ggsave(file.path(paste0(figures.folder,
                          "/sp.Mortality_rate_by_period_udi_aboveN",
                          n.threshold, "_diff_mean_2000-2005 - 1982-1985_Deci.jpeg")),
         plot = g1, height = 5, width = 7, units='in')

  #####--------large trees' udi vs. adult tree mortality------------
  ## plotting diff drought-induced - background:
  adult.mrate.long.bck <- adult.mrate.long %>%
    subset(censusint.m %in% c("1982-85") & sp %in% canopy.sp & !is.na(udi.best)) %>% # c("1990-95", "1995-00")
    mutate(sp = as.character(sp)) %>%
    group_by(avg.abund, sp, udi.best, Deciduousness) %>%
    summarise(mrate = mean(mrate, na.rm = TRUE)) %>%
    mutate(period = "1982-85")
  adult.mrate.long.drought <- adult.mrate.long %>%
    subset(censusint.m %in% c("2000-05") & sp %in% canopy.sp & !is.na(udi.best)) %>% # c("2000-05")
    group_by(avg.abund, sp, udi.best, Deciduousness) %>%
    summarise(mrate = mean(mrate, na.rm = TRUE)) %>%
    mutate(period = "2000-2005")
  adult.mrate.long.diff <- bind_rows(adult.mrate.long.drought, adult.mrate.long.bck) %>% arrange(sp, period) %>%
    group_by(sp) %>%
    mutate(mrate.diff = mrate - lag(mrate)) %>% ungroup(sp)
  adult.mrate.long.diff
  formula <- y ~ x
  g0 <- ggplot(adult.mrate.long.diff, #  %>% subset(udi.best < 0.8)
               aes(y = udi.best, x = mrate.diff, color = udi.best)) +
    scale_color_gradient(name = "UDI", trans = "rev_sqrt",
                         low = "red", high = "blue", breaks = c(0, 0.5, 2.5, 5, 7.5)) +
    geom_point() +
    geom_vline(xintercept = 0) +
    geom_smooth(method = "lm", color = "black") +
    scale_y_continuous(trans="rev_sqrt",
                       breaks = c(0.01, signif(soil.depths[c(8, 10, 12)], 1))) +
    stat_poly_eq(aes(label = paste(..rr.label..)),
                 npcx = 0.9, npcy = 0.2, rr.digits = 2,
                 formula = formula, parse = TRUE, size = 4) +
    stat_fit_glance(method = 'lm',
                    method.args = list(formula = formula),
                    geom = 'text_npc',
                    aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                    npcx = 0.9, npcy = 0.1, size = 4) +
    xlab("2000-2005 - 1982-1985\nDifference Mean Mortality Rate (% per year)") +
    ylab("Water Uptake Depth Index (m)") +
    ggtitle("Adult tree mortality (>= 10 cm DBH)") +
    theme(axis.title.x = element_text(vjust = 0, hjust = 0.5), plot.title = element_text(hjust = 0.5))
  g0.1 <- g0 + theme(plot.margin = margin(2, 2, 2, 2, "cm")) + theme(plot.title = element_text(vjust = 3))
  ggsave(file.path(paste0(figures.folder,
                          "/adult.Mortality_rate_by_period_udi_aboveN",
                          n.threshold, "_diff_mean_2000-2005 - 1982-1985.jpeg")),
         plot = g0.1, height = 5, width = 5.5, units='in')
  g1 <- g0 %+% subset(adult.mrate.long.diff, Deciduousness != "D") +
    facet_wrap(. ~ Deciduousness)
  ggsave(file.path(paste0(figures.folder,
                          "/adult.Mortality_rate_by_period_udi_aboveN",
                          n.threshold, "_diff_mean_2000-2005 - 1982-1985_Deci.jpeg")),
         plot = g1, height = 5, width = 7, units='in')
  ## by censuses
  formula = y ~ x
  p0 <- ggplot(adult.mrate.long %>% subset(!is.na(udi.best) & sp %in% canopy.sp),
               aes(y = udi.best, x = mrate, color = udi.best)) +
    scale_color_gradient(name = "UDI", trans = "rev_sqrt",
                         low = "red", high = "blue", breaks = c(0.5, 2.5, 5, 7.5, 10)) +
    facet_grid(censusint.m ~ .) +
    geom_point() +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_poly_eq(aes(label = paste(..rr.label..)),
                 npcx = 0.85, npcy = 0.15, rr.digits = 2,
                 formula = formula, parse = TRUE, size = 3) +
    stat_fit_glance(method = 'lm',
                    method.args = list(formula = formula),
                    geom = 'text_npc',
                    aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                    npcx = 0.85, npcy = 0.05, size = 3) +
    scale_y_continuous(trans="rev_sqrt", breaks = c(soil.depths[1], signif(soil.depths[c(8, 10, 12)], 1))) +
    xlab(expression("Mean Mortality Rate (% per year)")) + ylab("Water Uptake Depth Index (m)") +
    ggtitle("Adult tree mortality (>= 10 cm DBH)") +
    theme(plot.title = element_text(hjust = 0.3))
  ggsave(file.path(paste0(figures.folder, "/adult.Mortality_rate_by_period_udi_aboveN",
                          n.threshold, "_vertical_lm.jpeg")),
         plot = p0, height = 9, width = 4, units='in')
  p1 <- p0 +
    facet_grid(censusint.m ~ Deciduousness)
  ggsave(file.path(paste0(figures.folder, "/adult.Mortality_rate_by_period_udi_aboveN",
                          n.threshold, "_vertical_lm_Deci.jpeg")),
         plot = p1, height = 9, width = 10, units='in')

  adult.mrate.mean <- adult.mrate.long %>%
    group_by(sp, Deciduousness) %>%
    summarize_at(vars(mrate, avg.abund, udi.best), mean, na.rm = TRUE) %>%
    mutate(mrate = ifelse(!is.finite(mrate),
                          rep(NA, length(mrate)), mrate))
  pm.2 <- ggplot(adult.mrate.mean %>%
                     subset(avg.abund >= n.threshold & Deciduousness %in% c("DF", "Evg")),
                 aes(y = udi.best, x = mrate)) +
    geom_point() +
    facet_wrap(. ~ size, scales = "free_y") +
    scale_x_continuous(trans="sqrt", breaks = soil.depths[-c(2,3, 4, 6, 7, 9)]) +
    stat_poly_eq(aes(label = paste(..rr.label..)),
                 npcx = 0.1, npcy = 0.15, rr.digits = 2,
                 formula = formula, parse = TRUE, size = 4) +
    stat_fit_glance(method = 'lm',
                    method.args = list(formula = formula),
                    geom = 'text_npc',
                    aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                    npcx = 0.1, npcy = 0.05, size = 4) +
    scale_x_log10() + scale_y_reverse() +
    xlab(expression("Mean Mortality Rate (% per year)")) + ylab("Water Uptake Depth Index (m)") +
    facet_grid(. ~ Deciduousness)  +  geom_smooth(method = "lm", se = FALSE) +
    theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
    ggtitle("Adult tree mortality (>= 10 cm DBH)") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(paste0(figures.folder, "/adult_Mortality_vs_udi_with_outliers_avg.abund_above",
                          n.threshold, "_sp_size_deci.jpeg")), plot = pm.2, height = 4, width = 6, units='in')


  ####--
  ##### Is this because UDI has a relationship with mean growth rate?----
  demo.sp_size <- demo.sp_size  %>% left_join(deci, by = "sp") %>%
    mutate(Deciduousness = replace_na(Deciduousness, "Evg")) %>%
    mutate(DeciLvl = as.numeric(recode_factor(Deciduousness, !!!deci.level_key))) %>%
    droplevels()
  demo.sp_size.udi <- left_join(demo.sp_size, select(udi, sp_size, udi.best), by = "sp_size") %>%
    transform(sp_size = reorder(sp_size, -udi.best)) %>%
    separate(sp_size, c("sp", "size", sep = "_"), remove = FALSE, extra = "drop", fill = "right") %>%
    select(-"_") %>% subset(!size == "NA") %>% droplevels() %>%
    mutate(size = factor(size, levels = c("tiny", "small", "medium", "large"))) %>%
    subset(sp %in% canopy.sp) ## only canopy species

  pg.0 <- ggplot(demo.sp_size.udi, aes(y = udi.best, x = grate)) +
    geom_point() +
    facet_wrap(. ~ size, scales = "free_y") +
    stat_poly_eq(aes(label = paste(..rr.label..)),
                 npcx = 0.9, npcy = 0.19, rr.digits = 2,
                 formula = formula, parse = TRUE, size = 4) +
    stat_fit_glance(method = 'lm',
                    method.args = list(formula = formula),
                    geom = 'text_npc',
                    aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                    npcx = 0.9, npcy = 0.09, size = 4) +
    scale_x_log10() + scale_y_reverse() +
    xlab(expression("Mean Growth Rate (mm/yr)")) +  ylab("Water Uptake Depth Index (m)")
  pg.01 <- pg.0 + geom_smooth(method = "lm")
  ggsave(file.path(paste0(figures.folder, "/sp_mean_Growth_vs_udi_size.jpeg")), plot = pg.01, height = 6, width = 6, units='in')

  pg.1 <- pg.0 %+%  subset(demo.sp_size.udi, avg.abund >= n.threshold) +
    facet_grid(size ~ Deciduousness) +   geom_smooth(method = "lm", se = FALSE)
  ggsave(file.path(paste0(figures.folder, "/sp_mean_Growth_vs_udi_size_deci.jpeg")), plot = pg.1, height = 9, width = 10, units='in')

  pg.2 <- pg.0 %+%  subset(demo.sp_size.udi, avg.abund >= n.threshold & size == "large" &
                             Deciduousness %in% c("DF", "Evg")) +
    facet_grid(. ~ Deciduousness)  +  geom_smooth(method = "lm", se = FALSE) +
    theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
    ggtitle("Large Size class (>= 30 cm DBH)") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(paste0(figures.folder, "/sp_mean_Growth_vs_udi_size_deci_large.jpeg")), plot = pg.2, height = 4, width = 6, units='in')

  pm.0 <- ggplot(demo.sp_size.udi %>% subset(avg.abund >= n.threshold), aes(y = udi.best, x = mrate)) +
    geom_point() +
    facet_wrap(. ~ size, scales = "free_y") +
    scale_x_continuous(trans="sqrt", breaks = soil.depths[-c(2,3, 4, 6, 7, 9)]) +
    stat_poly_eq(aes(label = paste(..rr.label..)),
                 npcx = 0.1, npcy = 0.15, rr.digits = 2,
                 formula = formula, parse = TRUE, size = 4) +
    stat_fit_glance(method = 'lm',
                    method.args = list(formula = formula),
                    geom = 'text_npc',
                    aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                    npcx = 0.1, npcy = 0.05, size = 4) +
    scale_x_log10() + scale_y_reverse() +
    xlab(expression("Mean Mortality Rate (% per year)")) + ylab("Water Uptake Depth Index (m)")
  pm.01 <-   pm.0 + geom_smooth(method = "lm")
  ggsave(file.path(paste0(figures.folder, "/sp_mean_Mortality_vs_udi_with_outliers_avg.abund_above", n.threshold, "_sp_size.jpeg")),
         plot = pm.01, height = 5, width = 6.5, units='in')

  pm.1 <- pm.0 %+%  subset(demo.sp_size.udi, avg.abund >= n.threshold) +
    facet_grid(size ~ Deciduousness) +   geom_smooth(method = "lm", se = FALSE)
  ggsave(file.path(paste0(figures.folder, "/sp_mean_Mortality_vs_udi_with_outliers_avg.abund_above", n.threshold, "_sp_size_deci.jpeg")),
         plot = pm.1, height = 9, width = 10, units='in')

  pm.2 <- pm.0 %+%  subset(demo.sp_size.udi, avg.abund >= n.threshold & size == "large" &
                             Deciduousness %in% c("DF", "Evg")) +
    facet_grid(. ~ Deciduousness)  +  geom_smooth(method = "lm", se = FALSE) +
    theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
    ggtitle("Large Size class (>= 30 cm DBH)") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(paste0(figures.folder, "/sp_mean_Mortality_vs_udi_with_outliers_avg.abund_above",
                          n.threshold, "_sp_size_deci_large.jpeg")), plot = pm.2, height = 4, width = 6, units='in')
  ## Mortality vs growth rate
  p.d0 <- ggplot(demo.sp_size %>% subset(size != "NA"), aes(y = mrate, x = grate)) +
    geom_point() +
    facet_wrap(. ~ size, scales = "free_y") +
    geom_smooth(method = "lm") +
    scale_x_log10() + scale_y_log10() +
    xlab(expression("Mean Growth Rate (mm/yr)")) +  ylab(expression("Mean Mortality Rate (% per year)"))
  p.d0.1 <- p.d0 +  stat_poly_eq(npcx = 0.1, npcy = 0.2, size = 4, aes(label = paste(..rr.label..)), rr.digits = 2, formula = formula, parse = TRUE) +
    stat_fit_glance(npcx = 0.1, npcy = 0.1, size = 4, method = 'lm',  method.args = list(formula = formula), geom = 'text_npc', aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")))
  ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_mean_Growth_vs_mrate.jpeg")),
         plot = p.d0.1, height = 6, width = 6, units='in')

  p.d1 <-  p.d0 + facet_wrap(size ~ Deciduousness, scales = "free_y") + theme(axis.text.x = element_text(angle = 90)) +
    stat_poly_eq(npcx = 0.05, npcy = 0.2, size = 4, aes(label = paste(..rr.label..)), rr.digits = 2, formula = formula, parse = TRUE) +
    stat_fit_glance( npcx = 0.05, npcy = 0.1, size = 4, method = 'lm',  method.args = list(formula = formula), geom = 'text_npc', aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")))
  ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_mean_Growth_vs_mrate_Deci.jpeg")),
         plot = p.d1, height = 9, width = 10, units='in')
  p.d2 <-  p.d0 %+% subset(demo.sp_size, avg.abund >= n.threshold  &
                             !Deciduousness %in% c("DO")) +
    facet_wrap(. ~ Deciduousness, scales = "free_y") + theme(axis.text.x = element_text(angle = 90)) +
    stat_poly_eq(npcx = 0.85, npcy = 0.2, size = 4, aes(label = paste(..rr.label..)), rr.digits = 2, formula = formula, parse = TRUE) +
    stat_fit_glance( npcx = 0.85, npcy = 0.1, size = 4, method = 'lm',  method.args = list(formula = formula), geom = 'text_npc', aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")))+
    ggtitle("Large Size class (>= 30 cm DBH)") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_mean_Growth_vs_mrate_Deci_avg.abund_above", n.threshold,"_large.jpeg")),
         plot = p.d2, height = 6, width = 6, units='in')
}
