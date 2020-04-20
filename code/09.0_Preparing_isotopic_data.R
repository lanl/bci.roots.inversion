
#-----------------------------------------------------
# Title: Preparing isotopic data for inverse model evaluation
# Author : Rutuja Chitra-Tarak
# Original date: December 18, 2019
#-----------------------------------------------------

rm(list=ls())
gc()

iso.compare <- function(goodness.fit = goodness.fit,
                         dryseason = dryseason,
                         root.selection =  root.selection,
                         iso.subset = iso.subset, drop.months = drop.months){
  pacman::p_load(tidyverse, scales, ggpmisc)
  # graphics info
  theme_set(theme_bw())
  theme_update(text = element_text(size = 14),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.background = element_blank()
  )

  rev_sqrt_trans <- function() {
    scales::trans_new(
      name = "rev_sqrt",
      transform = function(x) -sqrt(abs(x)),
      inverse = function(x) x^2);
  }

  # load interval and working.iter
  load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
  load(file.path("results/GLUEsetup_part2_BCI.RData")) # has working.iter and growth and si matrix

  soil.depths <- unique(info$root.param.long$depth)
  intervals <- info$intervals
  n.ensembles <- growth_by_si.info$n.ensembles
  growth.type <- growth_by_si.info$growth.type
  growth.selection <- growth_by_si.info$growth.selection
  dbh.residuals <- growth_by_si.info$dbh.residuals
  si.type <- growth_by_si.info$si.type
  ##

  file.extension.base4 <- paste0("drop.months", drop.months, "_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type,
                                 "_", growth.selection, "_", dbh.residuals, "_", intervals,
                                 "_dryseason_", dryseason, "_iso.subset_", iso.subset, "_root.selection_", root.selection)

  load(file = paste0("results/splevel/ds.bestfit_", file.extension.base4, ".Rdata"))
  ds <- ds.bestfit
  load(file = paste0("results/commlevel/ds.bestfit_", file.extension.base4, ".Rdata"))
  ds <- ds.bestfit
  ds <- rbind(ds, ds.bestfit)

  ds <- ds %>% mutate(tlplevel = as.factor(tlplevel)) %>%
    subset(!is.na(udi.med.rsq)) %>% droplevels()

  iso <- read.csv("data-raw/traits/isotopes/Meinzer1999_Table1_Xylem_Sap_deltaD_Fig4.csv", na.strings = c("NA", ""), header = T, row.names = NULL, check.names = F)

  head(iso)
  # adding species codes
  load(file.path("data-raw/CTFScensuses/bci.spptable.rdata"))
  write.csv(bci.spptable, file.path("data-raw/CTFScensuses/bci.spptable.csv"), row.names = FALSE)

  iso$genus.sp <- paste(iso$Genus, iso$Species, sep = " ")
  bci.spptable$genus.sp <- paste(bci.spptable$Genus, bci.spptable$Species, sep = " ")
  matchrows = match(iso$genus.sp, bci.spptable$genus.sp)
  iso$genus.sp

  iso$sp <- bci.spptable$sp[matchrows]

  iso$sp
  ## some codes are not found, because of change of latin names, so searching in synonymous name-field should return the match
  match("Lonchocarpus", bci.spptable$Genus) #"Lonchocarpus latifolius"

  grep("Lonchocarpus", bci.spptable$syn)
  iso$sp[which(iso$Genus == "Lonchocarpus")] <- bci.spptable$sp[grep("Lonchocarpus", bci.spptable$syn)]
  head(iso)

  #"Guapira" ; iso$Family[11]
  grep("Guapira", bci.spptable$syn)
  iso$sp[which(iso$Genus == "Guapira")] <- bci.spptable$sp[grep("Guapira", bci.spptable$syn)]# "Guapira" "campochagres"
  iso$sp

  ## the following species code need to be matched with that in ds.bestfit
  # isotopic uptake depth is available for Quararibea asterolepis subsp. stenophylla, but
  # modelled uptake depth is for Quararibea asterolepis, so making them synonymous for comparison
  iso$sp[which(iso$sp == "quara1")] <- "quaras"
  ## adding phonology score

  iso$PhenoRank <- recode(iso$Phenology, `Evergreen` = 1, `Facultatively deciduous` = 2,
                          Brevideciduous = 3,  `Partially deciduous` = 4,  `Deciduous` = 5)

  iso$Phenology <- factor(iso$Phenology, levels = c("Evergreen", "Facultatively deciduous",
                                                    "Brevideciduous", "Partially deciduous", "Deciduous"))

  iso.sp <- iso$sp #[!is.na(iso$Xylem_sap_deltaD_permil)]
  save(iso.sp, file = "results/sp_with_isotopic_record.Rdata")
  load("results/sp_with_isotopic_record.Rdata")
  write.csv(iso, file.path(paste0("data-raw/traits/isotopes/Meinzer1999_Table1_Xylem_Sap_deltaD_Fig4_sp_code.csv")), row.names = FALSE)

  ### More data from March1997 from Fig 5B (includes data in Fig 5A)
  iso.2.raw <- read.csv("data-raw/traits/isotopes/Meinzer1999_Xylem_Sap_deltaD_March97_DBH_Fig5B.csv", na.strings = c("NA", ""), header = T, row.names = NULL, check.names = F)
  iso.change <- read.csv("data-raw/traits/isotopes/Meinzer1999_Xylem_Sap_deltaD_by_delta_SapFlow_&_monthly_LeafFall_Fig7AB.csv", na.strings = c("NA", ""), header = T, row.names = NULL, check.names = F)

  head(iso.2.raw)
  iso.2 <- iso.2.raw %>% group_by(sp) %>%
    summarise(Xylem_sap_deltaD_permil = mean(Xylem_sap_deltaD_permil, na.rm = TRUE),
              n = n(),
              SD = sd(Xylem_sap_deltaD_permil, na.rm = TRUE),
              SE = SD/sqrt(n),
              DBH = mean(DBH, na.rm = TRUE))
  iso.2 <- iso.2 %>% left_join(iso %>% select(sp, Phenology), by = "sp") %>%
    left_join(iso.change, by = "sp")

  #### Plotting soil isotopic variation from Meinzer et al and Jackson et al -----------
  iso.soil.1 <- read.csv("data-raw/traits/isotopes/Meinzer1999_Fig2A_soil_deltaD_BCI_data_Mar_Apr_1997.csv", na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)
  iso.soil.2 <- read.csv("data-raw/traits/isotopes/Oecologia 1995 Jackson_fig2_soil_deltaD_Gigante_data_dry_season_1992.csv", na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)
  depth.m1 <- lm(depth ~ soil.deltaD, data = iso.soil.1)
  # summ.depth.m1 <- summary(depth.m1)
  # soil.label <- expression('Soil '*delta~""^2*"H (\u2030)"*'')
  # ## expression(paste(delta^{2}, "H (\u2030)")))
  # depth.m1.label = paste0("Depth = ", round(depth.m1$coefficients[1], 2), " + ",
  #                         round(depth.m1$coefficients[2], 2),
  #                         " * Soil delta 2H\nR-squared = ", round(summ.depth.m1$r.squared, 2),
  #                         ", p-val = ", round(summ.depth.m1$coefficients[2, 4], 4))
  #
  # g1 <- ggplot(iso.soil.1, aes(y = depth, x = soil.deltaD)) +
  #   geom_point(size = 3) +
  #   xlab(soil.label) + ylab("Depth (cm)") +
  #   scale_y_reverse()
  # g1 +   geom_smooth(method = "lm", se = FALSE) +
  #   geom_errorbarh(aes(xmax = soil.deltaD + se, xmin = soil.deltaD - se), size = 0.5) +
  #   geom_text(aes(x = -40, y = 0, label = depth.m1.label), color = "blue")
  # ggsave(file.path(paste0("figures/UDI_confidence/Meinzer_etal_1999_Depth_vs_soil_deltaD.jpeg")), height = 5, width = 6, units='in')
  #
  # g1 +   geom_smooth(method = "loess", se = FALSE, span = 0.7) +
  #   geom_errorbarh(aes(xmax = soil.deltaD + se, xmin = soil.deltaD - se), size = 0.5)
  #
  # depth.m2 <- lm(depth ~ soil.deltaD, data = iso.soil.2)
  # summ.depth.m2 <- summary(depth.m2)
  # depth.m2.label = paste0("Depth = ", round(depth.m2$coefficients[1], 2), " + ",
  #                         round(depth.m2$coefficients[2], 2),
  #                         " * Soil delta 2H\nR-squared = ", round(summ.depth.m2$r.squared, 2),
  #                         ", p-val = ", round(summ.depth.m2$coefficients[2, 4], 4))
  # g1 %+% iso.soil.2 + geom_smooth(method = "glm", se = FALSE) +
  #   geom_text(aes(x = -35, y = 5, label = depth.m2.label), color = "blue")
  # ggsave(file.path(paste0("figures/UDI_confidence/Jackson_etal_1995_Depth_vs_soil_deltaD.jpeg")), height = 5, width = 6, units='in')


  ## adding Xylem sap's soil uptake depth ----
  # how to incorporate se in this?
  iso <- iso %>% mutate(depth = predict.lm(depth.m1, newdata = data.frame(soil.deltaD = Xylem_sap_deltaD_permil))/100,
                        depth.se = depth - predict.lm(depth.m1,
                                                      newdata = data.frame(soil.deltaD = Xylem_sap_deltaD_permil + SE))/100) # from cm to m
  head(iso)
  #### Joining isotopic and udi data-----------
  udi.all <- ds
  udi <- ds %>%
    subset(!is.na(udi.med.rsq) & size %in% c("large") &
             sp %in% iso$sp) %>% droplevels() #%>%
  # group_by(sp, tlplevel) %>%
  # summarise_at(vars(udi, udi.sd, sdi, sdi.sd, rsq.mean, rsq.min), mean, na.rm = TRUE)
  # udi <- udi.all %>% group_by(sp) %>% summarise(udi = mean(udi, na.rm = TRUE))
  #
  # udi.long <- udi.all %>% gather(key = key, value = R, -uptake.depth, -rsq, -sp_size)
  #
  # ggplot(udi.long, aes(x =  uptake.depth, y = R)) +
  #   geom_bar(aes(fill = key), stat  = "identity")  + theme_classic()
  # ggsave(file.path(paste0("figures/UDI_confidence/R_by_uptake_depth_all_sp_size", "_", model.type, "_", demand.var, ".tiff")), height = 5, width = 6, units ='in', compression = "lzw")
  #
  # ggsave(file.path(paste0("figures/UDI_confidence/R_by_uptake_depth_large_trees", "_", model.type, "_", demand.var, ".tiff")), height = 5, width = 6, units ='in', compression = "lzw")

  iso.udi <- left_join(iso, udi, by = "sp")
  iso.udi <- iso.udi[order(iso.udi$udi.med.rsq),]
  iso.udi <- iso.udi %>% subset(!is.na(iso.udi$udi.med.rsq) & !is.na(iso.udi$Xylem_sap_deltaD_permil))
  # View(iso.udi)
  iso.udi.sp <- iso.udi$sp
  save(iso.udi.sp, file = "results/sp_with_isotopic_record_&_udi.Rdata")

  tlp <- read.csv("data-raw/traits/HydraulicTraits_Kunert/tlp_sp_mean.csv", na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)
  iso.udi <- iso.udi %>% left_join(tlp %>% mutate(sp = as.character(sp)), by = "sp")
  #View(iso.udi)
  write.csv(iso.udi, file.path(paste0("results/iso.udi_", file.extension.base4, ".csv", sep = "")), row.names = FALSE)
