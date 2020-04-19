## Comparison with Meinzer et al. 1999 results.

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

  ## creating paths, if they don't exist-------
  if(!dir.exists(file.path("figures"))) {dir.create(file.path("figures"))}
  if(!dir.exists(file.path("figures", "UDI_confidence"))) {dir.create(file.path("figures","UDI_confidence"))}
  if(!dir.exists(file.path("figures", "UDI_confidence", growth.type))) {dir.create(file.path("figures","UDI_confidence", growth.type))}
  if(!dir.exists(file.path("figures", "UDI_confidence", growth.type, growth.selection))) {
    dir.create(file.path("figures","UDI_confidence", growth.type, growth.selection))}
  if(!dir.exists(file.path("figures", "UDI_confidence", growth.type, growth.selection,  paste0("dbh.residuals_", dbh.residuals)))) {
    dir.create(file.path("figures","UDI_confidence", growth.type, growth.selection,  paste0("dbh.residuals_", dbh.residuals)))}
  file.path.udi <- file.path("figures","UDI_confidence", growth.type, growth.selection,
                             paste0("dbh.residuals_", dbh.residuals), paste0("dryssn_", dryseason, "_root.selection_", root.selection, drop.months))
  if(!dir.exists(file.path.udi)) {dir.create(file.path.udi)}
  file.path.udi.rsq <- file.path(file.path.udi, "udi.med.rsq")
  file.path.udi.ll <- file.path(file.path.udi, "udi.med.ll")
  if(!dir.exists(file.path.udi.rsq)) {dir.create(file.path.udi.rsq)}
  if(!dir.exists(file.path.udi.ll)) {dir.create(file.path.udi.ll)}

  xylem.label <- expression('Xylem Sap '*delta~""^2*"H (\u2030)"*'')
  change.xylem.label <- expression('Change in Xylem Sap '*delta~""^2*"H (\u2030)"*'')

  iso.2.raw <- iso.2.raw %>% left_join(iso %>% select(sp, Phenology), by = "sp") %>%
    left_join(iso.change, by = "sp") %>% left_join(tlp, by = "sp")

  ggplot(iso.2.raw, aes(y =  Xylem_sap_deltaD_permil, x = DBH, color = tlp)) +
    geom_text(aes(y =  Xylem_sap_deltaD_permil + 1.5, x = DBH, label = sp),
              size = 2) +
    geom_point(aes(shape = Phenology), size = 3) +
    scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
    scale_y_continuous(limits = c(-70, -18), breaks = seq(-70, -20, by = 10)) +
    xlab("DBH (cm)") + ylab(xylem.label) +
    ggtitle("Species March 1997 Xylem Sap deltaD Vs DBH")
  ggsave(file.path(paste0("data-raw/traits/isotopes/Meinzer1999_Xylem_Sap_deltaD_March97_DBH_Fig5B.jpeg")), height = 5, width = 8, units ='in')

  ggplot(iso.2.raw, aes(y =  delta_sapflow_percent_per_day, x = delta_xylem_sap_deltaD_permil_per_day, color = tlp)) +
    geom_smooth(method = "lm", se = FALSE, size = 0.5, color = "black") +
    geom_text(aes(y =  delta_sapflow_percent_per_day + 0.05, label = sp),
              size = 4) +
    geom_point(aes(shape = Phenology), size = 3) +
    scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
    scale_x_continuous(limits = c(-0.6, 0.1)) +
    xlab(change.xylem.label) + ylab("Change in daily sap flow (% per day)") +
    ggtitle("Sap flow Vs. Xylem Sap deltaD\nWithin species variation in 1997 dry season")
  ggsave(file.path(paste0("data-raw/traits/isotopes/Meinzer1999_Xylem_Sap_deltaD_by_delta_SapFlow_&_monthly_LeafFall_Fig7A.jpeg")), height = 5, width = 7, units ='in')

  ggplot(iso.2.raw, aes(y = SE_mean_monthly_leaf_fall, x = delta_xylem_sap_deltaD_permil_per_day, color = tlp)) +
    geom_smooth(method = "lm", se = FALSE, size = 0.5, color = "black") +
    geom_text(aes(y = SE_mean_monthly_leaf_fall + 0.002, label = sp),
              size = 4) +
    geom_point(aes(shape = Phenology), size = 3) +
    scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
    scale_x_continuous(limits = c(-0.6, 0.1)) +
    xlab(change.xylem.label) + ylab("SE of mean monthly leaf fall") +
    ggtitle("Leaf fall Vs. Xylem Sap deltaD\nWithin species variation in 1997 dry season")
  ggsave(file.path(paste0("data-raw/traits/isotopes/Meinzer1999_Xylem_Sap_deltaD_by_delta_SapFlow_&_monthly_LeafFall_Fig7B.jpeg")), height = 5, width = 7, units ='in')


  # iso <- read.csv("data-raw/traits/isotopes/Oecologia 1995 Jackson _Fig3_Fig4.csv", na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)
  # iso.udi <- left_join(iso, udi, by = "sp")
  # iso.sp <- unique(iso.udi$sp[!is.na(iso.udiudi.med.rsq) & !is.na(iso.udi$Xylem_sap_deltaD_permil)])
  # iso.udi <- iso.udi[order(iso.udiudi.med.rsq),]

  tlplevels <- c("sp", "comm")
  subsetting <- c("with outliers", "without outliers")
  # udi.med.rsq -------

  for (i in 1: 2) {
    if(i == 2) {
      iso.udi.i <- iso.udi %>% subset(!sp %in% c("guapst"))
    } else {
      iso.udi.i <- iso.udi
    }
    for (j in 1: length(tlplevels)) {
      iso.udi.sub <- iso.udi.i %>% subset(tlplevel == tlplevels[j] & size == "large")
      formula <- y ~ x
      p0 <- ggplot(iso.udi.sub, aes(x =  Xylem_sap_deltaD_permil, y = udi.med.rsq, color = tlp)) +
        geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + SE, xmin = Xylem_sap_deltaD_permil - SE),
                       size = 0.5) +
        geom_errorbar(aes(ymax = udi.med.rsq + udi.sd.rsq, ymin = udi.med.rsq - udi.sd.rsq),
                      size = 0.5, width = 0.2) +
        scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
        geom_text(aes(x =  Xylem_sap_deltaD_permil + 1.5, y = udi.med.rsq + diff(range(iso.udi.sub$udi.med.rsq, na.rm = TRUE))/20, label = sp),
                  size = 4) +
        ylab(expression("Water Uptake Depth (m)")) + xlab(xylem.label) +
        # scale_y_continuous(trans="rev_sqrt", breaks = c(0.00001, soil.depths))
        scale_y_reverse() +
        stat_poly_eq(aes(label = paste(..rr.label..)),
                     npcx = 0.6, npcy = 0.1, rr.digits = 2,
                     formula = formula, parse = TRUE, size = 6) +
        stat_fit_glance(method = 'lm',
                        method.args = list(formula = formula),
                        geom = 'text_npc',
                        aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                        npcx = 0.85, npcy = 0.1, size = 6) +
        ggtitle(paste0("TLPlevel = ", tlplevels[j], ", ",subsetting[i], "\nSpecies Uptake Depth Vs Xylem Sap deltaD\n"))
      p0 +
        geom_point(size = 3, show.legend = TRUE)
      ggsave(file.path(file.path.udi.rsq, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
                                             goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 6, units ='in')
      p0 +
        geom_point(aes(shape = Phenology), size = 3)
      ggsave(file.path(file.path.udi.rsq, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_phenology_cor",
                                             goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 8, units ='in')
      ## only species with TLP data
      p0 %+% subset(iso.udi.sub, !is.na(tlp)) +
        geom_point(size = 3, show.legend = TRUE)
      ggsave(file.path(file.path.udi.rsq, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
                                             goodness.fit, "_", tlplevels[j], "_", subsetting[i], "_only_TLP.jpeg")), height = 5, width = 6, units ='in')

      ### Isotopic vs. Update depth
      p1 <- ggplot(iso.udi.sub, aes(x =  depth, y = udi.med.rsq, color = tlp)) +
        geom_errorbarh(aes(xmax = depth + depth.se, xmin = depth - depth.se),
                       size = 0.5) +
        geom_errorbar(aes(ymax = udi.med.rsq + udi.sd.rsq, ymin = udi.med.rsq - udi.sd.rsq),
                      size = 0.5, width = 0.01) +
        geom_point(size = 3, show.legend = TRUE) +
        scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
        geom_text(aes(x = depth, y = udi.med.rsq + diff(range(iso.udi.sub$udi.med.rsq, na.rm = TRUE))/20, label = sp),
                  size = 4) +
        stat_poly_eq(aes(label = paste(..rr.label..)),
                     npcx = 0.6, npcy = 0.1, rr.digits = 2,
                     formula = formula, parse = TRUE, size = 6) +
        stat_fit_glance(method = 'lm',
                        method.args = list(formula = formula),
                        geom = 'text_npc',
                        aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                        npcx = 0.85, npcy = 0.1, size = 6) +
        ggtitle(paste0("TLPlevel = ", tlplevels[j], ", ",subsetting[i],
                       "\nSpecies Uptake Depth: Modelled Vs. Isotopic\n")) +
        ylab(expression("Modelled Water Uptake Depth (m)")) + xlab("Isotopic Water Uptake Depth (m)") +
        scale_x_reverse() +   scale_y_reverse()
      #scale_y_continuous(trans="rev_sqrt", breaks = c(0.00001, soil.depths))
      p1
      ggsave(file.path(file.path.udi.rsq, paste0("Comparison_with_Meinzer1999_isotopic_vs.modelled_uptake.depth_cor",
                                             goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 6, units ='in')
    }
  }

  # udi.med.ll -------

  for (i in 1: 2) {
    if(i == 2) {
      iso.udi.i <- iso.udi %>% subset(!sp %in% c("guapst"))
    } else {
      iso.udi.i <- iso.udi
    }
    for (j in 1: length(tlplevels)) {
      iso.udi.sub <- iso.udi.i %>% subset(tlplevel == tlplevels[j] & size == "large")
      formula <- y ~ x
      p0 <- ggplot(iso.udi.sub, aes(x =  Xylem_sap_deltaD_permil, y = udi.med.ll, color = tlp)) +
        geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + SE, xmin = Xylem_sap_deltaD_permil - SE),
                       size = 0.5) +
        geom_errorbar(aes(ymax = udi.med.ll + udi.sd.ll, ymin = udi.med.ll - udi.sd.ll),
                      size = 0.5, width = 0.2) +
        scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
        geom_text(aes(x =  Xylem_sap_deltaD_permil + 1.5, y = udi.med.ll + diff(range(iso.udi.sub$udi.med.ll, na.rm = TRUE))/20, label = sp),
                  size = 4) +
        ylab(expression("Water Uptake Depth (m)")) + xlab(xylem.label) +
        # scale_y_continuous(trans="rev_sqrt", breaks = c(0.00001, soil.depths))
        scale_y_reverse() +
        stat_poly_eq(aes(label = paste(..rr.label..)),
                     npcx = 0.6, npcy = 0.1, rr.digits = 2,
                     formula = formula, parse = TRUE, size = 6) +
        stat_fit_glance(method = 'lm',
                        method.args = list(formula = formula),
                        geom = 'text_npc',
                        aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                        npcx = 0.85, npcy = 0.1, size = 6) +
        ggtitle(paste0("TLPlevel = ", tlplevels[j], ", ",subsetting[i], "\nSpecies Uptake Depth Vs Xylem Sap deltaD\n"))
      p0 +
        geom_point(size = 3, show.legend = TRUE)
      ggsave(file.path(file.path.sdi.ll, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
                                                 goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 6, units ='in')
      p0 +
        geom_point(aes(shape = Phenology), size = 3)
      ggsave(file.path(file.path.sdi.ll, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_phenology_cor",
                                                 goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 8, units ='in')
      ## only species with TLP data
      p0 %+% subset(iso.udi.sub, !is.na(tlp)) +
        geom_point(size = 3, show.legend = TRUE)
      ggsave(file.path(file.path.sdi.ll, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
                                                 goodness.fit, "_", tlplevels[j], "_", subsetting[i], "_only_TLP.jpeg")), height = 5, width = 6, units ='in')

      ### Isotopic vs. Update depth
      p1 <- ggplot(iso.udi.sub, aes(x =  depth, y = udi.med.ll, color = tlp)) +
        geom_errorbarh(aes(xmax = depth + depth.se, xmin = depth - depth.se),
                       size = 0.5) +
        geom_errorbar(aes(ymax = udi.med.ll + udi.sd.ll, ymin = udi.med.ll - udi.sd.ll),
                      size = 0.5, width = 0.01) +
        geom_point(size = 3, show.legend = TRUE) +
        scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
        geom_text(aes(x = depth, y = udi.med.ll + diff(range(iso.udi.sub$udi.med.ll, na.rm = TRUE))/20, label = sp),
                  size = 4) +
        stat_poly_eq(aes(label = paste(..rr.label..)),
                     npcx = 0.6, npcy = 0.1, rr.digits = 2,
                     formula = formula, parse = TRUE, size = 6) +
        stat_fit_glance(method = 'lm',
                        method.args = list(formula = formula),
                        geom = 'text_npc',
                        aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                        npcx = 0.85, npcy = 0.1, size = 6) +
        ggtitle(paste0("TLPlevel = ", tlplevels[j], ", ",subsetting[i],
                       "\nSpecies Uptake Depth: Modelled Vs. Isotopic\n")) +
        ylab(expression("Modelled Water Uptake Depth (m)")) + xlab("Isotopic Water Uptake Depth (m)") +
        scale_x_reverse() +   scale_y_reverse()
      #scale_y_continuous(trans="rev_sqrt", breaks = c(0.00001, soil.depths))
      p1
      ggsave(file.path(file.path.sdi.ll, paste0("Comparison_with_Meinzer1999_isotopic_vs.modelled_uptake.depth_cor",
                                                 goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 6, units ='in')
    }
  }
#
#   ### sdi.med.rsq -------
#   file.path.sdi.rsq <- file.path(file.path.udi, "sdi")
#   if(!dir.exists(file.path.sdi.rsq)) {dir.create(file.path.sdi.rsq)}
#   for (i in 1: 2) {
#     if(i == 2) {
#       iso.udi.i <- iso.udi %>% subset(!sp %in% c("guapst"))
#     } else {
#       iso.udi.i <- iso.udi
#     }
#     for (j in 1: length(tlplevels)) {
#       iso.udi.sub <- iso.udi.i %>% subset(tlplevel == tlplevels[j])
#       p0 <- ggplot(iso.udi.sub, aes(x =  Xylem_sap_deltaD_permil, y = sdi.med.rsq, color = tlp)) +
#         geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + SE, xmin = Xylem_sap_deltaD_permil - SE),
#                        size = 0.5) +
#         geom_errorbar(aes(ymax = sdi.med.rsq + sdi.sd.rsq, ymin = sdi.med.rsq - sdi.sd.rsq),
#                       size = 0.5, width = 0.2) +
#         scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
#         geom_text(aes(x =  Xylem_sap_deltaD_permil + 1.5, y = sdi.med.rsq + diff(range(iso.udi.sub$sdi.med.rsq, na.rm = TRUE))/20, label = sp),
#                   size = 4) +
#         stat_poly_eq(aes(label = paste(..rr.label..)),
#                      npcx = 0.6, npcy = 0.1, rr.digits = 2,
#                      formula = formula, parse = TRUE, size = 6) +
#         stat_fit_glance(method = 'lm',
#                         method.args = list(formula = formula),
#                         geom = 'text_npc',
#                         aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
#                         npcx = 0.85, npcy = 0.1, size = 6) +
#         ggtitle(paste0("TLPlevel = ", tlplevels[j], ", ",subsetting[i], "\nSpecies Uptake Depth Vs Xylem Sap deltaD\n")) +
#         ylab(expression("Water-Stress Depth (m)")) + xlab(xylem.label) + scale_y_reverse()
#       # scale_y_continuous(trans="rev_sqrt", breaks = c(0.00001, soil.depths))
#       p0 +
#         geom_point(size = 3, show.legend = TRUE)
#       #scale_color_continuous(name = "Rsq\nMean", trans = "reverse")
#       ggsave(file.path(file.path.sdi.rsq, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
#                                              goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 6, units ='in')
#     }
#   }
#

  ###### Best-----------
  file.path.udi.best.rsq <- file.path(file.path.udi, "udi.best.rsq")
  if(!dir.exists(file.path.udi.best.rsq)) {dir.create(file.path.udi.best.rsq)}
  file.path.udi.best.ll <- file.path(file.path.udi, "udi.best.ll")
  if(!dir.exists(file.path.udi.best.ll)) {dir.create(file.path.udi.best.ll)}

  ###### udi.best.rsq --------
  for (i in 1: 2) {
    if(i == 2) {
      iso.udi.i <- iso.udi %>% subset(!sp %in% c("guapst"))
    } else {
      iso.udi.i <- iso.udi
    }
    for (j in 1: length(tlplevels)) {
      iso.udi.sub <- iso.udi.i %>% subset(tlplevel == tlplevels[j])
      p0 <- ggplot(iso.udi.sub, aes(x =  Xylem_sap_deltaD_permil, y = udi.best.ll, color = tlp)) +
        geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + SE, xmin = Xylem_sap_deltaD_permil - SE),
                       size = 0.5) +
        # geom_errorbar(aes(ymax = udi.best.ll + udi.sd, ymin = udi.best.ll - udi.sd),
        #               size = 0.5, width = 0.2) +
        scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
        geom_text(aes(x =  Xylem_sap_deltaD_permil + 2.5, y = udi.best.ll + diff(range(iso.udi.sub$udi.best.ll, na.rm = TRUE))/20, label = sp),
                  size = 4) +
        ggtitle(paste0("TLPlevel = ", tlplevels[j], ", ",subsetting[i], "\nSpecies Uptake Depth Vs Xylem Sap deltaD\n")) +
        ylab(expression("Uptake Depth (m)")) + xlab(xylem.label) + scale_y_reverse()
      p0 +
        geom_point(size = 3, show.legend = TRUE) +
        stat_poly_eq(aes(label = paste(..rr.label..)),
                     npcx = 0.6, npcy = 0.1, rr.digits = 2,
                     formula = formula, parse = TRUE, size = 4) +
        stat_fit_glance(method = 'lm',
                        method.args = list(formula = formula),
                        geom = 'text_npc',
                        aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                        npcx = 0.87, npcy = 0.1, size = 4)
      ggsave(file.path(file.path.udi.best.ll, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
                                                  goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 6, units ='in')
      p0 +
        geom_point(aes(shape = Phenology), size = 3) +
        stat_poly_eq(aes(label = paste(..rr.label..)),
                     npcx = 0.6, npcy = 0.1, rr.digits = 2,
                     formula = formula, parse = TRUE, size = 4) +
        stat_fit_glance(method = 'lm',
                        method.args = list(formula = formula),
                        geom = 'text_npc',
                        aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                        npcx = 0.87, npcy = 0.1, size = 4)
      ggsave(file.path(file.path.udi.best.ll, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_phenology_cor",
                                                  goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 8, units ='in')
      if(i == 1) {
        p0 + geom_point(size = 3, show.legend = TRUE) +
          geom_smooth(method = "lm", se = FALSE, lty = "dotted", color = "darkgray", size = 0.5) +
          stat_poly_eq(aes(label = paste(..rr.label..)),
                       npcx = 0.6, npcy = 0.1, rr.digits = 2,
                       formula = formula, parse = TRUE, size = 4, color = "darkgray") +
          stat_fit_glance(method = 'lm',
                          method.args = list(formula = formula),
                          geom = 'text_npc',
                          aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                          npcx = 0.88, npcy = 0.1, size = 4, color = "darkgray") +
          geom_smooth(data = subset(iso.udi.sub, !sp %in% c("sponra", "guapst")),
                      method = "lm", se = FALSE, color = "black", size = 0.5) +
          stat_poly_eq(data = subset(iso.udi.sub, !sp %in% c("sponra", "guapst")), aes(label = paste(..rr.label..)),
                       npcx = 0.6, npcy = 0.2, rr.digits = 2,
                       formula = formula, parse = TRUE, size = 4) +
          stat_fit_glance(data = subset(iso.udi.sub, !sp %in% c("sponra", "guapst")), method = 'lm',
                          method.args = list(formula = formula),
                          geom = 'text_npc',
                          aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                          npcx = 0.9, npcy = 0.2, size = 4)
      }
      ggsave(file.path(file.path.udi.best.ll, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
                                                  goodness.fit, "_", tlplevels[j], "_", subsetting[i], "_lm.jpeg")), height = 5, width = 6, units ='in')
      ## without loncla
      if(j == 2) {
        p0 %+% subset(iso.udi.sub, !sp %in% c("loncla")) +
          geom_point(size = 3, show.legend = TRUE)
        ggsave(file.path(file.path.udi.best.ll, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
                                                    goodness.fit, "_", tlplevels[j], "_", subsetting[i], "_witout loncla.jpeg")), height = 5, width = 6, units ='in')
      }
      ## only species with TLP data
      p0 %+% subset(iso.udi.sub, !is.na(tlp)) +
        geom_point(size = 3, show.legend = TRUE)
      ggsave(file.path(file.path.udi.best.ll, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
                                             goodness.fit, "_", tlplevels[j], "_", subsetting[i], "_only_TLP.jpeg")), height = 5, width = 6, units ='in')
      p0 %+% subset(iso.udi.sub, !is.na(tlp)) +
        geom_point(aes(shape = Phenology), size = 3)
      ggsave(file.path(file.path.udi.best.ll, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_phenology_cor",
                                                  goodness.fit, "_", tlplevels[j], "_", subsetting[i], "_only_TLP.jpeg")), height = 5, width = 8, units ='in')

      p1 <- ggplot(iso.udi.sub, aes(x =  Xylem_sap_deltaD_permil, y = udi.best.ll)) +
        geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + SE, xmin = Xylem_sap_deltaD_permil - SE),
                       size = 0.5) +
        # geom_errorbar(aes(ymax = udi.best + udi.sd, ymin = udi.best - udi.sd),
        #               size = 0.5, width = 0.2) +
        geom_text(aes(x =  Xylem_sap_deltaD_permil + 2.2, y = udi.best.ll + diff(range(iso.udi.sub$udi.best.ll, na.rm = TRUE))/20, label = sp),
                  size = 4) +
        stat_poly_eq(aes(label = paste(..rr.label..)),
                     npcx = 0.6, npcy = 0.1, rr.digits = 2,
                     formula = formula, parse = TRUE, size = 6, color = "darkgray") +
        stat_fit_glance(method = 'lm',
                        method.args = list(formula = formula),
                        geom = 'text_npc',
                        aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                        npcx = 0.9, npcy = 0.1, size = 6, color = "darkgray") +
        ylab(expression("Water Uptake Depth (m)")) + xlab(xylem.label) + scale_y_reverse() +
        geom_point(size = 3, show.legend = TRUE) + theme(text = element_text(size = 22), axis.text = element_text(size = 22)) +
        geom_smooth(method = "lm", se = FALSE, lty = "dotted", color = "darkgray")
      if( i == 1) {
        p1 +
          geom_smooth(data = subset(iso.udi.sub, !sp %in% c("sponra", "guapst")),
                      method = "lm", se = FALSE, color = "black") +
          stat_poly_eq(data = subset(iso.udi.sub, !sp %in% c("sponra", "guapst")), aes(label = paste(..rr.label..)),
                       npcx = 0.6, npcy = 0.2, rr.digits = 2,
                       formula = formula, parse = TRUE, size = 6) +
          stat_fit_glance(data = subset(iso.udi.sub, !sp %in% c("sponra", "guapst")), method = 'lm',
                          method.args = list(formula = formula),
                          geom = 'text_npc',
                          aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                          npcx = 0.9, npcy = 0.2, size = 6)

      } else {
        p1
      }
      ggsave(file.path(file.path.udi.best.ll, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
                                                  goodness.fit, "_", tlplevels[j], "_", subsetting[i], "_no_color.jpeg")), height = 5, width = 5.25, units ='in')

      ### Isotopic vs. Update depth
      p1 <- ggplot(iso.udi.sub, aes(x =  depth, y = udi.best.ll, color = tlp)) +
        geom_errorbarh(aes(xmax = depth + depth.se, xmin = depth - depth.se),
                       size = 0.5) +
        # geom_errorbar(aes(ymax = udi.best.ll + udi.sd, ymin = udi.best.ll - udi.sd),
        #               size = 0.5, width = 0.01) +
        geom_point(size = 3, show.legend = TRUE) +
        scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
        geom_text(aes(x = depth, y = udi.best.ll + diff(range(iso.udi.sub$udi.best.ll, na.rm = TRUE))/20, label = sp),
                  size = 4) +
        ggtitle(paste0("TLPlevel = ", tlplevels[j], ", ",subsetting[i],
                       "\nSpecies Uptake Depth: Modelled Vs. Isotopic\n")) +
        ylab(expression("Modelled Uptake Depth (m)")) + xlab("Isotopic Uptake Depth (m)") +
        scale_x_reverse() + scale_y_reverse() +
        stat_poly_eq(aes(label = paste(..rr.label..)),
                     npcx = 0.57, npcy = 0.1, rr.digits = 2,
                     formula = formula, parse = TRUE, size = 6) +
        stat_fit_glance(method = 'lm',
                        method.args = list(formula = formula),
                        geom = 'text_npc',
                        aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                        npcx = 0.85, npcy = 0.1, size = 6)
      #scale_y_continuous(trans="rev_sqrt", breaks = c(0.00001, soil.depths))
      #scale_color_continuous(name = "Rsq\nMean", trans = "reverse")
      p1
      ggsave(file.path(file.path.udi.best.ll, paste0("Comparison_with_Meinzer1999_isotopic_vs.modelled_uptake.depth_cor",
                                                  goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 6, units ='in')
    }
  }

  ### udi.best.ll --------
  for (i in 1: 2) {
    if(i == 2) {
      iso.udi.i <- iso.udi %>% subset(!sp %in% c("guapst"))
    } else {
      iso.udi.i <- iso.udi
    }
    for (j in 1: length(tlplevels)) {
      iso.udi.sub <- iso.udi.i %>% subset(tlplevel == tlplevels[j])
      p0 <- ggplot(iso.udi.sub, aes(x =  Xylem_sap_deltaD_permil, y = udi.best.rsq, color = tlp)) +
        geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + SE, xmin = Xylem_sap_deltaD_permil - SE),
                       size = 0.5) +
        # geom_errorbar(aes(ymax = udi.best.rsq + udi.sd, ymin = udi.best.rsq - udi.sd),
        #               size = 0.5, width = 0.2) +
        scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
        geom_text(aes(x =  Xylem_sap_deltaD_permil + 2.5, y = udi.best.rsq + diff(range(iso.udi.sub$udi.best.rsq, na.rm = TRUE))/20, label = sp),
                  size = 4) +
        ggtitle(paste0("TLPlevel = ", tlplevels[j], ", ",subsetting[i], "\nSpecies Uptake Depth Vs Xylem Sap deltaD\n")) +
        ylab(expression("Uptake Depth (m)")) + xlab(xylem.label) + scale_y_reverse()
      p0 +
        geom_point(size = 3, show.legend = TRUE) +
        stat_poly_eq(aes(label = paste(..rr.label..)),
                     npcx = 0.6, npcy = 0.1, rr.digits = 2,
                     formula = formula, parse = TRUE, size = 4) +
        stat_fit_glance(method = 'lm',
                        method.args = list(formula = formula),
                        geom = 'text_npc',
                        aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                        npcx = 0.87, npcy = 0.1, size = 4)
      ggsave(file.path(file.path.udi.best.rsq, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
                                                      goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 6, units ='in')
      p0 +
        geom_point(aes(shape = Phenology), size = 3) +
        stat_poly_eq(aes(label = paste(..rr.label..)),
                     npcx = 0.6, npcy = 0.1, rr.digits = 2,
                     formula = formula, parse = TRUE, size = 4) +
        stat_fit_glance(method = 'lm',
                        method.args = list(formula = formula),
                        geom = 'text_npc',
                        aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                        npcx = 0.87, npcy = 0.1, size = 4)
      ggsave(file.path(file.path.udi.best.rsq, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_phenology_cor",
                                                      goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 8, units ='in')
      if(i == 1) {
        p0 + geom_point(size = 3, show.legend = TRUE) +
          geom_smooth(method = "lm", se = FALSE, lty = "dotted", color = "darkgray", size = 0.5) +
          stat_poly_eq(aes(label = paste(..rr.label..)),
                       npcx = 0.6, npcy = 0.1, rr.digits = 2,
                       formula = formula, parse = TRUE, size = 4, color = "darkgray") +
          stat_fit_glance(method = 'lm',
                          method.args = list(formula = formula),
                          geom = 'text_npc',
                          aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                          npcx = 0.88, npcy = 0.1, size = 4, color = "darkgray") +
          geom_smooth(data = subset(iso.udi.sub, !sp %in% c("sponra", "guapst")),
                      method = "lm", se = FALSE, color = "black", size = 0.5) +
          stat_poly_eq(data = subset(iso.udi.sub, !sp %in% c("sponra", "guapst")), aes(label = paste(..rr.label..)),
                       npcx = 0.6, npcy = 0.2, rr.digits = 2,
                       formula = formula, parse = TRUE, size = 4) +
          stat_fit_glance(data = subset(iso.udi.sub, !sp %in% c("sponra", "guapst")), method = 'lm',
                          method.args = list(formula = formula),
                          geom = 'text_npc',
                          aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                          npcx = 0.9, npcy = 0.2, size = 4)
      }
      ggsave(file.path(file.path.udi.best.rsq, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
                                                      goodness.fit, "_", tlplevels[j], "_", subsetting[i], "_lm.jpeg")), height = 5, width = 6, units ='in')
      ## without loncla
      if(j == 2) {
        p0 %+% subset(iso.udi.sub, !sp %in% c("loncla")) +
          geom_point(size = 3, show.legend = TRUE)
        ggsave(file.path(file.path.udi.best.rsq, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
                                                        goodness.fit, "_", tlplevels[j], "_", subsetting[i], "_witout loncla.jpeg")), height = 5, width = 6, units ='in')
      }
      ## only species with TLP data
      p0 %+% subset(iso.udi.sub, !is.na(tlp)) +
        geom_point(size = 3, show.legend = TRUE)
      ggsave(file.path(file.path.udi.best.rsq, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
                                                      goodness.fit, "_", tlplevels[j], "_", subsetting[i], "_only_TLP.jpeg")), height = 5, width = 6, units ='in')
      p0 %+% subset(iso.udi.sub, !is.na(tlp)) +
        geom_point(aes(shape = Phenology), size = 3)
      ggsave(file.path(file.path.udi.best.rsq, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_phenology_cor",
                                                      goodness.fit, "_", tlplevels[j], "_", subsetting[i], "_only_TLP.jpeg")), height = 5, width = 8, units ='in')

      p1 <- ggplot(iso.udi.sub, aes(x =  Xylem_sap_deltaD_permil, y = udi.best.rsq)) +
        geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + SE, xmin = Xylem_sap_deltaD_permil - SE),
                       size = 0.5) +
        # geom_errorbar(aes(ymax = udi.best + udi.sd, ymin = udi.best - udi.sd),
        #               size = 0.5, width = 0.2) +
        geom_text(aes(x =  Xylem_sap_deltaD_permil + 2.2, y = udi.best.rsq + diff(range(iso.udi.sub$udi.best.rsq, na.rm = TRUE))/20, label = sp),
                  size = 4) +
        stat_poly_eq(aes(label = paste(..rr.label..)),
                     npcx = 0.6, npcy = 0.1, rr.digits = 2,
                     formula = formula, parse = TRUE, size = 6, color = "darkgray") +
        stat_fit_glance(method = 'lm',
                        method.args = list(formula = formula),
                        geom = 'text_npc',
                        aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                        npcx = 0.9, npcy = 0.1, size = 6, color = "darkgray") +
        ylab(expression("Water Uptake Depth (m)")) + xlab(xylem.label) + scale_y_reverse() +
        geom_point(size = 3, show.legend = TRUE) + theme(text = element_text(size = 22), axis.text = element_text(size = 22)) +
        geom_smooth(method = "lm", se = FALSE, lty = "dotted", color = "darkgray")
      if( i == 1) {
        p1 +
          geom_smooth(data = subset(iso.udi.sub, !sp %in% c("sponra", "guapst")),
                      method = "lm", se = FALSE, color = "black") +
          stat_poly_eq(data = subset(iso.udi.sub, !sp %in% c("sponra", "guapst")), aes(label = paste(..rr.label..)),
                       npcx = 0.6, npcy = 0.2, rr.digits = 2,
                       formula = formula, parse = TRUE, size = 6) +
          stat_fit_glance(data = subset(iso.udi.sub, !sp %in% c("sponra", "guapst")), method = 'lm',
                          method.args = list(formula = formula),
                          geom = 'text_npc',
                          aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                          npcx = 0.9, npcy = 0.2, size = 6)

      } else {
        p1
      }
      ggsave(file.path(file.path.udi.best.rsq, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
                                                      goodness.fit, "_", tlplevels[j], "_", subsetting[i], "_no_color.jpeg")), height = 5, width = 5.25, units ='in')

      ### Isotopic vs. Update depth
      p1 <- ggplot(iso.udi.sub, aes(x =  depth, y = udi.best.rsq, color = tlp)) +
        geom_errorbarh(aes(xmax = depth + depth.se, xmin = depth - depth.se),
                       size = 0.5) +
        # geom_errorbar(aes(ymax = udi.best.rsq + udi.sd, ymin = udi.best.rsq - udi.sd),
        #               size = 0.5, width = 0.01) +
        geom_point(size = 3, show.legend = TRUE) +
        scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
        geom_text(aes(x = depth, y = udi.best.rsq + diff(range(iso.udi.sub$udi.best.rsq, na.rm = TRUE))/20, label = sp),
                  size = 4) +
        ggtitle(paste0("TLPlevel = ", tlplevels[j], ", ",subsetting[i],
                       "\nSpecies Uptake Depth: Modelled Vs. Isotopic\n")) +
        ylab(expression("Modelled Uptake Depth (m)")) + xlab("Isotopic Uptake Depth (m)") +
        scale_x_reverse() + scale_y_reverse() +
        stat_poly_eq(aes(label = paste(..rr.label..)),
                     npcx = 0.57, npcy = 0.1, rr.digits = 2,
                     formula = formula, parse = TRUE, size = 6) +
        stat_fit_glance(method = 'lm',
                        method.args = list(formula = formula),
                        geom = 'text_npc',
                        aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                        npcx = 0.85, npcy = 0.1, size = 6)
      #scale_y_continuous(trans="rev_sqrt", breaks = c(0.00001, soil.depths))
      #scale_color_continuous(name = "Rsq\nMean", trans = "reverse")
      p1
      ggsave(file.path(file.path.udi.best.rsq, paste0("Comparison_with_Meinzer1999_isotopic_vs.modelled_uptake.depth_cor",
                                                      goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 6, units ='in')
    }
  }

  # ### sdi.best --------
  # file.path.sdi.best.rsq <- file.path(file.path.sdi.rsq, "sdi.best.rsq")
  #
  # if(!dir.exists(file.path.sdi.best.rsq)) {dir.create(file.path.sdi.best.rsq)}
  # for (i in 1: 2) {
  #   if(i == 2) {
  #     iso.udi.i <- iso.udi %>% subset(!sp %in% c("guapst")) #c("sponra", "guapst")
  #   } else {
  #     iso.udi.i <- iso.udi
  #   }
  #   for (j in 1: length(tlplevels)) {
  #     iso.udi.sub <- iso.udi.i %>% subset(tlplevel == tlplevels[j])
  #     p0 <- ggplot(iso.udi.sub, aes(x =  Xylem_sap_deltaD_permil, y = sdi.med.rsq, color = tlp)) +
  #       geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + SE, xmin = Xylem_sap_deltaD_permil - SE),
  #                      size = 0.5) +
  #       geom_errorbar(aes(ymax = sdi.best.rsq + sdi.sd.rsq, ymin = sdi.best.rsq - sdi.sd.rsq),
  #                     size = 0.5, width = 0.2) +
  #       scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
  #       geom_text(aes(x =  Xylem_sap_deltaD_permil + 1.5, y = sdi.best.rsq + diff(range(iso.udi.sub$sdi.best.rsq, na.rm = TRUE))/20, label = sp),
  #                 size = 4) +
  #       stat_poly_eq(aes(label = paste(..rr.label..)),
  #                    npcx = 0.6, npcy = 0.1, rr.digits = 2,
  #                    formula = formula, parse = TRUE, size = 6) +
  #       stat_fit_glance(method = 'lm',
  #                       method.args = list(formula = formula),
  #                       geom = 'text_npc',
  #                       aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
  #                       npcx = 0.85, npcy = 0.1, size = 6) +
  #       ggtitle(paste0("TLPlevel = ", tlplevels[j], ", ",subsetting[i], "\nSpecies Uptake Depth Vs Xylem Sap deltaD\n")) +
  #       ylab(expression("Water-Stress Depth (m)")) + xlab(xylem.label) +
  #       # scale_y_continuous(trans="rev_sqrt", breaks = c(0.00001, soil.depths))
  #       scale_y_reverse()
  #     p0 +
  #       geom_point(size = 3, show.legend = TRUE)
  #     #scale_color_continuous(name = "Rsq\nMean", trans = "reverse")
  #     ggsave(file.path(file.path.sdi.best, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
  #                                                 goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 6, units = 'in')
  #   }
  # }
}
