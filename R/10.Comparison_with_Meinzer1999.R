## Comparison with Meinzer et al. 1999 results.

#
rm(list=ls())
gc()
library(tidyverse)
# graphics info
theme_set(theme_bw())
theme_update(text = element_text(size = 14),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)
require(scales)
rev_sqrt_trans <- function() {
  scales::trans_new(
    name = "rev_sqrt",
    transform = function(x) -sqrt(abs(x)),
    inverse = function(x) x^2);
}

# load interval and working.iter
load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
load(file.path("results/4.1GLUEsetup_part2_BCI.RData")) # has working.iter and growth and si matrix

intervals <- info$intervals
n.ensembles <- growth_by_si.info$n.ensembles
growth.type <- growth_by_si.info$growth.type
growth.selection <- growth_by_si.info$growth.selection
dbh.residuals <- growth_by_si.info$dbh.residuals
si.type <- growth_by_si.info$si.type
goodness.fit <- 0.3
soil.depths <- unique(info$root.param.long$depth)
dryseason <- "on"
##
load(file = paste("results/splevel/ds.bestfit_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_id_dryseason_", dryseason, ".Rdata", sep = ""))
ds <- ds.bestfit
load(file = paste("results/commlevel/ds.bestfit_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_id_dryseason_", dryseason, ".Rdata", sep = ""))
ds <- rbind(ds, ds.bestfit) %>% mutate(tlplevel = as.factor(tlplevel)) %>% subset(!is.na(udi)) %>% droplevels()

iso <- read.csv("data-raw/traits/isotopes/Meinzer1999_Table1.csv", na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)

head(iso)
colnames(iso)[1] <- "species"
# adding species codes
load(file.path("data-raw/CTFScensuses/bci.spptable.rdata"))
write.csv(bci.spptable, file.path("data-raw/CTFScensuses/bci.spptable.csv"), row.names = FALSE)

iso$genus.sp <- paste(iso$species, iso$Family, sep = " ")
bci.spptable$genus.sp <- paste(bci.spptable$Genus, bci.spptable$Species, sep = " ")
matchrows = match(iso$genus.sp, bci.spptable$genus.sp)
iso$genus.sp

iso$sp <- bci.spptable$sp[matchrows]

iso$sp
## some codes are not found, because of change of latin names, so searching in synonymous name-field should return the match
match("Lonchocarpus", bci.spptable$Genus) #"Lonchocarpus latifolius"

grep("Lonchocarpus", bci.spptable$syn)
iso$sp[which(iso$species == "Lonchocarpus")] <- bci.spptable$sp[grep("Lonchocarpus", bci.spptable$syn)]
head(iso)

#"Guapira" ; iso$Family[11]
grep("Guapira", bci.spptable$syn)
iso$sp[which(iso$species == "Guapira")] <- bci.spptable$sp[grep("Guapira", bci.spptable$syn)]# "Guapira" "campochagres"
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

iso.soil.1 <- read.csv("data-raw/traits/isotopes/Meinzer1999_Fig2A_soil_deltaD_BCI_data_Mar_Apr_1997.csv", na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)
iso.soil.2 <- read.csv("data-raw/traits/isotopes/Oecologia 1995 Jackson_fig2_soil_deltaD_Gigante_data_dry_season_1992.csv", na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)
depth.m1 <- lm(depth ~ soil.deltaD, data = iso.soil.1)
summ.depth.m1 <- summary(depth.m1)
depth.m1.label = paste0("Depth = ", round(depth.m1$coefficients[1], 2), " + ",
                         round(depth.m1$coefficients[2], 2),
                         " * Soil deltaD\nR-squared = ", round(summ.depth.m1$r.squared, 2),
                         ", p-val = ", round(summ.depth.m1$coefficients[2, 4], 4))

g1 <- ggplot(iso.soil.1, aes(y = depth, x = soil.deltaD)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  xlab(expression("Soil deltaD")) + ylab("Depth (cm)") +
  scale_y_reverse()
g1 + geom_errorbarh(aes(xmax = soil.deltaD + se, xmin = soil.deltaD - se), size = 0.5) +
geom_text(aes(x = -40, y = 0, label = depth.m1.label), color = "blue")

ggsave(file.path(paste0("figures/UDI_confidence/Meinzer_etal_1999_Depth_vs_soil_deltaD.jpeg")), height = 5, width = 5, units='in')

depth.m2 <- lm(depth ~ soil.deltaD, data = iso.soil.2)
summ.depth.m2 <- summary(depth.m2)
depth.m2.label = paste0("Depth = ", round(depth.m2$coefficients[1], 2), " + ",
                        round(depth.m2$coefficients[2], 2),
                        " * Soil deltaD\nR-squared = ", round(summ.depth.m2$r.squared, 2),
                        ", p-val = ", round(summ.depth.m2$coefficients[2, 4], 4))
g1 %+% iso.soil.2 + geom_smooth(method = "glm", se = FALSE) +
  geom_text(aes(x = -35, y = 5, label = depth.m2.label), color = "blue")
ggsave(file.path(paste0("figures/UDI_confidence/Jackson_etal_1995_Depth_vs_soil_deltaD.jpeg")), height = 5, width = 5, units='in')

## how to incorporate se in this?
iso <- iso %>% mutate(depth = predict.lm(depth.m1, newdata = data.frame(soil.deltaD = Xylem_sap_deltaD_permil))/100,
                      depth.se = depth - predict.lm(depth.m1,
                                                    newdata = data.frame(soil.deltaD = Xylem_sap_deltaD_permil + SE))/100) # from cm to m
head(iso)

udi.all <- ds
udi <- ds %>%
  subset(!is.na(udi) & size %in% c("large")) %>% droplevels() #%>%
  # group_by(sp, tlplevel) %>%
  # summarise_at(vars(udi, udi.se, sdi, sdi.se, rsq.mean, rsq.min), mean, na.rm = TRUE)
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
nrow(iso.udi)
iso.sp <- unique(iso.udi$sp[!is.na(iso.udi$udi) & !is.na(iso.udi$Xylem_sap_deltaD_permil)])
save(iso.sp, file = "results/sp_with_isotopic_record.Rdata")
iso.udi <- iso.udi[order(iso.udi$udi),]
tlplevels <- c("sp", "comm")
subsetting <- c("with outliers", "without outliers")
View(iso.udi)
write.csv(iso.udi, file.path(paste0("results/iso.udi_cor", goodness.fit, "_", si.type, "_",
                                    n.ensembles, "_", growth.type, "_", intervals, "_id_dryseason_", dryseason, ".csv")), row.names = FALSE)

file.path.udi <- file.path("figures/UDI_confidence", growth.type, growth.selection, paste0("dbh.residuals_", dbh.residuals))
if(!dir.exists(file.path.udi)) {dir.create(file.path.udi)}

# iso <- read.csv("data-raw/traits/isotopes/Oecologia 1995 Jackson _Fig3_Fig4.csv", na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)
# iso.udi <- left_join(iso, udi, by = "sp")
# iso.sp <- unique(iso.udi$sp[!is.na(iso.udi$udi) & !is.na(iso.udi$Xylem_sap_deltaD_permil)])
# iso.udi <- iso.udi[order(iso.udi$udi),]

for (i in 1: 2) {
  if(i == 2) {
    iso.udi.i <- iso.udi %>% subset(!sp %in% c("sponra", "guapst"))
  } else {
    iso.udi.i <- iso.udi
  }
  for (j in 1: length(tlplevels)) {
    iso.udi.sub <- iso.udi.i %>% subset(tlplevel == tlplevels[j])
    m2 <- lm(udi ~ Xylem_sap_deltaD_permil, data = iso.udi.sub)
    rank.corr <- cor.test(x=iso.udi.sub$udi, y = -iso.udi.sub$Xylem_sap_deltaD_permil,
                          method = 'spearman', exact = FALSE)
    summ.m2 <- summary(m2)
    m2.label = paste0("R-squared = ", round(summ.m2$r.squared, 3),
                      ", p-val = ", round(summ.m2$coefficients[2, 4], 2),
                      "\n Spearman's Rho = ", round(as.numeric(rank.corr$estimate), 2), ", p-val = ", round(as.numeric(rank.corr$p.value), 2))

    p0 <- ggplot(iso.udi.sub, aes(x =  Xylem_sap_deltaD_permil, y = udi)) +
      # geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + SE, xmin = Xylem_sap_deltaD_permil - SE),
      #                colour = "darkgray", size = 0.5) +
      geom_errorbar(aes(ymax = udi + udi.se, ymin = udi - udi.se),
                    colour = "darkgray", size = 0.5, width = 0.2) +
      scale_shape_manual(values = c(21, 24)) +
      geom_text(aes(x =  Xylem_sap_deltaD_permil, y = udi + 0.12, label = sp),
                colour = "black", size = 4) +
      ggtitle(paste0("TLPlevel = ", tlplevels[j], ", ",subsetting[i], "\nSpecies Uptake Depth Vs Xylem Sap deltaD")) +
      ylab(expression("Uptake Depth (m)")) + xlab("Xylem Sap deltaD (permil)") +
      scale_y_reverse() +
      geom_text(aes(x = -35, y = max(iso.udi.sub$udi, na.rm = TRUE) - 0.2, label = m2.label), color = "blue")
    p0 +
      geom_point(size = 3, show.legend = TRUE)
    ggsave(file.path(file.path.udi, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
                            goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 5, units ='in')
    p0 +
      geom_point(aes(fill = Phenology), size = 3, shape = 21, colour = "white")
    ggsave(file.path(file.path.udi, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_phenology_cor",
                            goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 8, units ='in')

    ### Isotopic vs. Update depth
    m3 <- lm(udi ~ depth, data = iso.udi.sub)
    summ.m3 <- summary(m3)
    m3.label = paste0("R-squared = ", round(summ.m3$r.squared, 3),
                      "\np-val = ", round(summ.m3$coefficients[2, 4], 2))
    p1 <- ggplot(iso.udi.sub, aes(x =  depth, y = udi)) +
      geom_text(aes(x = 1, y = max(iso.udi.sub$udi, na.rm = TRUE) + 0.2, label = m3.label), color = "blue") +
      geom_errorbarh(aes(xmax = depth, xmin = depth),
                     colour = "darkgray", size = 0.5) +
      geom_errorbar(aes(ymax = udi + udi.se, ymin = udi - udi.se),
                    colour = "darkgray", size = 0.5, width = 0.01) +
      geom_point(size = 3, show.legend = TRUE) +
      scale_shape_manual(values = c(21, 24)) +
      geom_text(aes(x = depth, y = udi + 0.12, label = sp),
                colour = "black", size = 4) +
      ggtitle(paste0("TLPlevel = ", tlplevels[j], ", ",subsetting[i],
                     "\nSpecies Uptake Depth: Modelled Vs. Isotopic")) +
      ylab(expression("Modelled Uptake Depth (m)")) + xlab("Isotopic Uptake Depth (m)") +
      scale_x_reverse() + scale_y_reverse()
      #scale_color_continuous(name = "Rsq\nMean", trans = "reverse")
    p1
    ggsave(file.path(file.path.udi, paste0("Comparison_with_Meinzer1999_isotopic_vs.modelled_uptake.depth_cor",
                            goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 6, units ='in')
  }
}

### sdi
file.path.sdi <- file.path("figures/UDI_confidence", growth.type, growth.selection, paste0("dbh.residuals_", dbh.residuals), "sdi")
if(!dir.exists(file.path.sdi)) {dir.create(file.path.sdi)}
for (i in 1: 2) {
  if(i == 2) {
    iso.udi.i <- iso.udi %>% subset(!sp %in% c("sponra", "guapst"))
  } else {
    iso.udi.i <- iso.udi
  }
  for (j in 1: length(tlplevels)) {
    iso.udi.sub <- iso.udi.i %>% subset(tlplevel == tlplevels[j])
    m2 <- lm(udi ~ Xylem_sap_deltaD_permil, data = iso.udi.sub)
    rank.corr <- cor.test(x=iso.udi.sub$sdi, y = -iso.udi.sub$Xylem_sap_deltaD_permil, method = 'spearman')
    summ.m2 <- summary(m2)
    m2.label = paste0("R-squared = ", round(summ.m2$r.squared, 3),
                      ", p-val = ", round(summ.m2$coefficients[2, 4], 2),
                      "\n Spearman's Rho = ", round(as.numeric(rank.corr$estimate), 2), ", p-val = ", round(as.numeric(rank.corr$p.value), 2))

    p0 <- ggplot(iso.udi.sub, aes(x =  Xylem_sap_deltaD_permil, y = sdi)) +
      geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + SE, xmin = Xylem_sap_deltaD_permil - SE),
                     colour = "darkgray", size = 0.5) +
      geom_errorbar(aes(ymax = sdi + sdi.se, ymin = sdi - sdi.se),
                    colour = "darkgray", size = 0.5, width = 0.2) +
      scale_shape_manual(values = c(21, 24)) +
      geom_text(aes(x =  Xylem_sap_deltaD_permil, y = sdi + 0.12, label = sp),
                colour = "black", size = 4) +
      ggtitle(paste0("TLPlevel = ", tlplevels[j], ", ",subsetting[i], "\nSpecies Uptake Depth Vs Xylem Sap deltaD")) +
      ylab(expression("Water-Stress Depth (m)")) + xlab("Xylem Sap deltaD (permil)") +
      scale_y_reverse() +
      geom_text(aes(x = -55, y = 0.2, label = m2.label), color = "blue")
    p0 +
      geom_point(size = 3, show.legend = TRUE)
      #scale_color_continuous(name = "Rsq\nMean", trans = "reverse")
    ggsave(file.path(file.path.sdi, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
                            goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 5, units ='in')
  }
}


###### Best
file.path.udi.best <- file.path("figures/UDI_confidence", growth.type, growth.selection, paste0("dbh.residuals_", dbh.residuals), "udi.best")
if(!dir.exists(file.path.udi.best)) {dir.create(file.path.udi.best)}


for (i in 1: 2) {
  if(i == 2) {
    iso.udi.i <- iso.udi %>% subset(!sp %in% c("sponra", "guapst"))
  } else {
    iso.udi.i <- iso.udi
  }
  for (j in 1: length(tlplevels)) {
    iso.udi.sub <- iso.udi.i %>% subset(tlplevel == tlplevels[j])
    m2 <- lm(udi.best ~ Xylem_sap_deltaD_permil, data = iso.udi.sub)
    rank.corr <- cor.test(x=iso.udi.sub$udi.best, y = -iso.udi.sub$Xylem_sap_deltaD_permil,
                          method = 'spearman', exact = FALSE)
    summ.m2 <- summary(m2)
    m2.label = paste0("R-squared = ", round(summ.m2$r.squared, 3),
                      ", p-val = ", round(summ.m2$coefficients[2, 4], 2),
                      "\n Spearman's Rho = ", round(as.numeric(rank.corr$estimate), 2), ", p-val = ", round(as.numeric(rank.corr$p.value), 2))

    p0 <- ggplot(iso.udi.sub, aes(x =  Xylem_sap_deltaD_permil, y = udi.best)) +
      # geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + SE, xmin = Xylem_sap_deltaD_permil - SE),
      #                colour = "darkgray", size = 0.5) +
      geom_errorbar(aes(ymax = udi.best + udi.se, ymin = udi.best - udi.se),
                    colour = "darkgray", size = 0.5, width = 0.2) +
      scale_shape_manual(values = c(21, 24)) +
      geom_text(aes(x =  Xylem_sap_deltaD_permil, y = udi.best + 0.12, label = sp),
                colour = "black", size = 4) +
      ggtitle(paste0("TLPlevel = ", tlplevels[j], ", ",subsetting[i], "\nSpecies Uptake Depth Vs Xylem Sap deltaD")) +
      ylab(expression("Uptake Depth (m)")) + xlab("Xylem Sap deltaD (permil)") +
      scale_y_reverse() +
      geom_text(aes(x = -35, y = max(iso.udi.sub$udi.best, na.rm = TRUE) - 0.2, label = m2.label), color = "blue")
    p0 +
      geom_point(size = 3, show.legend = TRUE)
    ggsave(file.path(file.path.udi.best, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
                                           goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 5, units ='in')
    p0 +
      geom_point(aes(fill = Phenology), size = 3, shape = 21, colour = "white")
    ggsave(file.path(file.path.udi.best, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_phenology_cor",
                                           goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 8, units ='in')

    ### Isotopic vs. Update depth
    m3 <- lm(udi.best ~ depth, data = iso.udi.sub)
    summ.m3 <- summary(m3)
    m3.label = paste0("R-squared = ", round(summ.m3$r.squared, 3),
                      "\np-val = ", round(summ.m3$coefficients[2, 4], 2))
    p1 <- ggplot(iso.udi.sub, aes(x =  depth, y = udi.best)) +
      geom_text(aes(x = 1, y = max(iso.udi.sub$udi.best, na.rm = TRUE) + 0.2, label = m3.label), color = "blue") +
      geom_errorbarh(aes(xmax = depth, xmin = depth),
                     colour = "darkgray", size = 0.5) +
      geom_errorbar(aes(ymax = udi.best + udi.se, ymin = udi.best - udi.se),
                    colour = "darkgray", size = 0.5, width = 0.01) +
      geom_point(size = 3, show.legend = TRUE) +
      scale_shape_manual(values = c(21, 24)) +
      geom_text(aes(x = depth, y = udi.best + 0.12, label = sp),
                colour = "black", size = 4) +
      ggtitle(paste0("TLPlevel = ", tlplevels[j], ", ",subsetting[i],
                     "\nSpecies Uptake Depth: Modelled Vs. Isotopic")) +
      ylab(expression("Modelled Uptake Depth (m)")) + xlab("Isotopic Uptake Depth (m)") +
      scale_x_reverse() + scale_y_reverse()
    #scale_color_continuous(name = "Rsq\nMean", trans = "reverse")
    p1
    ggsave(file.path(file.path.udi.best, paste0("Comparison_with_Meinzer1999_isotopic_vs.modelled_uptake.depth_cor",
                                           goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 6, units ='in')
  }
}

### sdi
file.path.sdi.best <- file.path("figures/UDI_confidence", growth.type, growth.selection, paste0("dbh.residuals_", dbh.residuals), "sdi.best")

if(!dir.exists(file.path.sdi)) {dir.create(file.path.sdi)}
for (i in 1: 2) {
  if(i == 2) {
    iso.udi.i <- iso.udi %>% subset(!sp %in% c("sponra", "guapst"))
  } else {
    iso.udi.i <- iso.udi
  }
  for (j in 1: length(tlplevels)) {
    iso.udi.sub <- iso.udi.i %>% subset(tlplevel == tlplevels[j])
    m2 <- lm(udi ~ Xylem_sap_deltaD_permil, data = iso.udi.sub)
    rank.corr <- cor.test(x=iso.udi.sub$sdi.best, y = -iso.udi.sub$Xylem_sap_deltaD_permil, method = 'spearman')
    summ.m2 <- summary(m2)
    m2.label = paste0("R-squared = ", round(summ.m2$r.squared, 3),
                      ", p-val = ", round(summ.m2$coefficients[2, 4], 2),
                      "\n Spearman's Rho = ", round(as.numeric(rank.corr$estimate), 2), ", p-val = ", round(as.numeric(rank.corr$p.value), 2))

    p0 <- ggplot(iso.udi.sub, aes(x =  Xylem_sap_deltaD_permil, y = sdi)) +
      geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + SE, xmin = Xylem_sap_deltaD_permil - SE),
                     colour = "darkgray", size = 0.5) +
      geom_errorbar(aes(ymax = sdi.best + sdi.se, ymin = sdi.best - sdi.se),
                    colour = "darkgray", size = 0.5, width = 0.2) +
      scale_shape_manual(values = c(21, 24)) +
      geom_text(aes(x =  Xylem_sap_deltaD_permil, y = sdi.best + 0.12, label = sp),
                colour = "black", size = 4) +
      ggtitle(paste0("TLPlevel = ", tlplevels[j], ", ",subsetting[i], "\nSpecies Uptake Depth Vs Xylem Sap deltaD")) +
      ylab(expression("Water-Stress Depth (m)")) + xlab("Xylem Sap deltaD (permil)") +
      scale_y_reverse() +
      geom_text(aes(x = -55, y = 0.2, label = m2.label), color = "blue")
    p0 +
      geom_point(size = 3, show.legend = TRUE)
    #scale_color_continuous(name = "Rsq\nMean", trans = "reverse")
    ggsave(file.path(file.path.sdi.best, paste0("Comparison_with_Meinzer1999_deltaD_vs.modelled_uptake.depth_cor",
                                           goodness.fit, "_", tlplevels[j], "_", subsetting[i], ".jpeg")), height = 5, width = 5, units ='in')
  }
}
