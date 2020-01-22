##
rm(list = ls())
gc()
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, scales)

# graphics info
theme_set(theme_bw())
theme_update(text = element_text(size = 14),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)

# load interval and n.ensembles
load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
load(file.path("results/4.1GLUEsetup_part2_BCI.RData")) # has n.ensembles and growth and si matrix

intervals <- info$intervals
n.ensembles <- growth_by_si.info$n.ensembles
growth.type <- growth_by_si.info$growth.type
growth.selection <- growth_by_si.info$growth.selection
dbh.residuals <- growth_by_si.info$dbh.residuals
si.type <- growth_by_si.info$si.type
runs <- info$root.param
goodness.fit <- 0.3
dryseason = "on"
soil.depths <- unique(info$root.param.long$depth)
if(!dir.exists(file.path("figures/UDI", growth.type))) {dir.create(file.path("figures/UDI", growth.type))}

##
load(file = paste("results/splevel/ds.bestfit_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_id_dryseason_", dryseason, ".Rdata", sep = ""))
load(file = paste("results/splevel/ds.bestfit.longer_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_id_dryseason_", dryseason, ".Rdata", sep = ""))
ds <- ds.bestfit
dss <- ds.bestfit.longer
load(file = paste("results/commlevel/ds.bestfit_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_id_dryseason_", dryseason, ".Rdata", sep = ""))
load(file = paste("results/commlevel/ds.bestfit.longer_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_id_dryseason_", dryseason, ".Rdata", sep = ""))

ds <- rbind(ds, ds.bestfit) %>% mutate(tlplevel = as.factor(tlplevel))
dss <- rbind(dss, ds.bestfit.longer) %>% mutate(tlplevel = as.factor(tlplevel))
length(ds$sp_size[!is.na(ds$udi)]) # 308 for cor.03
# #dss <- read.csv(paste("results/StressIndex.modelfits_", n.ensembles, growth.type,".csv", sep = ""), na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)
# # load(file = paste("results/ds.bestfit_", "si_light", "_", n.ensembles, "_", growth.type, "_", intervals, "_id.Rdata", sep = ""))
# # ds.bestfit1 <- ds.bestfit %>% mutate(si.type = "si_light")
# # load(file = paste("results/ds.bestfit_", "si_mean", "_", n.ensembles, "_", growth.type, "_", intervals, "_id.Rdata", sep = ""))
# ds.bestfit2 <- ds.bestfit %>% mutate(si.type = "si_mean")
# ds.bestfit.both <- bind_rows(ds.bestfit1, ds.bestfit2)
# ds.bestfit.both %>% group_by(si.type) %>% summarise(mean(rsq, na.rm = T))
# ggplot(ds.bestfit.both %>% subset(size == "large"),
#        aes(x = sp, y = -UDI, colour = rsq, shape = si.type)) +
#   geom_point(size = 4) + ylim(c(-3, 0)) + theme_bw() +
#   theme(legend.position = "top") + my.theme +
#   scale_color_continuous(name = "Model", trans = "reverse") +
#   theme(axis.text.x = element_text(size = 15, face = "plain", angle = 45, vjust = 1, hjust = 1)) +
#   ylab(expression("Uptake Depth Index")) + xlab("sp_size")
# ggsave(file.path(paste("figures/UDI/", growth.type,"/BCI/UDI_confidence/UDI_si_light_vs_SI_mean",
#                        growth.type, "_#ewlQti.jpeg")), height = 8.94, width = 12.8, units = 'in')

# ds.bestfit <- ds.bestfit %>% filter(!is.na(sp) & rsq >= 0.3)

# load(file = paste("results/best10.type.rsq_", si.type, "_", n.ensembles, "_", growth.type, ".Rdata", sep = ""), envir = parent.frame(), verbose = FALSE)
# load(file = paste("results/best10.type_", si.type, "_", n.ensembles, "_", growth.type, ".Rdata", sep = ""), envir = parent.frame(), verbose = FALSE)
# View(dss %>% group_by(size) %>% summarise_all(funs(mean(., na.rm = T))))
# dat.sp <- read.csv("results/dat.sp.csv", na.strings = c("NA",""), header = T, row.names = NULL, check.names = F)
# dss <- subset(dss, N >= 5)
# dss <- droplevels(dss)
# mean. Rsq Vs Uptake
require(scales)
rev_sqrt_trans <- function() {
  scales::trans_new(
    name = "rev_sqrt",
    transform = function(x) -sqrt(abs(x)),
    inverse = function(x) x^2);
}
dss.long.sub <- dss %>% # %>% subset(!is.na(root.frac))
  select(udi, depth, cum.root.frac, root.frac, rf.sam,
         sp, sp_size, size, sp_size_par.sam_n.best, n, max.root, tlplevel, rsq) %>%
  mutate(depth = as.numeric(depth)) %>% arrange(desc(n))

sp_n <- unique(select(dss, sp_size, n)) %>% arrange(desc(n))
dss.long.sub.large.2 <- dss.long.sub %>%
  subset(size == "large") %>%
  transform(sp = reorder(sp, udi))

ggplot(dss.long.sub,
       aes(y = depth, x = root.frac)) + xlab("Root Fraction") + ylab("Depth (m)") +
  geom_line(aes(group = sp_size_par.sam_n.best, color = sp), show.legend = FALSE) +
  facet_grid(size ~ tlplevel) +
  scale_y_continuous(trans="rev_sqrt", breaks = soil.depths)
ggsave(paste0("figures/UDI/", growth.type,"/root.frac_by_size_cor", goodness.fit, "_", si.type, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, ".jpeg"), width = 7, height = 7, units = 'in')

### for species with isotopic records
load(file = "results/sp_with_isotopic_record.Rdata")
dss.long.sub.isosp <- dss.long.sub %>% subset(sp %in% iso.sp & size %in% c("large")) %>%
  droplevels() %>% arrange(udi)
# ggplot(dss.long.sub.isosp,
#        aes(y = depth, x = cum.root.frac)) + xlab("Cumulative Root Fraction") + ylab("Depth (m)") +
#   geom_line(aes(group = sp_size_par.sam_n.best, color = rsq)) +
#   facet_wrap(sp ~ .) +
#   scale_y_continuous(trans="rev_sqrt", breaks = soil.depths) +
#   scale_color_continuous(trans = "reverse")
# ggsave(paste0("figures/UDI/", growth.type,"/cum.root.frac_iso.sp_cor", goodness.fit, "_", si.type, "_", growth.type, "_", intervals, ".jpeg"), width = 9, height = 9, units = 'in')

best.dss.long.sub.isosp <- dss.long.sub.isosp %>%
  unite("sp_size.tlplevel", sp_size, tlplevel, remove = FALSE) %>% group_by(sp_size.tlplevel) %>%
  mutate(max.rsq = rsq == max(rsq, na.rm = TRUE)) %>% subset(max.rsq == TRUE)
tlplevels <- c("sp", "comm")
rsq.thresh <- 0.9
for (i in 1: length(tlplevels)) {
  ggplot(dss.long.sub.isosp %>% subset(rsq >= rsq.thresh & tlplevel == tlplevels[i]),
         aes(y = depth, x = cum.root.frac)) + xlab("Cumulative Root Fraction") + ylab("Depth (m)") +
    geom_line(aes(group = sp_size_par.sam_n.best, color = rsq)) +
    geom_line(data = best.dss.long.sub.isosp %>% subset(tlplevel == tlplevels[i]),
              aes(group = sp_size_par.sam_n.best), color = "red", size = 1) +
    facet_wrap(sp ~ .) +
    scale_y_continuous(trans="rev_sqrt", breaks = c(0.00001, soil.depths)) +
    scale_color_continuous(trans = "reverse")
  ggsave(paste0("figures/UDI/", growth.type,"/cum.root.frac_iso.sp_rsq", rsq.thresh ,"_cor", goodness.fit, "_",
                tlplevels[i], "_", si.type, "_", growth.type, "_", intervals,".jpeg"), width = 9, height = 9, units = 'in')
}

ggplot(dss.long.sub.isosp,
       aes(y = depth, x = root.frac)) + xlab("Root Fraction") + ylab("Depth (m)") +
  geom_line(aes(group = sp_size_par.sam_n.best, color = sp), show.legend = FALSE) +
  facet_wrap(sp ~ .) +
  scale_y_continuous(trans="rev_sqrt", breaks = soil.depths)
ggsave(paste0("figures/UDI/", growth.type,"/root.frac_iso.sp_cor", goodness.fit, "_", si.type, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, ".jpeg"), width = 7, height = 7, units = 'in')

ggplot(dss.long.sub,
       aes(y = depth, x = cum.root.frac)) + xlab("Cumulative Root Fraction") + ylab("Depth (m)") +
  geom_line(aes(group = sp_size_par.sam_n.best, color = sp), show.legend = FALSE) +
  facet_grid(size ~ tlplevel) +
  scale_y_continuous(trans="rev_sqrt", breaks = soil.depths)
ggsave(paste0("figures/UDI/", growth.type,"/cum.root.frac_by_size_cor", goodness.fit, "_", si.type, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, ".jpeg"), width = 7, height = 9, units = 'in')

dss <- transform(dss, sp = reorder(sp, -udi))
# this plots rsq for individual  sp_size_par.sam_n.best
ggplot(dss, aes(x = udi, y = rsq, colour = size)) +
  geom_point(size = 2, alpha = 0.7) +
  facet_grid(. ~ tlplevel) +
  guides(colour = guide_legend(title = "sp_size")) +
  ggtitle("Size class: Rsq Vs Uptake Depth Index") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, by = 0.2)) +
  xlab(expression("Uptake Depth Index (m)")) + ylab("Rsq")
ggsave(paste0("figures/UDI/", growth.type,"/UDI_rsq_size_class_cor", goodness.fit, "_", si.type, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_rsq_vs_UDI.jpeg"), width = 7, height = 5, units = 'in')
# this plots mean rsq for sp_size
ggplot(ds, aes(x = udi, y = rsq.mean, colour = size)) +
  geom_point(size = 2, alpha = 0.7) +
  facet_grid(. ~ tlplevel) +
  guides(colour = guide_legend(title = "sp_size")) +
  ggtitle("Size class: Rsq Vs Uptake Depth Index") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, by = 0.2)) +
  xlab(expression("Uptake Depth Index (m)")) + ylab("Mean Rsq")
ggsave(paste0("figures/UDI/", growth.type,"/UDI_mean.rsq_size_class_cor", goodness.fit, "_", si.type, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_rsq_vs_UDI.jpeg"), width = 7, height = 5, units = 'in')

range(ds$n, na.rm = TRUE)
# 9 56512 # 5 835
select.sp <- unique(ds$sp)
# save(select.sp_size, file = paste("results/sp_size_>Rsq0.5_>5N_non-borderline", n.ensembles, growth.type, ".Rdata", sep = ""))
## making a sp_size names table
load(file = paste("data-raw/CTFScensuses/bci.spptable.Rdata", sep = ""), envir = parent.frame(), verbose = FALSE)
head(bci.spptable)
sp.info <- subset(bci.spptable, sp %in% select.sp, select = -c(wsg, wsglevel))
head(sp.info)
ds.sub <- subset(ds, select = c(sp_size, sp, size, udi, udi.se, n))
# colnames(ds.sub) <- c("sp", "UDmean_best_fit_of_top10_parameter_ensembles", "UDmean_mean__of_top10_parameter_ensembles", "UDmean_SE_of_top10_parameter_ensembles", "UD.demand", "sample_size")
ds.sub$sp <- as.character(ds.sub$sp)
library(dplyr)
sp.info.ud <- left_join(sp.info, ds.sub, by = "sp")
head(sp.info.ud)
# write.csv(sp.info.ud, file = "results/sp_size_Names_UDmean_#XKGnjY.csv", row.names = FALSE)
# write.table(sp.info.ud, file = "results/sp_size_Names_UDmean_#XKGnjY.txt", row.names = FALSE)
string(1, 6)
# only  one sp_size has border line equifinality issue
nrow(ds) # 493 # 57
nrow(dss) # 580 # 65 sp_size
range(dss$n)
# 5 835
# with dss$best.type.rsq >= 0.7; N > 5 -> 53 sp,
# with dss$best.type.rsq >= 0.8; N > 5 -> 46 sp,
# with dss$best.type.rsq >= 0.9; N > 5 -> 31 sp,
ds <- ds %>% transform(sp = reorder(sp, udi)) %>%
  mutate(size = factor(size, levels = c("tiny", "small", "medium", "large"))) %>%
  droplevels()
## to mantain same color for each sp_size
# library(scales)
# spnum <- length(unique(ds$sp_size))
# myColors <- hue_pal()(spnum) # t ovisualise try: show_col(hue_pal()(spnum))
# names(myColors) <- levels(ds$sp_size)
# fillScale <- scale_colour_manual(name = "sp_size", values = myColors)
# huh doesnt seem to control colours after sp_size order transformed?
###-------to plot few sp_size, selecting sp with N > 100------
sp_n <- unique(select(ds, sp_size, n)) %>% arrange(desc(n))
#few_sp <- sp_n$sp_size[1:20]
# dss <- subset(dss, sp_size %in% few_sp)
# ds <- subset(ds, sp_size %in% few_sp)
#growth.type <- paste(growth.type, "top_sp", sep = "_")
###---------------
ds.large <- ds %>% subset(size == "large" & !is.na(udi)) %>% transform(sp = reorder(sp, udi)) %>% droplevels()
g <- ggplot(ds.large, aes(x = as.factor(sp), y = -udi)) +
  geom_errorbar(aes(ymax = -udi + udi.se, ymin = -udi - udi.se), colour = "black", width = 0.001, size = 0.5) +
  #geom_point(colour = "black", fill = "white", aes(shape = tlplevel, size = n)) + # ,
  scale_size(range = c(1, 5), "Sample Size", breaks = c(10, 100, 500)) +
  theme(legend.position = c(0.85, 0.75)) +
  scale_y_continuous(trans="rev_sqrt", breaks = soil.depths) +
  theme(axis.text.x = element_text(size = 10, face = "plain", angle = 45, vjust = 1, hjust = 1)) +
  xlab("Species") + ylab(expression("Uptake Depth Index (m)"))
g0 <- g + geom_point(colour = "darkgray", aes(shape = tlplevel, fill = sp, size = n)) +
  scale_shape_manual(values = c(21, 24)) + guides(fill = "none") +
  theme(axis.text.x = element_text(size = 8, face = "plain", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12))
ggsave(plot = g0, file.path(paste0("figures/UDI/", growth.type,"/sp_uptake.depth_large_cor", goodness.fit, "_", growth.type, "_skIqAp.jpeg")), height = 5, width = 8, units ='in')

# ggsave(file.path(paste("figures/UDI/", growth.type,"/BCI/UDI_confidence/sp_UDI_R-sq>0.5 & N>5_by_UDmean", growth.type, "_#skIqAp.tiff")), height = 5, width = 8, units ='in', compression = "lzw")
# string(1, 6)
# "skIqAp"
# x <- factor(c("apple", "bear", "banana", "dear"))
# levels <- c(fruit = "apple", fruit = "banana")
# fct_recode(x, !!!levels)

g1 <- ggplot(ds %>% subset(!is.na(size) & !is.na(udi)) %>% transform(sp = reorder(sp, udi)) %>% droplevels(),
         aes(x = as.factor(sp), y = -udi)) +
  geom_errorbar(aes(ymax = -udi + udi.se, ymin = -udi - udi.se), colour = "black", width = 0.001, size = 0.5) +
  scale_y_continuous(trans="rev_sqrt", breaks = soil.depths) +
  theme(axis.text.x = element_text(size = 10, face = "plain", angle = 45, vjust = 1, hjust = 1)) +
  xlab("Species") + ylab(expression("Uptake Depth Index (m)")) +
  geom_point(colour = "darkgray", aes(shape = tlplevel, fill = sp, size = n)) +
  scale_shape_manual(values = c(21, 24)) + guides(fill = "none") +
  theme(axis.text.x = element_text(size = 4, face = "plain", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 20),
        strip.text = element_text(size = 20), legend.text = element_text(size = 16), legend.title = element_text(size = 20)) +
  facet_wrap(. ~ size) +  theme(legend.position = c(0.9, 0.9)) +
  scale_size(range = c(1, 3), "Sample Size", breaks = c(10, 100, 500))
ggsave(plot = g1, file.path(paste0("figures/UDI/", growth.type,"/sp_uptake.depth_all_sizes_cor", goodness.fit, "_", growth.type, growth.selection, "_", dbh.residuals, "_", "_skIqAp.jpeg")), height = 8, width = 12, units ='in') #compression = "lzw"

ggplot(ds, aes(x = size, y = udi)) +
  geom_boxplot(fill = "grey80") + xlab(expression("Size Class")) +
  facet_wrap(. ~ tlplevel) + scale_y_reverse() +
  ylab(expression("Uptake Depth (m)"))
ggsave(file.path(paste0("figures/UDI/", growth.type,"/boxplot_all_sizes_cor", goodness.fit, "_",growth.type, "_", growth.selection, "_", dbh.residuals, "_", "_skIqAp.jpeg")), height = 5, width = 7, units ='in') #compression = "lzw"

ggplot(ds, aes(x = size, y = udi)) +
  geom_violin(fill = "grey80", draw_quantiles = c(0.25, 0.5, 0.75)) + xlab(expression("Size Class")) +
  facet_wrap(. ~ tlplevel) + scale_y_reverse() +
  ylab(expression("Uptake Depth (m)"))
ggsave(file.path(paste0("figures/UDI/", growth.type,"/violinplot_all_sizes_cor", goodness.fit, "_",growth.type, "_skIqAp.jpeg")), height = 5, width = 7, units ='in') #compression = "lzw"

sp.classes <- with(ds, table(sp, size))
classes <- data.frame(sp = row.names(data.frame(rowSums(sp.classes))),
                      n = data.frame(rowSums(sp.classes))[,1])

ggplot(ds %>% subset(sp %in% classes$sp[classes$n >= 2]), aes(sp, as.factor(size), fill = -udi)) +
  geom_tile() + xlab(expression("Species")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab(expression("Size Class")) +
  scale_fill_viridis_c("Uptake\nDepth\nIndex\n(0-1)", option = "plasma") +
  ggtitle("Water Uptake Depth Index by Size Class")
ggsave(file.path(paste0("figures/UDI/", growth.type,"/sp_uptake.depth_all_sizes_", "n2", goodness.fit, "_", growth.type, "heatmap_skIqAp.jpeg")), height = 8.94, width = 8.94, units='in')

ggsave(file.path(paste0("figures/UDI/", growth.type,"/sp_uptake.depth_all_sizes_", "n3", goodness.fit, "_", growth.type, "heatmap_skIqAp.jpeg")), height = 8.94, width = 8.94, units='in')
