
#-----------------------------------------------------
# Title: Sensitivity analysis for deciduous species:
#         Effect of dry season month exclusion on UDI
# Author : Rutuja Chitra-Tarak
# Original date: April 11, 2020
#-----------------------------------------------------

rm(list=ls())
gc()

pacman::p_load(tidyverse, readxl)
sensitivity.analysis <- function(goodness.fit = goodness.fit,
                        dryseason = dryseason,
                        root.selection =  root.selection,
                        iso.subset = iso.subset, drop.months = drop.months){

deci <- read_excel(file.path("data-raw/traits/nomenclature_R_20190524_Rready_Osvaldo Calderon & JoeWright_expert_opinion.xlsx")) %>%
  mutate(sp = tolower(sp6)) %>% select(sp, deciduous, tree)
tlp <- read.csv(file.path("data-raw/traits/HydraulicTraits_Kunert/tlp_sp_mean.csv")) %>%
  subset(!is.na(tlp))

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
goodness.fit <- 0.3
##
level.folder <- c("splevel", "commlevel")

max.ll.list <- vector(mode = "list", length = length(level.folder))
names(max.ll.list) <- level.folder
load(file.path("results/GLUEsetup_part1.2_BCI.RData"))
drop.months.vec <- names(info.2$sp.si)

for (j in 1: length(level.folder)){
  for (i in 1: length(drop.months.vec)) {
    file.extension.base1 <- paste0("drop.months", drop.months.vec[i], "_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_",
                                   growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals)
    if( !(j == 1 & i == 1)) {rm(GLUE.negLL); rm(GLUE.rsq)}
    load(file = paste0("results/", level.folder[j], "/GLUE.negLL_", file.extension.base1 , ".Rdata"), envir = parent.frame(), verbose = FALSE)
    load(file = paste0("results/", level.folder[j], "/GLUE.rsq_", file.extension.base1 , ".Rdata"), envir = parent.frame(), verbose = FALSE)
    load(file = paste0("results/", level.folder[j], "/GLUE.matches_", file.extension.base1 , ".Rdata"), envir = parent.frame(), verbose = FALSE)
    # retaining only those Rsq and negLL for which at least three data points pass through 95% CI of osbervations
    matrix.likelihood <- matrix(exp(-unlist(GLUE.negLL)), nrow = nrow(GLUE.negLL), ncol = ncol(GLUE.negLL))
    matrix.rsq <- matrix(unlist(GLUE.rsq), nrow = nrow(GLUE.negLL), ncol = ncol(GLUE.negLL))
    matrix.matches <- matrix(unlist(GLUE.matches), nrow = nrow(GLUE.negLL), ncol = ncol(GLUE.negLL))
    matrix.likelihood[matrix.matches < 3] <- NA
    matrix.rsq[matrix.matches < 3] <- NA

    col.max.likelihood <- as.numeric(apply(matrix.likelihood, 1, function(x) {which(x == max(x, na.rm = TRUE))}))
    rsq.at.max.likelihood <- vector()
    for (k in 1: length(col.max.likelihood)) {rsq.at.max.likelihood[k] <-  as.numeric(matrix.rsq[k, col.max.likelihood[k]])}
    ### load udi as well
    file.extension.base4 <- paste0("drop.months", drop.months.vec[i], "_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type,
                                   "_", growth.selection, "_", dbh.residuals, "_", intervals,
                                   "_dryseason_", dryseason, "_iso.subset_", iso.subset, "_root.selection_", root.selection)
    # rm(ds.bestfit);
    load(file = paste0("results/", level.folder[j], "/ds.bestfit_", file.extension.base4, ".Rdata"))

    max.ll.list[[j]][[i]] <- data.frame(sp_size = row.names(GLUE.negLL),
                                   max.likelihood = as.numeric(apply(matrix.likelihood, 1, max, na.rm = TRUE)),
                                   rsq.at.max.likelihood = rsq.at.max.likelihood,
                                   max.rsq = as.numeric(apply(matrix.rsq, 1, max, na.rm = TRUE)),
                                   drop.months = drop.months.vec[i],
                                   tlplevel = level.folder[j]) %>%
    left_join(ds.bestfit %>% select(sp_size, udi.best.rsq, udi.best.ll,
    udi.med.rsq, udi.med.ll, udi.upr.rsq, udi.lwr.rsq, udi.upr.ll, udi.lwr.ll), by = "sp_size")
  }
  names(max.ll.list[[j]]) <- drop.months.vec
  max.ll.list[[j]] <- do.call(rbind, max.ll.list[[j]])
}

max.ll <- do.call(rbind, max.ll.list)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
max.ll <- max.ll %>%
  separate(sp_size, into = c("sp", "size"), remove = FALSE) %>%
  #unite(col = drop.months_tlplevel, drop.months, tlplevel, remove = FALSE) %>%
  group_by(sp_size) %>%
  mutate(max.likelihood.norm = range01(max.likelihood),
         max.rsq.norm = range01(max.rsq)) %>% ungroup(sp_size) %>%
  group_by(sp_size, tlplevel) %>%
  mutate(max.likelihood.norm.tlp = range01(max.likelihood),
         max.rsq.norm.tlp = range01(max.rsq)) %>% ungroup(sp_size, tlplevel) %>%
  mutate(drop.months = factor(drop.months, levels =
                                   c("None", "Jan", "Feb", "Mar", "JanFeb", "FebMar", "JanFebMar"))) %>%
  left_join(deci, by = "sp") %>%
  left_join(tlp, by = "sp") %>%
  mutate(tlp.avail = case_when(is.na(tlp) ~ "NoTLP", !is.na(tlp) ~ "TLP"))
max.ll <- max.ll %>%
  arrange(desc(deciduous), desc(sp)) %>%
  unite(col = deciduous_sp, deciduous, sp, remove = FALSE) %>%
  mutate(deciduous_sp = fct_relevel(deciduous_sp, function(x) {sort(x, decreasing = TRUE)})) %>%
  subset(!deciduous %in% c("NA", "?")) %>%
  unite(col = tlp_deci, tlp.avail, deciduous, remove = FALSE) %>%
  unite(col = tlp_sp, tlp_deci, sp, remove = FALSE) %>%
  mutate(tlp_sp = fct_relevel(tlp_sp, function(x) {sort(x, decreasing = TRUE)}))

file.extension.base5 <- paste0("_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type,
                               "_", growth.selection, "_", dbh.residuals, "_", intervals,
                               "_dryseason_", dryseason, "_iso.subset_", iso.subset, "_root.selection_", root.selection)

save(max.ll, file = paste0("results/max.ll_", file.extension.base5, ".Rdata"))

### For figures, creating folders if not present

load(file.path("results/GLUEsetup_part2_BCI.RData")) # has working.iter and growth and si matrix

growth.type <- growth_by_si.info$growth.type
growth.selection <- growth_by_si.info$growth.selection
dbh.residuals <- growth_by_si.info$dbh.residuals

if(!dir.exists(file.path("figures"))) {dir.create(file.path("figures"))}
if(!dir.exists(file.path("figures", "UDI_confidence"))) {dir.create(file.path("figures","UDI_confidence"))}
if(!dir.exists(file.path("figures", "UDI_confidence", growth.type))) {dir.create(file.path("figures","UDI_confidence", growth.type))}
if(!dir.exists(file.path("figures", "UDI_confidence", growth.type, growth.selection))) {
  dir.create(file.path("figures","UDI_confidence", growth.type, growth.selection))}
if(!dir.exists(file.path("figures", "UDI_confidence", growth.type, growth.selection,  paste0("dbh.residuals_", dbh.residuals)))) {
  dir.create(file.path("figures","UDI_confidence", growth.type, growth.selection,  paste0("dbh.residuals_", dbh.residuals)))}
if(!dir.exists(file.path("figures", "UDI_confidence", growth.type, growth.selection,  paste0("dbh.residuals_", dbh.residuals), paste0("dryssn_", dryseason, "_root.selection_", root.selection)))) {
  dir.create(file.path("figures","UDI_confidence", growth.type, growth.selection,  paste0("dbh.residuals_", dbh.residuals), paste0("dryssn_", dryseason, "_root.selection_", root.selection)))}
if(!dir.exists(file.path("figures", "UDI_confidence", growth.type, growth.selection,  paste0("dbh.residuals_", dbh.residuals), paste0("dryssn_", dryseason, "_root.selection_", root.selection), "sensitivity_analysis"))) {
  dir.create(file.path("figures","UDI_confidence", growth.type, growth.selection,  paste0("dbh.residuals_", dbh.residuals), paste0("dryssn_", dryseason, "_root.selection_", root.selection), "sensitivity_analysis"))}

file.path.ll <- file.path("figures","UDI_confidence", growth.type, growth.selection,
                          paste0("dbh.residuals_", dbh.residuals), paste0("dryssn_", dryseason, "_root.selection_", root.selection), "sensitivity_analysis")

## maximum likelihood sensitivity-------
heat.ll <- ggplot(max.ll %>% subset(size == "large") %>% droplevels(),
                     aes(y = deciduous_sp, x = as.factor(drop.months))) +
  ylab("Species") + xlab("Months Dropped from Btran Calculation") +
  facet_wrap(. ~ tlplevel) +
  scale_fill_viridis_c("Normalised\nMaximum\nLikelihood", direction = -1, option = "plasma") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

heat.ll.abs <- heat.ll + geom_tile(aes(fill = max.likelihood)) +
  ggtitle("Absolute Maximum Likelihood")
scale_fill_viridis_c("Absolute\nMaximum\nLikelihood", direction = -1, option = "plasma")
ggsave("max.likelihood_within_sp_across_drop.months_absolute.jpeg", plot = heat.ll.abs, path =
         file.path(file.path.ll), height = 8.94, width = 8.94, units='in')

heat.ll.one <- heat.ll + geom_tile(aes(fill = max.likelihood.norm)) +
  ggtitle("Maximum likelihood normalised\ntogether for both tlplevels")
ggsave("norm.likelihood_within_sp_across_drop.months_together_norm_tlplevel.jpeg", plot = heat.ll.one, path =
         file.path(file.path.ll), height = 8.94, width = 8.94, units='in')

heat.ll.two <- heat.ll + geom_tile(aes(fill = max.likelihood.norm.tlp)) +
  ggtitle("Maximum likelihood normalised\nseparately for the two tlplevels")
ggsave("norm.likelihood_within_sp_across_drop.months_separately_norm_tlplevel.jpeg", plot = heat.ll.two, path =
         file.path(file.path.ll), height = 8.94, width = 8.94, units='in')
## Rsq-------
heat.rsq <- ggplot(max.ll %>% subset(size == "large") %>% droplevels(),
                  aes(y = deciduous_sp, x = as.factor(drop.months))) +
  ylab("Species") + xlab("Months Dropped from Btran Calculation") +
  facet_wrap(. ~ tlplevel) +
  scale_fill_viridis_c("R-squared", direction = -1, option = "plasma") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
heat.rsq.one <- heat.rsq + geom_tile(aes(fill = max.rsq.norm)) +
  ggtitle("Maximum Rsq normalised together\nfor both tlplevels")
ggsave("rsq_within_sp_across_drop.months_together_norm_tlplevel.jpeg", plot = heat.rsq.one, path =
         file.path(file.path.ll), height = 8.94, width = 8.94, units='in')
heat.rsq.two <- heat.rsq + geom_tile(aes(fill = max.rsq.norm.tlp)) +
  ggtitle("Maximum Rsq normalised separately\nfor the tlplevels")
ggsave("rsq_within_sp_across_drop.months_separately_norm_tlplevel.jpeg", plot = heat.rsq.two, path =
         file.path(file.path.ll), height = 8.94, width = 8.94, units='in')
heat.rsq.abs <- heat.rsq + geom_tile(aes(fill = max.rsq)) +
  ggtitle("Absolute Maximum Rsq")
ggsave("rsq_within_sp_across_drop.months_absolute_rsq.jpeg", plot = heat.rsq.abs, path =
         file.path(file.path.ll), height = 8.94, width = 8.94, units='in')

heat.rsq.abs <- heat.rsq + geom_tile(aes(fill = rsq.at.max.likelihood)) +
  ggtitle("Absolute Rsq at Maximum likelihood")
ggsave("rsq_within_sp_across_drop.months_absolute_rsq at max likelihood.jpeg", plot = heat.rsq.abs, path =
         file.path(file.path.ll), height = 8.94, width = 8.94, units='in')

## UDI sensitivity to dropped months------
heat.udi <- ggplot(max.ll %>% subset(size == "large") %>% droplevels(),
                  aes(y = deciduous_sp, x = as.factor(drop.months))) +
  ylab("Species") + xlab("Months Dropped from Btran Calculation") +
  facet_wrap(. ~ tlplevel) +
  scale_fill_viridis_c("Uptake Depth Index", direction = -1, option = "viridis") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

best.udi.ll <- heat.udi + geom_tile(aes(fill = udi.best.ll)) +
  ggtitle("UDI for model with maximum likelihood")
ggsave("best.UDI_within_sp_across_drop.months_max.likelihood.jpeg", plot = best.udi.ll, path =
         file.path(file.path.ll), height = 8.94, width = 8.94, units='in')

best.udi.rsq <- heat.udi + geom_tile(aes(fill = udi.best.rsq)) +
  ggtitle("UDI for model with maximum R-squared")
ggsave("best.UDI_within_sp_across_drop.months_max.rsq.jpeg", plot = best.udi.rsq, path =
         file.path(file.path.ll), height = 8.94, width = 8.94, units='in')

top.udi.ll <- heat.udi + geom_tile(aes(fill = udi.med.ll)) +
  ggtitle("Median UDI for top ranking models\nbased on maximum likelihood")
ggsave("med.top.UDI_within_sp_across_drop.months_max.likelihood.jpeg", plot = top.udi.ll, path =
         file.path(file.path.ll), height = 8.94, width = 8.94, units='in')

top.udi.ll.upr <- heat.udi + geom_tile(aes(fill = udi.upr.ll)) +
  ggtitle("97.5 quantile UDI for top ranking models\nbased on maximum likelihood")
ggsave("upr.top.UDI_within_sp_across_drop.months_max.likelihood.jpeg", plot = top.udi.ll.upr, path =
         file.path(file.path.ll), height = 8.94, width = 8.94, units='in')

top.udi.rsq <- heat.udi + geom_tile(aes(fill = udi.med.rsq)) +
  ggtitle("Median UDI for top ranking models\nbased on maximum R-squared")
ggsave("med.top.UDI_within_sp_across_drop.months_max.rsq.jpeg", plot = top.udi.rsq, path =
         file.path(file.path.ll), height = 8.94, width = 8.94, units='in')
}

