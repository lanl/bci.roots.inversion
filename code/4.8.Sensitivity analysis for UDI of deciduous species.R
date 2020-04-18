#-------------------------------
# Title: Sensitivity analysis for deciduous species:
#         Effect of dry season month exclusion on UDI
# Author : Rutuja Chitra-Tarak
# Original date: April 11, 2020
#-------------------------------

rm(list=ls())
gc()

pacman::p_load(tidyverse, readxl)

deci <- read_excel(file.path("data-raw/traits/nomenclature_R_20190524_Rready_Osvaldo Calderon & JoeWright_expert_opinion.xlsx")) %>%
  mutate(sp = tolower(sp6)) %>% select(sp, deciduous, tree)
tlp <- read.csv(file.path("data-raw/traits/HydraulicTraits_Kunert/tlp_sp_mean.csv")) %>%
  subset(!is.na(tlp))

# load interval and working.iter
load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
load(file.path("results/4.1GLUEsetup_part2_BCI.RData")) # has working.iter and growth and si matrix

soil.depths <- unique(info$root.param.long$depth)
intervals <- info$intervals
n.ensembles <- growth_by_si.info$n.ensembles
growth.type <- growth_by_si.info$growth.type
growth.selection <- growth_by_si.info$growth.selection
dbh.residuals <- growth_by_si.info$dbh.residuals
si.type <- growth_by_si.info$si.type
goodness.fit <- 0.3 # rsq0.3
n.best.fit <- 100
dryseason = "on"
iso.subset = "off"
root.selection = "exp"
##
level.folder <- c("splevel", "commlevel")

min.ll.list <- vector(mode = "list", length = length(level.folder))
names(min.ll.list) <- level.folder
drop.months.vec <- names(growth_by_si.info$si)

for (j in 1: length(level.folder)){
  if(level.folder[j] == "commlevel"){
    level.drop.months.vec <- c("Feb", "Mar")
  } else {
    level.drop.months.vec <- drop.months.vec
  }

  for (i in 1: length(level.drop.months.vec)) {
    file.extension.base1 <- paste0("drop.months", level.drop.months.vec[i], "_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_",
                                   growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals)
    rm(GLUE.negLL); rm(GLUE.rsq)
    load(file = paste0("results/", level.folder[j], "/GLUE.negLL_", file.extension.base1 , ".Rdata"), envir = parent.frame(), verbose = FALSE)
    load(file = paste0("results/", level.folder[j], "/GLUE.rsq_", file.extension.base1 , ".Rdata"), envir = parent.frame(), verbose = FALSE)
    load(file = paste0("results/", level.folder[j], "/GLUE.matches_", file.extension.base1 , ".Rdata"), envir = parent.frame(), verbose = FALSE)
    # retaining only those Rsq and negLL for which at least three data points pass through 95% CI of osbervations
    matrix.negLL <- matrix(unlist(GLUE.negLL), nrow = nrow(GLUE.negLL), ncol = ncol(GLUE.negLL))
    matrix.rsq <- matrix(unlist(GLUE.rsq), nrow = nrow(GLUE.negLL), ncol = ncol(GLUE.negLL))
    matrix.matches <- matrix(unlist(GLUE.matches), nrow = nrow(GLUE.negLL), ncol = ncol(GLUE.negLL))
    matrix.negLL[matrix.matches < 3] <- NA
    matrix.rsq[matrix.matches < 3] <- NA

    col.min.negLL <- as.numeric(apply(matrix.negLL, 1, function(x) {which(x == min(x, na.rm = TRUE))}))
    rsq.at.min.negLL <- vector()
    for (k in 1: length(col.min.negLL)) {rsq.at.min.negLL[k] <-  as.numeric(matrix.rsq[k, col.min.negLL[k]])}
    ### load udi as well
    file.extension.base4 <- paste0("drop.months", level.drop.months.vec[i], "_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type,
                                   "_", growth.selection, "_", dbh.residuals, "_", intervals,
                                   "_dryseason_", dryseason, "_iso.subset_", iso.subset, "_root.selection_", root.selection)
    # rm(ds.bestfit);
    load(file = paste0("results/", level.folder[j], "/ds.bestfit_", file.extension.base4, ".Rdata"))

    min.ll.list[[j]][[i]] <- data.frame(sp_size = row.names(GLUE.negLL),
                                   min.negLL = as.numeric(apply(matrix.negLL, 1, min, na.rm = TRUE)),
                                   max.negLL = as.numeric(apply(matrix.negLL, 1, max, na.rm = TRUE)),
                                   rsq.at.min.negLL = rsq.at.min.negLL,
                                   max.rsq = as.numeric(apply(matrix.rsq, 1, max, na.rm = TRUE)),
                                   drop.months = level.drop.months.vec[i],
                                   tlplevel = level.folder[j]) %>%
    left_join(ds.bestfit %>% select(sp_size, udi.best.rsq, udi.best.ll,
    udi.med.rsq, udi.med.ll, udi.upr.rsq, udi.lwr.rsq, udi.upr.ll, udi.lwr.ll), by = "sp_size")
  }
  names(min.ll.list[[j]]) <- level.drop.months.vec
  min.ll.list[[j]] <- do.call(rbind, min.ll.list[[j]])
}

min.ll <- do.call(rbind, min.ll.list)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
min.ll <- min.ll %>%
  separate(sp_size, into = c("sp", "size"), remove = FALSE) %>%
  #unite(col = drop.months_tlplevel, drop.months, tlplevel, remove = FALSE) %>%
  group_by(sp_size) %>%
  mutate(min.negLL.norm = range01(min.negLL),
         max.rsq.norm = range01(max.rsq)) %>% ungroup(sp_size) %>%
  group_by(sp_size, tlplevel) %>%
  mutate(min.negLL.norm.tlp = range01(min.negLL),
         max.rsq.norm.tlp = range01(max.rsq)) %>% ungroup(sp_size, tlplevel) %>%
  mutate(drop.months = factor(drop.months, levels =
                                   c("None", "Jan", "Feb", "Mar", "JanFeb", "FebMar", "JanFebMar"))) %>%
  left_join(deci, by = "sp") %>%
  left_join(tlp, by = "sp") %>%
  mutate(tlp.avail = case_when(is.na(tlp) ~ "NoTLP", !is.na(tlp) ~ "TLP"))
min.ll <- min.ll %>%
  arrange(desc(deciduous), desc(sp)) %>%
  unite(col = deciduous_sp, deciduous, sp, remove = FALSE) %>%
  mutate(deciduous_sp = fct_relevel(deciduous_sp, function(x) {sort(x, decreasing = TRUE)})) %>%
  subset(!deciduous %in% c("NA", "?")) %>%
  unite(col = tlp_deci, tlp.avail, deciduous, remove = FALSE) %>%
  unite(col = tlp_sp, tlp_deci, sp, remove = FALSE) %>%
  mutate(tlp_sp = fct_relevel(tlp_sp, function(x) {sort(x, decreasing = TRUE)}))

# ## Minimum likelihood sensitivity-------
# heat.ll <- ggplot(min.ll %>% subset(size == "large") %>% droplevels(),
#                      aes(y = deciduous_sp, x = as.factor(drop.months))) +
#   ylab("Species") + xlab("Months Dropped from Btran Calculation") +
#   facet_wrap(. ~ tlplevel) +
#   scale_fill_viridis_c("Normalised\nMinimum\nNegative\nLikelihood", direction = -1, option = "plasma") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
#
# heat.ll.abs <- heat.ll + geom_tile(aes(fill = min.negLL)) +
#   ggtitle("Absolute minimum negative Log Likelihood")
# ggsave("min.negLL_within_sp_across_drop.months_absolute.jpeg", plot = heat.ll.abs, path =
#          file.path("figures/UDI_confidence/sensitivity_analysis"), height = 8.94, width = 8.94, units='in')
#
# heat.ll.abs.max <- heat.ll + geom_tile(aes(fill = max.negLL)) +
#   ggtitle("Absolute maximum negative Log Likelihood")
# ggsave("max.negLL_within_sp_across_drop.months_absolute.jpeg", plot = heat.ll.abs.max, path =
#          file.path("figures/UDI_confidence/sensitivity_analysis"), height = 8.94, width = 8.94, units='in')
#
# heat.ll.one <- heat.ll + geom_tile(aes(fill = min.negLL.norm)) +
#   ggtitle("minimum negative Log Likelihood normalised\ntogether for both tlplevels")
# ggsave("norm.negLL_within_sp_across_drop.months_together_norm_tlplevel.jpeg", plot = heat.ll.one, path =
#          file.path("figures/UDI_confidence/sensitivity_analysis"), height = 8.94, width = 8.94, units='in')
#
# heat.ll.two <- heat.ll + geom_tile(aes(fill = min.negLL.norm.tlp)) +
#   ggtitle("minimum negative Log Likelihood normalised\n separately for the two tlplevels")
# ggsave("norm.negLL_within_sp_across_drop.months_separately_norm_tlplevel.jpeg", plot = heat.ll.two, path =
#          file.path("figures/UDI_confidence/sensitivity_analysis"), height = 8.94, width = 8.94, units='in')
# ## Rsq-------
# heat.rsq <- ggplot(min.ll %>% subset(size == "large") %>% droplevels(),
#                   aes(y = deciduous_sp, x = as.factor(drop.months))) +
#   ylab("Species") + xlab("Months Dropped from Btran Calculation") +
#   facet_wrap(. ~ tlplevel) +
#   scale_fill_viridis_c("R-squared", direction = -1, option = "plasma") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# heat.rsq.one <- heat.rsq + geom_tile(aes(fill = max.rsq.norm)) +
#   ggtitle("maximum Rsq normalised together\nfor both tlplevels")
# ggsave("rsq_within_sp_across_drop.months_together_norm_tlplevel.jpeg", plot = heat.rsq.one, path =
#          file.path("figures/UDI_confidence/sensitivity_analysis"), height = 8.94, width = 8.94, units='in')
# heat.rsq.two <- heat.rsq + geom_tile(aes(fill = max.rsq.norm.tlp)) +
#   ggtitle("maximum Rsq normalised separately\nfor the tlplevels")
# ggsave("rsq_within_sp_across_drop.months_separately_norm_tlplevel.jpeg", plot = heat.rsq.two, path =
#          file.path("figures/UDI_confidence/sensitivity_analysis"), height = 8.94, width = 8.94, units='in')
# heat.rsq.abs <- heat.rsq + geom_tile(aes(fill = max.rsq)) +
#   ggtitle("Absolute maxRsq")
# ggsave("rsq_within_sp_across_drop.months_absolute_rsq.jpeg", plot = heat.rsq.abs, path =
#          file.path("figures/UDI_confidence/sensitivity_analysis"), height = 8.94, width = 8.94, units='in')
#
# heat.rsq.abs <- heat.rsq + geom_tile(aes(fill = rsq.at.min.negLL)) +
#   ggtitle("Absolute Rsq at min negLL")
# ggsave("rsq_within_sp_across_drop.months_absolute_rsq at min negLL.jpeg", plot = heat.rsq.abs, path =
#          file.path("figures/UDI_confidence/sensitivity_analysis"), height = 8.94, width = 8.94, units='in')

## UDI sensitivity to dropped months------
udi.ll <- ggplot(min.ll %>% subset(size == "large") %>% droplevels(),
                  aes(y = deciduous_sp, x = as.factor(drop.months))) +
  ylab("Species") + xlab("Months Dropped from Btran Calculation") +
  facet_wrap(. ~ tlplevel) +
  scale_fill_viridis_c("Uptake Depth Index", direction = -1, option = "plasma") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

best.udi.ll <- udi.ll + geom_tile(aes(fill = udi.best.ll)) +
  ggtitle("UDI for model with minimum negative Log Likelihood")
ggsave("best.UDI_within_sp_across_drop.months_min.negLL.jpeg", plot = best.udi.ll, path =
         file.path("figures/UDI_confidence/sensitivity_analysis"), height = 8.94, width = 8.94, units='in')

best.udi.rsq <- udi.ll + geom_tile(aes(fill = udi.best.rsq)) +
  ggtitle("UDI for model with maximum R-squared")
ggsave("best.UDI_within_sp_across_drop.months_max.rsq.jpeg", plot = best.udi.rsq, path =
         file.path("figures/UDI_confidence/sensitivity_analysis"), height = 8.94, width = 8.94, units='in')

top.udi.ll <- udi.ll + geom_tile(aes(fill = min.negLL)) +
  ggtitle("UDI for model with minimum negative Log Likelihood")
ggsave("best10.UDI_within_sp_across_drop.months_min.negLL.jpeg", plot = best.udi.ll, path =
         file.path("figures/UDI_confidence/sensitivity_analysis"), height = 8.94, width = 8.94, units='in')

top.udi.rsq <- udi.ll + geom_tile(aes(fill = min.negLL)) +
  ggtitle("UDI for model with maximum R-squared")
ggsave("best10.UDI_within_sp_across_drop.months_max.rsq.jpeg", plot = best.udi.ll, path =
         file.path("figures/UDI_confidence/sensitivity_analysis"), height = 8.94, width = 8.94, units='in')
