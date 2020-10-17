
#-----------------------------------------------------
# Title: Sensitivity analysis for deciduous species:
#         Effect of dry season month exclusion on UDI
# Author : Rutuja Chitra-Tarak
# Original date: April 11, 2020
#-----------------------------------------------------

rm(list=ls())
gc()

pacman::p_load(tidyverse, readxl, gridExtra, magrittr, multipanelfigure)
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
    inverse = function(x) x^2)
}
reverselog_trans <- function(base = exp(1)) {
  scales::trans_new(name = paste0("reverselog-", format(base)),
                    log_breaks(base = base),
                    domain = c(1e-100, Inf),
                    transform = function(x) -log(x, base),
                    inverse = function(x) base^(-x))
}
sensitivity.analysis <- function(goodness.fit = goodness.fit,
                        dryseason.vec = dryseason.vec,
                        root.selection = root.selection,
                        iso.subset = iso.subset) {

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
##
level.folder <- c("splevel", "commlevel")

max.ll.list <- vector(mode = "list", length = length(dryseason.vec))
names(max.ll.list) <- dryseason.vec

drop.months.vec <- names(growth_by_si.info$si.param.rel)

for (m in 1: length(dryseason.vec)) {
  max.ll.list[[m]] <- vector(mode = "list", length = length(level.folder))
  names(max.ll.list[[m]]) <- level.folder
  for (j in 1: length(level.folder)) {
    for (i in 1: length(drop.months.vec)) {
      file.extension.base1 <- paste0("drop.months", drop.months.vec[i], "_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_",
                                     growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals)
      if(!(j == 1 & i == 1)) {rm(GLUE.negLL); rm(GLUE.rsq)}
      load(file = paste0("results/", level.folder[j], "/GLUE.negLL_", file.extension.base1 , ".Rdata"), envir = parent.frame(), verbose = FALSE)
      load(file = paste0("results/", level.folder[j], "/GLUE.rsq_", file.extension.base1 , ".Rdata"), envir = parent.frame(), verbose = FALSE)
      load(file = paste0("results/", level.folder[j], "/GLUE.matches_", file.extension.base1 , ".Rdata"), envir = parent.frame(), verbose = FALSE)
      # retaining only those Rsq and negLL for which at least three data points pass through 95% CI of osbervations
      matrix.likelihood <- matrix(exp(-unlist(GLUE.negLL)), nrow = nrow(GLUE.negLL), ncol = ncol(GLUE.negLL))
      matrix.rsq <- matrix(unlist(GLUE.rsq), nrow = nrow(GLUE.negLL), ncol = ncol(GLUE.negLL))
      matrix.matches <- matrix(unlist(GLUE.matches), nrow = nrow(GLUE.negLL), ncol = ncol(GLUE.negLL))
      matrix.likelihood[matrix.matches < 3] <- NA
      matrix.rsq[matrix.matches < 3] <- NA
      row.names(matrix.matches) <- row.names(matrix.likelihood) <- row.names( matrix.rsq) <- row.names(GLUE.negLL)
      ### load udi as well
      file.extension.base4 <- paste0("drop.months", drop.months.vec[i], "_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type,
                                     "_", growth.selection, "_", dbh.residuals, "_", intervals,
                                     "_dryseason_", dryseason.vec[m], "_iso.subset_", iso.subset, "_root.selection_", root.selection)
      # rm(ds.bestfit);
      load(file = paste0("results/", level.folder[j], "/ds.bestfit_", file.extension.base4, ".Rdata"))
      matrix.likelihood <- matrix.likelihood[row.names(matrix.likelihood) %in% ds.bestfit$sp_size,]
      matrix.rsq <- matrix.rsq[row.names(matrix.rsq) %in% ds.bestfit$sp_size,]
      matrix.matches <- matrix.matches[row.names(matrix.matches) %in% ds.bestfit$sp_size,]

      col.max.likelihood <- as.numeric(apply(matrix.likelihood, 1, function(x) {which(x == max(x, na.rm = TRUE))}))
      rsq.at.max.likelihood <- vector()
      for (k in 1: length(col.max.likelihood)) {
        rsq.at.max.likelihood[k] <-  as.numeric(matrix.rsq[k, col.max.likelihood[k]])
        }
      max.ll.list[[m]][[j]][[i]] <- data.frame(sp_size = ds.bestfit$sp_size,
                                          max.likelihood = as.numeric(apply(matrix.likelihood, 1, max, na.rm = TRUE)),
                                          rsq.at.max.likelihood = rsq.at.max.likelihood,
                                          max.rsq = as.numeric(apply(matrix.rsq, 1, max, na.rm = TRUE)),
                                          drop.months = drop.months.vec[i],
                                          tlplevel = level.folder[j],
                                          dryseason = dryseason.vec[m]) %>%
        left_join(ds.bestfit %>% select(sp_size, udi.best.rsq, udi.best.ll,
                                        udi.med.rsq, udi.med.ll, udi.upr.rsq, udi.lwr.rsq, udi.upr.ll, udi.lwr.ll), by = "sp_size")
    }
    names(max.ll.list[[m]][[j]]) <- drop.months.vec
    max.ll.list[[m]][[j]] <- do.call(rbind, max.ll.list[[m]][[j]])
  }
  max.ll.list[[m]] <- do.call(rbind, max.ll.list[[m]])
}

max.ll <- do.call(rbind, max.ll.list)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
max.ll <- max.ll %>%
  separate(sp_size, into = c("sp", "size"), remove = FALSE) %>%
  #unite(col = drop.months_tlplevel, drop.months, tlplevel, remove = FALSE) %>%
  ## difference between early and late dryseason
  mutate(dryseason = recode_factor(dryseason, on = "Mar-Apr", jan = "Jan")) %>%
  mutate(time = recode(dryseason, `Mar-Apr` = 2, `Jan` = 1)) %>%
  arrange(sp_size, drop.months, tlplevel, time) %>%
  group_by(sp_size, drop.months, tlplevel) %>%
  mutate(change.udi.best.rsq = udi.best.rsq - lag(udi.best.rsq, default = NA),
         change.udi.best.ll = udi.best.ll - lag(udi.best.ll, default = NA),
         change.udi.med.rsq = udi.med.rsq - lag(udi.med.rsq, default = NA),
         change.udi.med.ll = udi.med.ll - lag(udi.med.ll, default = NA)) %>%
  group_by(dryseason, sp_size) %>%
  mutate(max.likelihood.norm = range01(max.likelihood),
         max.rsq.norm = range01(max.rsq)) %>% ungroup(sp_size) %>%
  group_by(dryseason, sp_size, tlplevel) %>%
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
                               "_dryseason_", paste0(dryseason.vec, collapse = "-"), "_iso.subset_", iso.subset, "_root.selection_", root.selection)

# load(file = paste0("results/max.ll", file.extension.base5, ".Rdata"))
load(file = "results/all_isotopic_record.Rdata")

iso.all.sort <- iso.all %>% group_by(sp) %>%
  summarise(mean.iso = mean(Xylem_sap_deltaD_permil, na.rm = TRUE)) %>%
  mutate(sp = as.character(sp)) %>%
  subset(!is.na(mean.iso)) %>% droplevels()
plot(mean.iso ~ sp, data = iso.all.sort)

max.ll <- max.ll %>%
  mutate(sp = as.character(sp)) %>%
  left_join(iso.all.sort, by = "sp") %>%
  transform(sp = reorder(sp, mean.iso),
            deciduous_sp = reorder(deciduous_sp, mean.iso))
# mutate(sp = as.factor(sp),
#        deciduous_sp = as.factor(deciduous_sp)) %>%
plot(mean.iso ~ sp, data = max.ll)
save(max.ll, file = paste0("results/max.ll", file.extension.base5, ".Rdata"))

### For figures, creating folders if not present

if(!dir.exists(file.path("figures"))) {dir.create(file.path("figures"))}
if(!dir.exists(file.path("figures", "UDI_confidence"))) {dir.create(file.path("figures","UDI_confidence"))}
if(!dir.exists(file.path("figures", "UDI_confidence", growth.type))) {dir.create(file.path("figures","UDI_confidence", growth.type))}
if(!dir.exists(file.path("figures", "UDI_confidence", growth.type, growth.selection))) {
  dir.create(file.path("figures","UDI_confidence", growth.type, growth.selection))}
if(!dir.exists(file.path("figures", "UDI_confidence", growth.type, growth.selection,  paste0("dbh.residuals_", dbh.residuals)))) {
  dir.create(file.path("figures","UDI_confidence", growth.type, growth.selection,  paste0("dbh.residuals_", dbh.residuals)))}
if(!dir.exists(file.path("figures", "UDI_confidence", growth.type, growth.selection,  paste0("dbh.residuals_", dbh.residuals), paste0("dryssn_", paste0(dryseason.vec, collapse = "-"), "_root.selection_", root.selection)))) {
  dir.create(file.path("figures","UDI_confidence", growth.type, growth.selection,  paste0("dbh.residuals_", dbh.residuals), paste0("dryssn_", paste0(dryseason.vec, collapse = "-"), "_root.selection_", root.selection)))}
if(!dir.exists(file.path("figures", "UDI_confidence", growth.type, growth.selection,  paste0("dbh.residuals_", dbh.residuals), paste0("dryssn_", paste0(dryseason.vec, collapse = "-"), "_root.selection_", root.selection), "sensitivity_analysis"))) {
  dir.create(file.path("figures","UDI_confidence", growth.type, growth.selection,  paste0("dbh.residuals_", dbh.residuals), paste0("dryssn_", paste0(dryseason.vec, collapse = "-"), "_root.selection_", root.selection), "sensitivity_analysis"))}

file.path.ll <- file.path("figures","UDI_confidence", growth.type, growth.selection,
                          paste0("dbh.residuals_", dbh.residuals), paste0("dryssn_", paste0(dryseason.vec, collapse = "-"), "_root.selection_", root.selection), "sensitivity_analysis")

## maximum likelihood sensitivity-------
heat.ll <- ggplot(max.ll %>% subset(size == "large" & dryseason == "Mar-Apr") %>% droplevels(),
                     aes(y = deciduous_sp, x = as.factor(drop.months))) +
  ylab("Species") + xlab("Months Dropped from Btran Calculation") +
  facet_wrap(. ~ tlplevel) +
  scale_fill_viridis_c("Normalised\nMaximum\nLikelihood", direction = -1, option = "plasma") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

heat.ll.abs <- heat.ll + geom_tile(aes(fill = max.likelihood)) +
  ggtitle("Absolute Maximum Likelihood") +
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
heat.rsq <- ggplot(max.ll %>% subset(size == "large" & dryseason == "Mar-Apr") %>% droplevels(),
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
  facet_wrap(dryseason ~ tlplevel) +
  scale_fill_viridis_c("Uptake\nDepth\nIndex", direction = -1, option = "viridis") +
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

## Change in UDI as dry season progresses------
change.heat.udi <- ggplot(max.ll %>% subset(size == "large") %>% droplevels(),
                   aes(y = deciduous_sp, x = as.factor(drop.months))) +
  ylab("Species") + xlab("Months Dropped from Btran Calculation") +
  facet_wrap(. ~ tlplevel) +
  scale_fill_viridis_c("Uptake\nDepth\nIndex", direction = -1, option = "viridis") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

change.best.udi.ll <- change.heat.udi + geom_tile(aes(fill = change.udi.best.ll)) +
  ggtitle("Jan to Mar change in UDI for model with maximum likelihood")
ggsave("change.best.UDI_within_sp_across_drop.months_max.likelihood.jpeg", plot = change.best.udi.ll, path =
         file.path(file.path.ll), height = 8.94, width = 8.94, units='in')

change.best.udi.rsq <- change.heat.udi + geom_tile(aes(fill = change.udi.best.rsq)) +
  ggtitle("Jan to Mar change in UDI for model with maximum R-sq")
ggsave("change.best.UDI_within_sp_across_drop.months_max.rsq.jpeg", plot = change.best.udi.rsq, path =
         file.path(file.path.ll), height = 8.94, width = 8.94, units='in')

change.top.udi.ll <- change.heat.udi + geom_tile(aes(fill = change.udi.med.ll)) +
  ggtitle("Jan to Mar change in Median UDI for top ranking models\nbased on maximum likelihood")
ggsave("change.med.top.UDI_within_sp_across_drop.months_max.ll.jpeg", plot = change.top.udi.ll, path =
         file.path(file.path.ll), height = 8.94, width = 8.94, units='in')

change.top.udi.rsq <- change.heat.udi + geom_tile(aes(fill = change.udi.med.rsq)) +
  ggtitle("Jan to Mar change in Median UDI for top ranking models\nbased on maximum R-squared")
ggsave("change.med.top.UDI_within_sp_across_drop.months_max.rsq.jpeg", plot = change.top.udi.rsq, path =
         file.path(file.path.ll), height = 8.94, width = 8.94, units='in')
### UDI versus isotopic records
iso.udi <- max.ll %>% left_join(iso.all, by = "sp")
iso.udi.sub <- iso.udi %>% subset(deciduous == "E" & drop.months == "None") %>%
  rbind(iso.udi %>% subset(deciduous == "DF" & drop.months == "Jan")) %>%
  rbind(iso.udi %>% subset(deciduous == "DB" & drop.months == "JanFeb")) %>%
  rbind(iso.udi %>% subset(deciduous == "DO" & drop.months == "JanFebMar"))

ggplot(iso.udi.sub %>% subset(tlplevel == "splevel" & dryseason == "Mar-Apr"),
       aes(x = mean.iso, y = udi.best.rsq, color = deciduous_sp)) + # drop.months
  # geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + se, xmin = Xylem_sap_deltaD_permil - se),
  #                size = 0.5) +
  # geom_errorbar(aes(ymax = udi.med.rsq + udi.upr.rsq, ymin = udi.med.rsq - udi.lwr.rsq),
  #               size = 0.5, width = 0.2) +
  geom_jitter(aes(shape = deciduous), size = 3, width = 0.25) +

  # scale_color_viridis_c("TLP [MPa]", option = "plasma", direction = -1) +
  # geom_text(aes(x =  Xylem_sap_deltaD_permil + 2.5, y = udi.best.rsq + sqrt(udi.med.rsq*diff(range(iso.udi.sub$udi.best.rsq, na.rm = TRUE)))/20, label = sp),
  #           size = 4) +
  scale_y_continuous(trans="rev_sqrt", breaks = signif(soil.depths, 2))
  # scale_y_continuous(trans=reverselog_trans(10), breaks = signif(soil.depths, 2))

g1 <- ggplot(iso.udi %>% subset(tlplevel == "splevel" & dryseason == "Mar-Apr" &
                            max.likelihood.norm.tlp == 1 & size == "large"), # &
       aes(x = mean.iso, color = deciduous_sp)) + # drop.months
  scale_y_continuous(trans="rev_sqrt", breaks = signif(soil.depths, 2),
                     limits = c(range(iso.udi$udi.best.ll, na.rm = TRUE)[2], range(iso.udi$udi.best.ll, na.rm = TRUE)[1])) +
  theme(legend.position = "None") +
  xlim(c(range(iso.udi$mean.iso, na.rm = TRUE)[1] - 5, range(iso.udi$mean.iso, na.rm = TRUE)[2] + 5))

g1.all <- g1 %+% subset(iso.udi, tlplevel == "commlevel" & dryseason == "Mar-Apr"
                        & size == "large") +
  geom_text(aes(y = udi.best.ll - sqrt(udi.best.ll*diff(range(iso.udi$udi.best.ll, na.rm = TRUE)))/20, label = deciduous_sp)) +
  ylab("udi.best.ll")
g2.all <- g1 %+% subset(iso.udi, tlplevel == "commlevel" & dryseason == "Jan"
                          & size == "large") +
  geom_text(aes(y = udi.best.ll - sqrt(udi.best.ll*diff(range(iso.udi$udi.best.ll, na.rm = TRUE)))/20, label = deciduous_sp)) +
  ylab("udi.best.ll")
grid.arrange(g2.all, g1.all, ncol = 2)

## splevel, best.max.ll
g1.0 <- g1 + # geom_jitter(aes(shape = deciduous, y = udi.best.ll), size = 3, width = 0.25) +
geom_text(aes(y = udi.best.ll - sqrt(udi.best.ll*diff(range(iso.udi$udi.best.ll, na.rm = TRUE)))/20, label = deciduous_sp)) +
  ylab("udi.best.ll")
g2.0 <- g1 %+% subset(iso.udi, tlplevel == "splevel" & dryseason == "Jan" &
                      max.likelihood.norm.tlp == 1 & size == "large") +
  geom_text(aes(y = udi.best.ll - sqrt(udi.best.ll*diff(range(iso.udi$udi.best.ll, na.rm = TRUE)))/20, label = deciduous_sp)) +
  ylab("udi.best.ll")
grid.arrange(g2.0, g1.0, ncol = 2)
## commlevel
g1.1 <- g1 %+% subset(iso.udi, tlplevel == "commlevel" & dryseason == "Mar-Apr" &
                                  max.likelihood.norm.tlp == 1 & size == "large") +
  geom_text(aes(y = udi.best.ll - sqrt(udi.best.ll*diff(range(iso.udi$udi.best.ll, na.rm = TRUE)))/20, label = deciduous_sp)) +
  ylab("udi.best.ll")

g2.1 <- g1 %+% subset(iso.udi, tlplevel == "commlevel" & dryseason == "Jan" &
                        max.likelihood.norm.tlp == 1 & size == "large") +
  geom_text(aes(y = udi.best.ll - sqrt(udi.best.ll*diff(range(iso.udi$udi.best.ll, na.rm = TRUE)))/20, label = deciduous_sp)) +
  ylab("udi.best.ll")

grid.arrange(g2.1, g1.1, ncol = 2)
## including medium size
g1.2 <- g1 %+% subset(iso.udi, tlplevel == "splevel" & dryseason == "Mar-Apr" &
                        max.likelihood.norm.tlp == 1) +
  geom_text(aes(y = udi.best.ll - sqrt(udi.best.ll*diff(range(iso.udi$udi.best.ll, na.rm = TRUE)))/20, label = deciduous_sp))

g2.2 <- g1 %+% subset(iso.udi, tlplevel == "splevel" & dryseason == "Jan" &
                        max.likelihood.norm.tlp == 1) +
  geom_text(aes(y = udi.best.ll - sqrt(udi.best.ll*diff(range(iso.udi$udi.best.ll, na.rm = TRUE)))/20, label = deciduous_sp))

grid.arrange(g2.2, g1.2, ncol = 2)
## max.rsq, splevel
g1.3 <- g1 %+% subset(iso.udi, tlplevel == "splevel" & dryseason == "Mar-Apr" &
                        max.rsq.norm.tlp == 1 & size == "large") +
  geom_text(aes(y = udi.best.rsq - sqrt(udi.best.rsq*diff(range(iso.udi$udi.best.rsq, na.rm = TRUE)))/20, label = deciduous_sp)) +
  ylab("udi.best.rsq")

g2.3 <- g1 %+% subset(iso.udi, tlplevel == "splevel" & dryseason == "Jan" &
                        max.rsq.norm.tlp == 1 & size == "large") +
  geom_text(aes(y = udi.best.rsq - sqrt(udi.best.rsq*diff(range(iso.udi$udi.best.rsq, na.rm = TRUE)))/20, label = deciduous_sp)) +
  ylab("udi.best.rsq")
grid.arrange(g2.3, g1.3, ncol = 2)
## max.rsq, commlevel
g1.4 <- g1 %+% subset(iso.udi, tlplevel == "commlevel" & dryseason == "Mar-Apr" &
                        max.rsq.norm.tlp == 1 & size == "large") +
  geom_text(aes(y = udi.best.rsq - sqrt(udi.best.rsq*diff(range(iso.udi$udi.best.rsq, na.rm = TRUE)))/20, label = deciduous_sp)) +
  ylab("udi.best.rsq")

g2.4 <- g1 %+% subset(iso.udi, tlplevel == "commlevel" & dryseason == "Jan" &
                        max.rsq.norm.tlp == 1 & size == "large") +
  geom_text(aes(y = udi.best.rsq - sqrt(udi.best.rsq*diff(range(iso.udi$udi.best.rsq, na.rm = TRUE)))/20, label = deciduous_sp)) +
  ylab("udi.best.rsq")
grid.arrange(g2.4, g1.4, ncol = 2)


# figure2 <- multi_panel_figure(columns = 5, rows = 2, panel_label_type = "upper-roman")
# figure2 %<>%
#   fill_panel(g2, column = 1:2, row = 1:2) %<>%
#   fill_panel(g1, column = 3:5, row = 1:2)
# figure2
}

