rm(list = ls())
gc()
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, data.table)

summarise.udi <- function(splevel, dryseason, rsq.thresh) {
  load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
  load(file.path("results/4.1GLUEsetup_part2_BCI.RData")) # has n.ensembles and growth and si matrix

  intervals <- info$intervals
  n.ensembles <- growth_by_si.info$n.ensembles
  growth.type <- growth_by_si.info$growth.type ## this is not data but type such as "individual"
  si.type <- growth_by_si.info$si.type
  growth.selection <- growth_by_si.info$growth.selection
  dbh.residuals <- growth_by_si.info$dbh.residuals
  goodness.fit <- 0.3 # rsq0.3
  rsq.thresh <-  rsq.thresh
  n.best.fit <- 100
  if (splevel == "on") {
    level.folder <- "splevel"
  } else {
    level.folder <- "commlevel"
  }

  load(file = paste("results/", level.folder, "/ds.bestfit.all_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_id_dryseason_", dryseason, ".Rdata", sep = ""))
  ds.bestfit <- ds.bestfit.all %>% group_by(sp_size) %>%
    arrange(sp_size, rsq) %>%
    # limiting udi to the cases that are top 100
    mutate(sdi.test = ifelse(rsq < rsq.thresh, NA, sdi),
           udi.test = ifelse(rsq < rsq.thresh, NA, udi),
           rsq.test = ifelse(rsq < rsq.thresh, NA, rsq),
           udi.best = ifelse(rsq == max(rsq, na.rm = TRUE), udi, NA),
           sdi.best = ifelse(rsq == max(rsq, na.rm = TRUE), sdi, NA)) %>%
    summarise(udi.best = mean(udi.best, na.rm = TRUE),
              udi = median(udi, na.rm = TRUE), #sum(udi.test*rsq.test, na.rm = TRUE)/sum(rsq.test, na.rm = TRUE),
              udi.sd = sd(udi.test, na.rm = TRUE),
              n = n(),
              udi.se = udi.sd/(n^0.5),
              sdi.best = mean(udi.best, na.rm = TRUE),
              sdi = sum(sdi.test*rsq.test, na.rm = TRUE)/sum(rsq.test, na.rm = TRUE),
              sdi.sd = sd(sdi.test, na.rm = TRUE),
              sdi.se = sdi.sd/(n^0.5),
              rsq.max = max(rsq.test, na.rm = TRUE),
              rsq.mean = mean(rsq.test, na.rm = TRUE),
              rsq.min = min(rsq.test, na.rm = TRUE),
              root.95 = mean(root.95, na.rm = TRUE),
              max.root = mean(max.root, na.rm = TRUE)) %>%
    ungroup(sp_size) %>%
    separate(sp_size, c("sp", "size"), remove = FALSE) %>%
    mutate(size = factor(size, levels = c("tiny", "small", "medium", "large")))

  ds.bestfit.longer <- ds.bestfit.all %>% pivot_longer(cols = starts_with("depth."),
                                                       names_to = "depth",
                                                       names_prefix = "depth.",
                                                       values_to = "root.frac") %>%
    mutate(depth = as.numeric(depth),
           sp_size_par.sam_rf.sam = paste(sp_size, par.sam, rf.sam, sep = "_")) %>%
    group_by(sp_size_par.sam_rf.sam) %>%
    arrange(sp_size_par.sam_rf.sam, depth) %>%
    mutate(cum.root.frac = cumsum(root.frac)) %>%
    ungroup(sp_size_par.sam_rf.sam)
  head(ds.bestfit.longer)
  summary(ds.bestfit.longer)

  if (splevel == "on") {
    ds.bestfit$tlplevel <- ds.bestfit.longer$tlplevel <- "sp"
  } else {
    ds.bestfit$tlplevel <- ds.bestfit.longer$tlplevel <- "comm"
  }
  save(ds.bestfit, file = paste("results/", level.folder, "/ds.bestfit_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_id_dryseason_", dryseason, ".Rdata", sep = ""))
  save(ds.bestfit.longer, file = paste("results/", level.folder, "/ds.bestfit.longer_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_id_dryseason_", dryseason, ".Rdata", sep = ""))
  # load(file = paste("results/splevel/ds.bestfit_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", intervals, "_id.Rdata", sep = ""))
  # load(file = paste("results/splevel/ds.bestfit.longer_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", intervals, "_id.Rdata", sep = ""))
}
