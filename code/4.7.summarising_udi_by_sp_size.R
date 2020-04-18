rm(list = ls())
gc()
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, data.table)

summarise.udi <- function(splevel = splevel, dryseason = dryseason, rsq.thresh = rsq.thresh,
                          root.selection = root.selection, iso.subset = iso.subset,
                          n.rank = n.rank, drop.months = drop.months) {
  load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
  load(file.path("results/4.1GLUEsetup_part2_BCI.RData")) # has n.ensembles and growth and si matrix

  intervals <- info$intervals
  rm(info)
  n.ensembles <- growth_by_si.info$n.ensembles
  growth.type <- growth_by_si.info$growth.type ## this is not data but type such as "individual"
  si.type <- growth_by_si.info$si.type
  growth.selection <- growth_by_si.info$growth.selection
  dbh.residuals <- growth_by_si.info$dbh.residuals
  goodness.fit <- 0.3 # rsq0.3
  rsq.thresh <-  rsq.thresh
  if (splevel == "on") {
    level.folder <- "splevel"
  } else {

    level.folder <- "commlevel"
  }

  file.extension.base3 <- paste0("drop.months", drop.months, "_cor", goodness.fit, "_",
                                 si.type, "_", n.ensembles, "_", growth.type, "_",
                                 growth.selection, "_", dbh.residuals, "_", intervals,
                                 "_dryseason_", dryseason, "_iso.subset_", iso.subset)
  file.extension.base4 <- paste0(file.extension.base3, "_root.selection_", root.selection)

  load(file = paste0("results/", level.folder, "/ds.bestfit.all_", file.extension.base3, ".Rdata"))

  if (root.selection == "on") {
    select.rf.sam <- read.csv(file = file.path("results/rf.sam_power.threshold_0.5.csv"), header = TRUE)
    ds.bestfit.all <- ds.bestfit.all %>% subset(rf.sam %in% select.rf.sam$x)
  } # 702 root profiles selected
    ds.bestfit.select <- ds.bestfit.all %>%
      group_by(sp_size) %>%
    arrange(sp_size, rsq) %>%
    # need to find uncertainty in top n.rank avaiable udi
      ## but UDI is not normally distributed but lognormally distributed
    mutate(rsq.rank = rank(-rsq, ties.method = "average"),
           ll.rank = rank(neg.likelihood, ties.method = "average"),
           udi.best.rsq = ifelse(rsq.rank == 1, udi, NA),
           sdi.best.rsq = ifelse(rsq.rank == 1, sdi, NA),
           udi.best.ll = ifelse(ll.rank == 1, udi, NA),
           sdi.best.ll = ifelse(ll.rank == 1, sdi, NA),
           udi.tops.rsq = ifelse(rsq < rsq.thresh | rsq.rank > n.rank, NA, udi),
           sdi.tops.rsq = ifelse(rsq < rsq.thresh | rsq.rank > n.rank, NA, sdi),
           rsq.tops.rsq = ifelse(rsq < rsq.thresh | rsq.rank > n.rank, NA, rsq),
           udi.tops.ll = ifelse(ll.rank > n.rank, NA, udi),
           sdi.tops.ll = ifelse(ll.rank > n.rank, NA, sdi),
           rsq.tops.ll = ifelse(ll.rank > n.rank, NA, rsq),
           root.95.tops.rsq = ifelse(rsq < rsq.thresh | rsq.rank > n.rank, NA, root.95),
           root.75.tops.rsq = ifelse(rsq < rsq.thresh | rsq.rank > n.rank, NA, root.75),
           max.root.tops.rsq = ifelse(rsq < rsq.thresh | rsq.rank > n.rank, NA, max.root),
           root.95.tops.ll = ifelse(ll.rank > n.rank, NA, root.95),
           root.75.tops.ll = ifelse(ll.rank > n.rank, NA, root.75),
           max.root.tops.ll = ifelse(ll.rank > n.rank, NA, max.root)) %>%
      subset(!is.na(rsq.tops.rsq)) %>%
      arrange(sp_size, ll.rank, rsq.rank)
    ds.bestfit <- ds.bestfit.select %>%
    summarise(udi.best.rsq = mean(udi.best.rsq, na.rm = TRUE),
              udi.best.ll = mean(udi.best.ll, na.rm = TRUE),
              udi.med.rsq = median(udi.tops.rsq, na.rm = TRUE), #sum(udi.tops*rsq.tops, na.rm = TRUE)/sum(rsq.tops, na.rm = TRUE),
              udi.upr.rsq = quantile(udi.tops.rsq, probs = 0.95, na.rm = TRUE),
              udi.lwr.rsq = quantile(udi.tops.rsq, probs = 0.05, na.rm = TRUE),
              udi.sd.rsq = sd(log(udi.tops.rsq), na.rm = TRUE),
              udi.med.ll = median(udi.tops.ll, na.rm = TRUE), #sum(udi.tops*rsq.tops, na.rm = TRUE)/sum(rsq.tops, na.rm = TRUE),
              udi.upr.ll = quantile(udi.tops.ll, probs = 0.975, na.rm = TRUE),
              udi.lwr.ll = quantile(udi.tops.ll, probs = 0.025, na.rm = TRUE),
              udi.sd.ll = sd(log(udi.tops.ll), na.rm = TRUE),
              n = n(),
              sdi.best.rsq = mean(sdi.best.rsq, na.rm = TRUE),
              sdi.med.rsq = median(sdi.tops.rsq, na.rm = TRUE), ## sum(sdi.tops*rsq.tops, na.rm = TRUE)/sum(rsq.tops, na.rm = TRUE),
              sdi.sd.rsq = sd(log(sdi.tops.rsq), na.rm = TRUE),
              sdi.best.ll = mean(sdi.best.ll, na.rm = TRUE),
              sdi.med.ll = median(sdi.tops.ll, na.rm = TRUE), ## sum(sdi.tops*rsq.tops, na.rm = TRUE)/sum(rsq.tops, na.rm = TRUE),
              sdi.sd.ll = sd(log(sdi.tops.ll), na.rm = TRUE),

              root.95.rsq = median(root.95.tops.rsq, na.rm = TRUE),
              root.75.rsq = median(root.75.tops.rsq, na.rm = TRUE),
              max.root.rsq = median(max.root.tops.rsq, na.rm = TRUE),
              rsq.max.ll = max(rsq.tops.ll, na.rm = TRUE),
              rsq.mean.ll = mean(rsq.tops.ll, na.rm = TRUE),
              rsq.min.ll = min(rsq.tops.ll, na.rm = TRUE),
              root.95.ll = median(root.95.tops.ll, na.rm = TRUE),
              root.75.ll = median(root.75.tops.ll, na.rm = TRUE),
              max.root.ll = median(max.root.tops.ll, na.rm = TRUE)) %>%
    ungroup(sp_size) %>%
    separate(sp_size, c("sp", "size"), remove = FALSE) %>%
    mutate(size = factor(size, levels = c("tiny", "small", "medium", "large")))

  ds.bestfit.longer <- ds.bestfit.select %>%
    ungroup(sp_size) %>%
    mutate(depth.0.00 = 0) %>%
    pivot_longer(cols = starts_with("depth."), names_to = "depth",
                 names_prefix = "depth.", values_to = "root.frac") %>%
    mutate(depth = as.numeric(depth),
           sp_size_par.sam_rf.sam = paste(sp_size, par.sam, rf.sam, sep = "_")) %>%
    group_by(sp_size_par.sam_rf.sam) %>%
    arrange(sp_size_par.sam_rf.sam, depth) %>%
    mutate(cum.root.frac = cumsum(root.frac)) %>%
    ungroup(sp_size_par.sam_rf.sam) %>%
    subset(!is.na(udi))
  head(ds.bestfit.longer)
  summary(ds.bestfit.longer)


  ds.bestfit.select <- ds.bestfit.select %>% mutate(udi.tops.rsq = list(udi.tops.rsq),
                                                    udi.tops.ll = list(udi.tops.ll))
    if (splevel == "on") {
      ds.bestfit$tlplevel <- ds.bestfit.longer$tlplevel <- "sp"
    } else {
      ds.bestfit$tlplevel <- ds.bestfit.longer$tlplevel <- "comm"
    }
  save(ds.bestfit, file = paste0("results/", level.folder, "/ds.bestfit_", file.extension.base4, ".Rdata"))
  save(ds.bestfit.longer, file = paste0("results/", level.folder, "/ds.bestfit.longer_", file.extension.base4, ".Rdata"))
  # load(file = paste0("results/splevel/ds.bestfit_", file.extension.base4, ".Rdata"))
  # load(file = paste0("results/splevel/ds.bestfit.longer_", file.extension.base4, ".Rdata"))
}
