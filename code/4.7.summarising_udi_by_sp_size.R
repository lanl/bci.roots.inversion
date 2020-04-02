rm(list = ls())
gc()
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, data.table)

summarise.udi <- function(splevel, dryseason, rsq.thresh, root.selection, iso.subset) {
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
  n.best.fit <- 100
  if (splevel == "on") {
    level.folder <- "splevel"
  } else {
    level.folder <- "commlevel"
  }

  file.extension.base3 <- paste0(goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_dryseason_", dryseason, "_iso.subset_", iso.subset)
  file.extension.base4 <- paste0(file.extension.base3, "_root.selection_", root.selection)

  load(file = paste0("results/", level.folder, "/ds.bestfit.all_cor", file.extension.base3, ".Rdata"))

  if (root.selection == "on") {
    select.rf.sam <- read.csv(file = file.path("results/rf.sam_power.threshold_0.5.csv"), header = TRUE)
    ds.bestfit.all <- ds.bestfit.all %>% subset(rf.sam %in% select.rf.sam$x)
  } # 702 root profiles selected
    ds.bestfit.select <- ds.bestfit.all %>%
      #subset(!is.na(udi)) %>%
      group_by(sp_size) %>%
    arrange(sp_size, rsq) %>%
    # need to find uncertainty in top-100 avaiable udi
    mutate(rsq.rank = order(order(rsq, decreasing = TRUE)),
           sdi.tops = ifelse(rsq < rsq.thresh | rsq.rank > 500, NA, sdi),
           udi.tops = ifelse(rsq < rsq.thresh | rsq.rank > 500, NA, udi),
           rsq.tops = ifelse(rsq < rsq.thresh | rsq.rank > 500, NA, rsq),
           root.95.tops = ifelse(rsq < rsq.thresh | rsq.rank > 500, NA, root.95),
           root.75.tops = ifelse(rsq < rsq.thresh | rsq.rank > 500, NA, root.75),
           max.root.tops = ifelse(rsq < rsq.thresh | rsq.rank > 500, NA, max.root),
           udi.best = ifelse(rsq.rank == 1, udi, NA),
           sdi.best = ifelse(rsq.rank == 1, sdi, NA)) %>%
      subset(!is.na(rsq.tops)) %>%
      arrange(sp_size, rsq.rank)
    ds.bestfit <- ds.bestfit.select %>%
    summarise(udi.best = mean(udi.best, na.rm = TRUE),
              udi = median(udi.tops, na.rm = TRUE), #sum(udi.tops*rsq.tops, na.rm = TRUE)/sum(rsq.tops, na.rm = TRUE),
              udi.sd = sd(udi.tops, na.rm = TRUE),
              n = n(),
              udi.se = udi.sd/(n^0.5),
              sdi.best = mean(sdi.best, na.rm = TRUE),
              sdi = median(sdi.tops, na.rm = TRUE), ## sum(sdi.tops*rsq.tops, na.rm = TRUE)/sum(rsq.tops, na.rm = TRUE),
              sdi.sd = sd(sdi.tops, na.rm = TRUE),
              sdi.se = sdi.sd/(n^0.5),
              rsq.max = max(rsq.tops, na.rm = TRUE),
              rsq.mean = mean(rsq.tops, na.rm = TRUE),
              rsq.min = min(rsq.tops, na.rm = TRUE),
              root.95 = median(root.95.tops, na.rm = TRUE),
              root.75 = median(root.75.tops, na.rm = TRUE),
              max.root = median(max.root.tops, na.rm = TRUE)) %>%
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
  # if (splevel == "on") {
  #   ds.bestfit$tlplevel <- "sp"
  # } else {
  #   ds.bestfit$tlplevel <- "comm"
  # }
    if (splevel == "on") {
      ds.bestfit$tlplevel <- ds.bestfit.longer$tlplevel <- "sp"
    } else {
      ds.bestfit$tlplevel <- ds.bestfit.longer$tlplevel <- "comm"
    }
  save(ds.bestfit, file = paste0("results/", level.folder, "/ds.bestfit_cor", file.extension.base4, ".Rdata"))
  save(ds.bestfit.longer, file = paste0("results/", level.folder, "/ds.bestfit.longer_cor", file.extension.base4, ".Rdata"))
  # load(file = paste0("results/splevel/ds.bestfit_cor", file.extension.base4, ".Rdata"))
  # load(file = paste0("results/splevel/ds.bestfit.longer_cor", file.extension.base4, ".Rdata"))
}
