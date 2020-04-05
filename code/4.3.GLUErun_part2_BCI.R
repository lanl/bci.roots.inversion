rm(list = ls())
gc()
if (!require("pacman"))
  install.packages("pacman")
library(pacman)
pacman::p_load(tidyverse, doParallel, foreach)

load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
intervals <- info$intervals
load("results/GLUEsetup_part1.2_BCI.RData")
load("results/GLUEsetup_part1.3_BCI.RData")
load(file.path("results/4.1GLUEsetup_part2_BCI.RData"))
n.ensembles <- growth_by_si.info$n.ensembles
source(file.path("code/4.2.GLUErun_part2_func_BCI.R"))

GLUE_run <- function(splevel, goodness.fit, drop.months) {

  if (splevel == "on") {
    level.folder <- "splevel"
  } else {
    level.folder <- "commlevel"
  }
  if(!dir.exists(file.path("results", level.folder))) {dir.create(file.path("results", level.folder))}

  file.suffix <-
    paste0("drop.months", drop.months, "_cor", goodness.fit, "_", growth_by_si.info$si.type, "_", n.ensembles, "_", growth_by_si.info$growth.type, "_",
           growth_by_si.info$growth.selection, "_",
           growth_by_si.info$dbh.residuals, "_",
           info$intervals, ".Rdata")
  ###-----------------
  ###---------GLUE.rsq----
  ###-----------------
  beg <- Sys.time()
  beg
  # Function from source(file.path("code/4.2.GLUErun_part2_func_BCI.R"))
  GLUE.rsq <-
    growth_by_si.func(
      fun.version = "rsq",
      growth_by_si.info = growth_by_si.info,
      statistic = "rsq",
      splevel = splevel,
      goodness.fit = goodness.fit,
      drop.months = drop.months
    )
  # fun.version = "rsq"; growth_by_si.info = growth_by_si.info; statistic = "rsq"; splevel = splevel;
  # goodness.fit = goodness.fit; drop.months = drop.months
  end = Sys.time()
  (end - beg) * nrow(growth_by_si.info$growth) / 60 / 60 / 100
  end - beg
  # 24 min on macbook
  # 26 secs at both levels ; 47.21 sec for parallel at lower levels # 10 n.ensembles & 12 growth
  ## took 15 hours on the server for 1k n.ensembles and nrow(growth) = 158107
  # takes 11 min on server for ~580 growth ts n.ensembles = 10000 n.best = 100
  save(GLUE.rsq,
       file = paste0("results/", level.folder, "/GLUE.rsq_", file.suffix))
  rm(GLUE.rsq)
  ##---------------------------------------------------------------------------------
  ##---------GLUE.matches----No. of growth predictions that are within 95% CI of median growth
  ##---------------------------------------------------------------------------------
  # Function from source(file.path("code/4.2.GLUErun_part2_func2_BCI.R"))

  GLUE.matches <-
    growth_by_si.func(
      fun.version = "likelihood",
      growth_by_si.info = growth_by_si.info,
      statistic = "rsq", ## this is actually needed for GLUE.rsq
      splevel = splevel,
      goodness.fit = goodness.fit,
      drop.months = drop.months
    )
  # fun.version = "likelihood"; growth_by_si.info = growth_by_si.info; statistic = "rsq"; splevel = splevel;
  # goodness.fit = goodness.fit; drop.months = drop.months
  save(GLUE.matches,
       file = paste0("results/", level.folder, "/GLUE.matches_", file.suffix))
  rm(GLUE.matches)
  }

# for profiling:
# Rprof(append = FALSE)
# growth_by_si.aprof <- aprof("code/4.2GLUErun_part2_BCI.R","growth_by_si.out")
#
# plot(growth_by_si.aprof)
