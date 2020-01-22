
  rm(list = ls())
  gc()
  if (!require("pacman")) install.packages("pacman"); library(pacman)
  pacman::p_load(tidyverse, doParallel, foreach)

  load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
  intervals <- info$intervals
  load("results/GLUEsetup_part1.5_BCI.RData")
  load(file.path("results/4.1GLUEsetup_part2_BCI.RData"))
  n.ensembles <- growth_by_si.info$n.ensembles
  source(file.path("R/4.2.GLUErun_part2_func_BCI.R"))
  GLUE_run <- function(splevel, goodness.fit){
    if(splevel == "on") {
      sp.GLUE.rsq <- growth_by_si.func(growth_by_si.info = growth_by_si.info,
                                       statistic = "rsq", splevel = "on", goodness.fit)

      save(sp.GLUE.rsq, file = paste0("results/sp.GLUE.rsq_",
                                      growth_by_si.info$si.type, "_", n.ensembles, "_",
                                      growth_by_si.info$growth.type, "_",
                                      growth_by_si.info$growth.selection, "_",
                                      growth_by_si.info$dbh.residuals, "_",
                                      info$intervals, "_cor", goodness.fit, ".Rdata"))
      rm(sp.GLUE.rsq)

      # sp.GLUE.cor <- growth_by_si.func(growth_by_si.info = growth_by_si.info,
      #                                  statistic = "cor", splevel = "on", goodness.fit)
      # save(sp.GLUE.cor, file = paste0("results/sp.GLUE.cor_", growth_by_si.info$si.type, "_", n.ensembles, "_",
      #                                 growth_by_si.info$growth.type, "_",
      #                                 info$intervals, "_cor", goodness.fit, ".Rdata"))

    } else {
    # Rprof(file = "growth_by_si.out", interval = 0.02,
    #       line.profiling = TRUE)
      beg <- Sys.time()
      beg
      GLUE.rsq <- growth_by_si.func(growth_by_si.info = growth_by_si.info,
                                    statistic = "rsq", splevel = "off", goodness.fit)
      end = Sys.time();  (end - beg)*nrow(growth_by_si.info$growth)/60/60/100
      end - beg
      # 24 min on macbook
      # 26 secs at both levels ; 47.21 sec for parallel at lower levels # 10 n.ensembles & 12 growth
      ## took 15 hours on the server for 1k n.ensembles and nrow(growth) = 158107
      # takes 11 min on server for ~580 growth ts n.ensembles = 10000 n.best =100
      save(GLUE.rsq, file = paste0("results/GLUE.rsq_", growth_by_si.info$si.type, "_", n.ensembles, "_",
                                   growth_by_si.info$growth.type, "_",
                                   growth_by_si.info$growth.selection, "_",
                                   growth_by_si.info$dbh.residuals, "_",
                                   info$intervals, "_cor", goodness.fit, ".Rdata"))
      rm(GLUE.rsq)

      # GLUE.cor <- growth_by_si.func(growth_by_si.info = growth_by_si.info,
      #                               statistic = "cor", splevel = "off", goodness.fit)
      # save(GLUE.cor, file = paste0("results/GLUE.cor_", growth_by_si.info$si.type, "_", n.ensembles, "_",
      #                              growth_by_si.info$growth.type, "_",
      #                              info$intervals, "_cor", goodness.fit, ".Rdata"))
    }
  }

  # for profiling:
  # Rprof(append = FALSE)
  # growth_by_si.aprof <- aprof("R/4.2GLUErun_part2_BCI.R","growth_by_si.out")
  #
  # plot(growth_by_si.aprof)
