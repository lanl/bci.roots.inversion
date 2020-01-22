## Finding best fit water-stress/growth model parameters for species growth time series

# .rs.restartR()
# for 50 ha obs species groups
rm(list = ls())
gc()
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, doParallel, foreach, data.table)

prep.udi.calculator <- function(splevel) {
  # load interval and n.ensembles
  load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
  load(file.path("results/4.1GLUEsetup_part2_BCI.RData")) # has n.ensembles and growth and si matrix

  intervals <- info$intervals
  n.ensembles <- growth_by_si.info$n.ensembles
  growth.type <- growth_by_si.info$growth.type ## this is not data but type such as "individual"
  growth.selection <- growth_by_si.info$growth.selection
  dbh.residuals <- growth_by_si.info$dbh.residuals
  si.type <- growth_by_si.info$si.type
  goodness.fit <- 0.3 # rsq0.3
  rm(growth_by_si.info)

  ######----------------------------------------------
  ###### Preparing data tables for uptake depth index calculation
  ######----------------------------------------------

  load(file = file.path("data-raw/paw.rda"))
  paw.df <- paw %>% data.table()
  paw.df <- paw.df[, date := as.factor(date)]
  rm(paw)
  rdt <- info$root.param %>%
    pivot_longer(cols = starts_with("depth."), names_to = "depth",
                 names_prefix = "depth.", values_to = "root.frac") %>%
    mutate(depth = as.numeric(depth)) %>%
    data.table(key = "depth")
  depth.thick <- data.frame(depth = as.numeric(unique(rdt$depth))) %>%
    mutate(thickness = depth - lag(depth, default = 0)) %>% data.table()

  if (splevel == "on") {
    load(file.path("results/sp.hydro.output.chosen.RData"))
    sp.hydro.output.chosen <- lapply(sp.hydro.output.chosen, function(x){
      data.table(x, key = "depth")[, date := as.factor(date)]
    })
    sdt.list <- lapply(sp.hydro.output.chosen, function(x){
      setkey(x[depth.thick, on = "depth", `:=`(thickness = i.thickness)][
        paw.df, on = c("date","depth","par.sam"), `:=`(paw = i.paw)][,
                                                                     tw := paw*thickness*1000][, !c("paw", "thickness")], depth)
    })
    save(sdt.list, file = paste("results/splevel/sdt.list_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_dryseason_off.Rdata", sep = ""))
    rm(sp.hydro.output.chosen)
    ### Subsetting to dryseason----
    sdt.list <- lapply(sdt.list, function(x) {
      x[, ':='(month = format(as.Date(date), "%m"), year = format(as.Date(date), "%Y"))][month %chin% c("03", "04") & year == "1997"] ## Choosing Feb, March and April
    })
    sdt.list <- lapply(sdt.list,  droplevels)
    save(sdt.list, file = paste("results/splevel/sdt.list_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_dryseason_on.Rdata", sep = ""))
  } else {
    hydro <- info$hydro %>% mutate(date = as.factor(date)) %>% data.table()
    sdt <- setkey(hydro[depth.thick, on = "depth", `:=`(thickness = i.thickness)][
      paw.df, on = c("date","depth","par.sam"), `:=`(paw = i.paw)][,
                                                                   tw := paw*thickness*1000][, !c("paw", "thickness")], depth)
    rm(hydro)
    save(sdt, file = paste("results/commlevel/sdt_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_dryseason_off.Rdata", sep = ""))
    ### Subsetting to dryseason----
    sdt <- sdt[, ':='(month = format(as.Date(date), "%m"), year = format(as.Date(date), "%Y"))][month %chin% c("03", "04") & year == "1997"]
    sdt <- sdt %>% droplevels()
    save(sdt, file = paste("results/commlevel/sdt_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_dryseason_on.Rdata", sep = ""))
  }
  rm(info, paw.df)
  save(rdt, file = "results/rdt.Rdata")
}
