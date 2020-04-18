#----------------------------------------------
# Title: Collate growth and btran data and metadata
# Author : Rutuja Chitra-Tarak
# Original date: Feb 18, 2019
#----------------------------------------------

rm(list = ls())
gc()
graphics.off()
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse)

GLUEsetup_part2 <- function(growth.type, dbh.residuals, solar.residuals, growth.selection) { #n.ensembles = 5000
  # n.ensembles = 100000,  #which ensemble set to use
    #growth.selection  <- "size_class_predefined_cc_scaled"
    load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
    intervals <- info$intervals
    ## load Result of 3.GLUErun_part1_BCI.R:
    load("results/GLUEsetup_part1.3_BCI.RData")
    n.ensembles <- length(info.3$btran.int.sam[[1]])
    ## Observed growth
    ## individual, species or sp by size level
    growth.name <- load(file =  paste0("results/sp_size.", growth.type, "_growth_dbh.residuals_", dbh.residuals, "_ci_", intervals, "_", growth.selection, ".Rdata"))
    growth <- get(growth.name); rm(growth.name)
    sp_size.name <- load(file = paste0("results/sp_size.", growth.type, ".names_", intervals, "_", growth.selection, ".Rdata"))
    growth.meta <- data.frame(sp_size = get(sp_size.name)) %>%
      separate(sp_size, into = c("sp", "size"), sep = "_", remove = FALSE)
    rm(sp_size.name)

    # growth.type <- "med_growth_sp_by_size"
    # load(paste0("results/", growth.type, "_intervals_", intervals ,".Rdata"))
    # growth.long <- med_growth_sp_by_size %>% unite("sp_size", sp, size, remove = TRUE) %>% select(-interval.2) %>%
    #   pivot_longer(cols = starts_with("interval."), names_to = "interval", names_prefix = "interval.", values_to = "growth") %>% mutate(interval = as.numeric(interval))
    # growth.prep <- split(growth.long, growth.long$sp_size)
    # growth <- lapply(growth.prep, function(x) {x %>% select(-sp_size)})
    #
    # sp_size.median.names <- names(growth)
    # growth.meta <- data.frame(sp_size = sp_size.median.names) %>%
    #   separate(sp_size, into = c("sp", "size"), sep = "_", remove = FALSE)


    # growth.type <- "med_growth_sp_by_size"
    ## species or size level
    ## these have centered and scaled median growth data for six intervals,
    # with those species retained for which at least 3 intervals in the last
    # (5) "intervals had at least 5 trees/interval; median growth is NA for intervals with insufficient sample sizes

    # load(paste0("results/", growth.type, "__intervals_", intervals ,".Rdata"))
    # ## only one of the following two will work at a time
    # growth <- med_growth_sp_by_size
    # growth <- med_growth_by_size

    # activecols.growth <- (ncol(growth)-intervals + 1):ncol(growth)
    # activecols.btran <- (ncol(btran[[1]])-intervals + 1):ncol(btran[[1]])
    ## convert growth with activecols to vector; then to matrix
    # growth.matrix <- as.matrix(growth[, activecols.growth], nrow = nrow(growth),
    #                            ncol = length(activecols.growth), byrow = T)

    growth_by_si.info <- list(
      n.ensembles = n.ensembles,
      dbh.residuals = dbh.residuals,
      solar.residuals = solar.residuals,
      growth.selection = growth.selection,
      growth.type = growth.type,
      growth.meta = growth.meta,
      growth = growth,
      si = info.3$btran.matrix,
      si.param.rel = info.3$si.param.rel,
      drop.months.vec = info.3$drop.months.vec,
      si.type = info$si.type
    )
    save(file = "results/4.1GLUEsetup_part2_BCI.RData", growth_by_si.info)
    # load(file.path("results/4.1GLUEsetup_part2_BCI.RData"))
}
