## Best-fit_param_BCI_fun
## load btran
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
    btran <- info$btran
    n.ensembles <- length(btran)
    ## Observed growth
    ## individual, species or sp by size level
    growth.name <- load(file =  paste0("results/sp_size.", growth.type, "_growth_dbh.residuals_", dbh.residuals, "_", intervals, "_", growth.selection, ".Rdata"))
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
    Sys.time()
    beg <- Sys.time()
    btran.mat.list <- lapply(btran, function(x) {
      btran.wide <- x %>% pivot_wider(names_from = interval, values_from = btran) %>% as.matrix()
      btran.wide[, -1]
    })
    Sys.time()
    end <- Sys.time()
    (end - beg)
    btran.matrix <- do.call(rbind, btran.mat.list)
    ### record relational info for si to rf.sam and par.sam
    btran.wide.1 <- btran[[1]] %>% pivot_wider(names_from = interval, values_from = btran) %>% as.matrix()
    n.best <- nrow(btran.wide.1); par.sam.vec <- btran.wide.1[, "par.sam"]
    si.param.rel <- data.frame(row.num = 1:nrow(btran.matrix),
                         rf.sam = rep(c(1:n.ensembles), each = n.best),
                         par.sam = rep(par.sam.vec, times = n.ensembles))

    growth_by_si.info <- list(
      n.ensembles = n.ensembles,
      dbh.residuals = dbh.residuals,
      solar.residuals = solar.residuals,
      growth.selection = growth.selection,
      growth.type = growth.type,
      growth.meta = growth.meta,
      growth = growth,
      si = btran.matrix,
      si.param.rel = si.param.rel,
      si.type = current.folder
    )
    save(file = "results/4.1GLUEsetup_part2_BCI.RData", growth_by_si.info)
    # load(file.path("results/4.1GLUEsetup_part2_BCI.RData"))
}
