## Finding best fit water-stress/growth model parameters for species growth time series

# .rs.restartR()
# for 50 ha obs species groups
rm(list = ls())
for (i in 1:5) {gc()}

if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, doParallel, foreach, data.table)

udi.calculator <- function(splevel, dryseason,  rsq.thresh) {
  # load interval and n.ensembles
  load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
  load(file.path("results/4.1GLUEsetup_part2_BCI.RData")) # has n.ensembles and growth and si matrix

  intervals <- info$intervals
  n.ensembles <- growth_by_si.info$n.ensembles
  growth.type <- growth_by_si.info$growth.type ## this is not data but type such as "individual"
  growth.selection <- growth_by_si.info$growth.selection
  dbh.residuals <- growth_by_si.info$dbh.residuals
  si.type <- growth_by_si.info$si.type
  si.param.rel <- growth_by_si.info$si.param.rel
  goodness.fit <- 0.3 # rsq0.3
  ncor <- detectCores() - 2

  if (splevel == "on") {
    level.folder <- "splevel"
    } else {
    level.folder <- "commlevel"
    }

  rf.sam <- unique(si.param.rel$rf.sam)
  n.par.sam <- length(unique(si.param.rel$par.sam)) ## this is the sequence repeated in growth_by_si.info$si.param.rel
  par.sam <- si.param.rel$par.sam[1:n.par.sam] ## this is the sequence used when GLUE.rsq was constructed
  n.rf.sam <- length(rf.sam)
  rm(info); rm(growth_by_si.info)
  ## sdt is chosen based on dryseason on/off
  load(file = paste0("results/", level.folder, "/GLUE.rsq.vector_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, ".Rdata"), envir = parent.frame(), verbose = FALSE)
  sp_size <- row.names(GLUE.rsq.vector)
  sp <- stringr::str_split(sp_size, "_", simplify = TRUE)[, 1]
  size <- stringr::str_split(sp_size, "_", simplify = TRUE)[, 2]

  if (splevel == "on") {
    load(file = paste("results/splevel/sdt.list_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_dryseason_", dryseason, ".Rdata", sep = ""), envir = parent.frame(), verbose = FALSE)
  } else {
    load(file = paste("results/commlevel/sdt_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_dryseason_", dryseason, ".Rdata", sep = ""), envir = parent.frame(), verbose = FALSE)
  }
  load(file = "results/rdt.Rdata", envir = parent.frame(), verbose = FALSE)

  ######----------------------------------------------
  ###### Calculating uptake depth index
  ######----------------------------------------------

  vec.in.list <- rep(-999, ncol(GLUE.rsq.vector))
  best10.type.udi.list <- best10.type.sdi.list <- rep(list(vec.in.list), nrow(GLUE.rsq.vector))
  #growth_by_si.info$si.param.rel
  Sys.time()
  rows.par.sam <- list()
  for (j in 1: length(par.sam)) {
    rows.par.sam[[j]] <- which(si.param.rel$par.sam == par.sam[j])
  }
  GLUE.rsq.vector.thresh <- as.matrix(GLUE.rsq.vector)
  GLUE.rsq.vector.thresh <- as.matrix(ifelse(GLUE.rsq.vector.thresh < rsq.thresh, NA, GLUE.rsq.vector.thresh))

  beg <- Sys.time()
  if (splevel == "on") { ## foreach at splevel exhausts memory, hence a for loop
    for (ii in 1:  length(sp_size)) {
        ## for each row of sp_size in best10.type, choose sp.hydro.output.chosen for the corresponding sp
        sdt <- sdt.list[sp[ii]][[1]]
      for (j in 1: length(par.sam)) {
        rows.on <- rows.par.sam[[j]]
        ## selecting rf.sam for which GLUE.rsq is above the threshold
        ## For this first selecting values in GLUE.rsq that belong to par.sam; sequence refer to rf.sam
        ## within those selecting values that are above threshold
        rf.sam.j <- rf.sam[!is.na(GLUE.rsq.vector.thresh[ii, rows.on])]
        for (k in 1: length(rf.sam.j)) {
          hydro.rf <- merge(sdt[par.sam == par.sam[j]],
                         rdt[rf.sam %in% rf.sam.j[k]])
          best10.type.sdi.list[[ii]][rows.on][rf.sam.j[k]] <- hydro.rf[, ## location works with [rf.sam.j[k]] because rf.sam values are sequential integers
            .(gr.mean = mean(gfac*root.frac, na.rm = TRUE)), by = depth][,
              .(sdi  = sum(depth*gr.mean, na.rm = TRUE))]
          ### same  as best10.type.sdi.list[[ii]][rows.on[k]] <- hydro.rf[,
          # gfac.rootfrac := gfac*root.frac][,
          #  gr := sum(gfac.rootfrac, na.rm = TRUE), by = .(date)][, ## does not sum to 1
          #    .(gr.mean = mean(gr, na.rm = TRUE)), by = depth][,
          #       .(sdi  = sum(depth*gr, na.rm = TRUE))]
          fw.mean <- hydro.rf[,
            w := root.frac*gfac*tw][,
              wsum := sum(w, na.rm = TRUE), by = .(date)][,
                fw := w/wsum][, ## water uptake profile
                  .(fw.mean = mean(fw, na.rm = TRUE)), by = depth]
          best10.type.udi.list[[ii]][rows.on][rf.sam.j[k]] <- fw.mean[,
            .(udi = sum(fw.mean*depth, na.rm = TRUE))]
        }
      }
      best10.type.sdi.list[[ii]] <- unlist(best10.type.sdi.list[[ii]])
      best10.type.udi.list[[ii]] <- unlist(best10.type.udi.list[[ii]])
    } # tool 2.3 hrs on macbook for ~60 sp
  } else {
    ###++++++++++++++++++++++++
    ### Calculating SDI--------
    ###++++++++++++++++++++++++
    doParallel::registerDoParallel(ncor)
    best10.type.sdi.list <- foreach(ii = 1: length(sp_size)) %dopar% {
      ## each of the numbers in best10.type.list[[ii]]
      for (j in 1: length(par.sam)) {
        rows.on <- rows.par.sam[[j]]
        ## selecting rf.sam for which GLUE.rsq is above the threshold
        ## For this first selecting values in GLUE.rsq that belong to par.sam; sequence refer to rf.sam
        ## within those selecting values that are above threshold
        rf.sam.j <- rf.sam[!is.na(GLUE.rsq.vector.thresh[ii, rows.on])]
        for (k in 1: length(rf.sam.j)) {
          hydro.rf <- merge(sdt[par.sam == par.sam[j]],
                         rdt[rf.sam %in% rf.sam.j[k]])
          best10.type.sdi.list[[ii]][rows.on][rf.sam.j[k]] <- hydro.rf[,
              .(gr.mean = mean(gfac*root.frac, na.rm = TRUE)), by = depth][,
                .(sdi  = sum(depth*gr.mean, na.rm = TRUE))]
        }
      }
      best10.type.sdi.list[[ii]] <- unlist(best10.type.sdi.list[[ii]])
      # 20 min on imac
    }
    stopImplicitCluster()
    #parallel::stopCluster(cl)

    ###++++++++++++++++++++++++
    ### Calculating UDI--------
    ###++++++++++++++++++++++++
    #cl <- parallel::makeForkCluster(ncor) # slower with makeForkCluster
    doParallel::registerDoParallel(ncor)
    best10.type.udi.list <- foreach(ii = 1: length(sp_size)) %dopar% {
      ## each of the numbers in best10.type.list[[ii]]
      for (j in 1: length(par.sam)) {
        rows.on <- rows.par.sam[[j]]
        ## selecting rf.sam for which GLUE.rsq is above the threshold
        ## For this first selecting values in GLUE.rsq that belong to par.sam; sequence refer to rf.sam
        ## within those selecting values that are above threshold
        rf.sam.j <- rf.sam[!is.na(GLUE.rsq.vector.thresh[ii, rows.on])]
        for (k in 1: length(rf.sam.j)) {
          hydro.rf <- merge(sdt[par.sam == par.sam[j]],
                            rdt[rf.sam %in% rf.sam.j[k]])
          fw.mean <- hydro.rf[,
            w := root.frac*gfac*tw][,
              wsum := sum(w, na.rm = TRUE), by = .(date)][,
              fw := w/wsum][, ## water uptake profile
                .(fw.mean = mean(fw, na.rm = TRUE)), by = depth]
          best10.type.udi.list[[ii]][rows.on][rf.sam.j[k]] <- fw.mean[,
            .(udi = sum(fw.mean*depth, na.rm = TRUE))]
        }
      }
      best10.type.udi.list[[ii]] <- unlist(best10.type.udi.list[[ii]])
    }
    stopImplicitCluster()

    # 30 min on imac
  }
  Sys.time()
  end <- Sys.time()
  (end - beg)
  best10.type.sdi <- data.frame(do.call(rbind, best10.type.sdi.list))
  best10.type.sdi[best10.type.sdi == -999] <- NA
  row.names(best10.type.sdi) <- sp_size
  rm(best10.type.sdi.list)
  best10.type.udi <- data.frame(do.call(rbind, best10.type.udi.list))
  best10.type.udi[best10.type.udi == -999] <- NA
  row.names(best10.type.udi) <- sp_size
  rm(best10.type.udi.list)

  ## 5 hrs for 134 sp_size in splevel == "on" on imac ~1500 rf.sam
  ## 4.18 for ~800 sp_size in splevel == "off" on imac ~1500 rf.sam
  ## 2.3 hrs on macbook
  ######----------------------------------------------

  ## the top-10 largest Rsq in decreasing order
  #best10.type.rsq <- t(apply(GLUE.rsq, 1, function(x) {sort(x, decreasing = T)[1:10]}))
  # > length(which(rowSums(is.na(best10.type.rsq)) == ncol(best10.type.rsq)))
  # [1] 732; so 732 spe-size have not found a match with r-sq >= 0.3

  save(best10.type.sdi, file = paste0("results/", level.folder, "/best10.type.sdi_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_dryseason_", dryseason, ".Rdata"))
  save(best10.type.udi, file = paste0("results/", level.folder, "/best10.type.udi_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_dryseason_", dryseason, ".Rdata"))
  # if (splevel == "on") {
  #   load(file = paste("results/splevel/best10.type.sdi_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_dryseason_", dryseason, ".Rdata", sep = ""), envir = parent.frame(), verbose = FALSE)
  #   load(file = paste("results/splevel/best10.type.udi_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_dryseason_", dryseason, ".Rdata", sep = ""), envir = parent.frame(), verbose = FALSE)
  # } else {
  #   load(file = paste("results/commlevel/best10.type.sdi_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_dryseason_", dryseason, ".Rdata", sep = ""), envir = parent.frame(), verbose = FALSE)
  #   load(file = paste("results/commlevel/best10.type.udi_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_dryseason_", dryseason, ".Rdata", sep = ""), envir = parent.frame(), verbose = FALSE)
  # }
  ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ## associating growth time series to their best fit parameters and model outputs and uptake depths----
  ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  # load root params
  load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
  # load(file = paste0("results/", level.folder, "/best10.type.rsq_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, ".Rdata"), envir = parent.frame(), verbose = FALSE)
  # load(file = paste0("results/", level.folder, "/best10.type_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, ".Rdata"), envir = parent.frame(), verbose = FALSE)
  root.param <- info$root.param

  pro.df <- read.csv(file = file.path("results/root.profiles.long.csv"), header = TRUE)

  root.param$root.95 <- as.numeric(by(pro.df, pro.df$rf.sam, function(X) min(X$depth[X$cum.root.frac >= 0.95])))
  root.param$max.root <- as.numeric(by(pro.df, pro.df$rf.sam, function(X) min(X$depth[X$cum.root.frac == 1.0000000])))
  summary(root.param)

  sp.n <- read.csv(file.path(paste("results/sp.n_med_growth_sp_by_size_6.csv")), row.names = 1)
  ### trial
  doParallel::registerDoParallel(ncor)
  ds.bestfit.list <- list()
  ds.bestfit.list <- foreach(ii  = 1 : length(sp_size)) %dopar% {
    rf.sam.rows <- si.param.rel$rf.sam
    ds.bestfit.list[[ii]] <- root.param[rf.sam.rows,] %>% ## this works because rf.sam in root.param is consecutive
      mutate(par.sam = si.param.rel$par.sam,
             rsq = as.numeric(GLUE.rsq.vector[ii,]),
             sdi = as.numeric(GLUE.rsq.vector[ii,]),
             udi = as.numeric(GLUE.rsq.vector[ii,]),
             sp = sp[ii],
             size = size[ii],
             size = factor(size, levels = c("tiny", "small", "medium", "large"))) %>%
      ## reordering levels in Species based on number of trees in each size class
      unite("sp_size", c("sp", "size"), remove = FALSE) %>%
      left_join(select(sp.n %>% mutate(sp_size = as.character(sp_size)), sp_size, n), by = "sp_size") %>%
      subset(!is.na(rf.sam)) %>% droplevels()
  }
  stopImplicitCluster()
  ds.bestfit.all <- do.call(rbind, ds.bestfit.list) %>%
    select(par.sam, sp_size, rsq, udi, root.95, max.root, everything())
  head(ds.bestfit.all)
  if (splevel == "on") {
    ds.bestfit.all$tlplevel <- "sp"
   } else {
    ds.bestfit.all$tlplevel <- "comm"
   }
  save(ds.bestfit.all, file = paste0("results/", level.folder, "/ds.bestfit.all_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_id_dryseason_", dryseason, ".Rdata"))
}


