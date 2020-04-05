# This file imports the data and sets up the 'Extra' variable

rm(list=ls())
for (i in 1:5) {gc()}
graphics.off()
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, doParallel, foreach, data.table)

GLUEsetup_part1 <- function(current.folder = current.folder, intervals = intervals) {
  # current.folder = "2019-10-14_5000"; intervals = 5
  if(!dir.exists(file.path("results"))) {dir.create(file.path("results"))}
  if(!dir.exists(file.path("results", current.folder))) {dir.create(file.path("results", current.folder))}

  ####-------------------------------------
  #### 1. To enable hydro time series into intervals
  ####-------------------------------------
  ##
  # # Storing census med and interval start and end dates--------------
  # load("data-raw/CTFScensuses/bci.tree1.Rdata")
  # load("data-raw/CTFScensuses/bci.tree2.Rdata")
  # load("data-raw/CTFScensuses/bci.tree3.Rdata")
  # load("data-raw/CTFScensuses/bci.tree4.Rdata")
  # load("data-raw/CTFScensuses/bci.tree5.Rdata")
  # load("data-raw/CTFScensuses/bci.tree6.Rdata")
  # load("data-raw/CTFScensuses/bci.tree7.Rdata")
  # load("data-raw/CTFScensuses/bci.tree8.Rdata")
  # head(bci.tree1)
  #
  # census.years <- c(1982, 1985, 1990, 1995, 2000, 2005, 2010, 2015)
  # bci.tree1$census <- census.years[1]
  # bci.tree2$census <- census.years[2]
  # bci.tree3$census <- census.years[3]
  # bci.tree4$census <- census.years[4]
  # bci.tree5$census <- census.years[5]
  # bci.tree6$census <- census.years[6]
  # bci.tree7$census <- census.years[7]
  # bci.tree8$census <- census.years[8]
  #
  # bci.tree <- rbind.data.frame(bci.tree1, bci.tree2, bci.tree3, bci.tree4, bci.tree5, bci.tree6, bci.tree7, bci.tree8)
  # bci.tree$ExactDate <- as.Date(bci.tree$ExactDate)
  # ggplot(bci.tree, aes(x = ExactDate)) +
  #   geom_histogram() +
  #   facet_wrap(~ census)
  # mediandates <- summarise(bci.tree %>% group_by(census), mediandate = median(ExactDate, na.rm = T))
  # census.meds <- mediandates$mediandate
  # saveRDS(census.meds, "results/census.mediandates.rds") #-------
  ## Use Stored census med dates # ------
  census.meds <- readr::read_rds("results/census.mediandates.rds")
  # can set intervals: to remove first & second census interval; 5; only to remove first, 6
  ##-------------******------------------
  intervals <- intervals # last 5
  ##-------------******------------------
  census <- (length(census.meds) - intervals):length(census.meds) # 3:length(census.meds); if removing 82-85 85-90
  census.meds.chosen <- census.meds[census]
  # [1] "1981-12-27" "1985-05-04" "1990-08-22" "1995-05-18" "2000-05-10" "2005-06-01"
  # [7] "2010-05-17" "2015-07-02"
  census.start <- census.meds[(length(census.meds) - intervals):(length(census.meds) - 1)]
  census.end <- c(census.meds[(length(census.meds) - intervals + 1):length(census.meds) ])
  interval.mean <- data.frame(census.start = census.start, census.end = census.end) %>%
    rowwise %>%
    mutate(interval.mean = mean.Date(c(census.start, census.end)))
  ## adding interval label to each date
  cut.breaks <- census.meds
  cut.labels <- paste("interval", 1:(length(census.meds) -1), sep = ".")
  ## ----------- End of 1 --------------

  ##-------------------------------------
  ## 2. Getting rooting profiles root.frac for root.sam
  ##-------------------------------------
  #root.nsam.long <- read.csv(file = file.path("data", current.folder, "params.root_fraction_by_depth_long.csv"), header = TRUE)
  # Generated from R/3.090_root_parameters_b.R
  pro.df <- read.csv(file = file.path("results/rf.sam_exponentially_decreasing.csv"), header = TRUE)

  root.nsam.long.sub <- pro.df %>% select(rf.sam, depth, root.frac) %>%
    mutate(depth = signif(round(depth, 2), 2))
  # Root.frac does some to 1 for each rf.sam:
  # root.nsam.long.sub %>% group_by(rf.sam) %>% summarise(root.frac = sum(root.frac)) %>% summary()

  root.nsam <- root.nsam.long.sub %>% pivot_wider(names_from = depth, names_prefix = "depth.", values_from = root.frac)
  ## ----------- End of 2 --------------

  ##----------------------------------------------------------------------------
  ## 3. Get Btran by census interval for each best-fit par.sam by root profiles
  ##----------------------------------------------------------------------------
  ##*************************************
  ## At both community and species level tlp
  ##*************************************
  ##--------------
  ## common code
  ##--------------
  nsam <- nrow(root.nsam)
  ncor <- detectCores() - 2
  cl <- parallel::makeForkCluster(ncor)
  doParallel::registerDoParallel(cl)
  rdt <- data.table(root.nsam.long.sub, key = "depth")

  ##*************************************
  ## 3.b  With community level tlp
  ##*************************************
  load(file = file.path("data-raw/psi.rda"))

  traits.indi <- read.csv("data-raw/traits/HydraulicTraits_Kunert/hydraulic_traits_panama_kunert.csv") # Nobby's data
  tlp <- traits.indi %>% group_by(sp) %>% select(-idividual, -ind_ID) %>%
    dplyr::summarise(tlp = mean(mean_TLP_Mpa, na.rm = TRUE))
  write.csv(tlp, file.path("data-raw/traits/HydraulicTraits_Kunert/tlp_sp_mean.csv"), row.names = FALSE)
  n.sp <- length(tlp$sp)
  tlp.mean.positive <- -mean(tlp$tlp, na.rm = TRUE)

  swp.gfac <- psi %>% mutate(psi.positive = -psi) %>%
    mutate(gfac = if_else(psi.positive >= 0 & psi.positive < 0.5, 1,
                          if_else(psi.positive == 0.5, 1,
                                  if_else(psi.positive > tlp.mean.positive, 0,
                                          (tlp.mean.positive - psi.positive)/(tlp.mean.positive - 0.5)))))

  ## depths are: signif(round(depth, 2), 2) in bci.elm.fates.hydro/R/6.0_ELM-FATES_full best-fit_swp.R
  ## so root.frac depths should also be formatted similarly

  hydro.output <- swp.gfac %>% select(date, gfac, par.sam, depth) %>% mutate(date = as.Date(date))
  rm(swp.gfac)

  hydro.output <- hydro.output %>%
    mutate(interval = cut(date, include.lowest = TRUE, breaks = cut.breaks,
                          labels = cut.labels, right = FALSE))
  hydro.output.chosen <- hydro.output %>%
    filter(date %in% min(census.meds.chosen):max(census.meds.chosen)) %>%
    mutate_at(vars(interval), droplevels)

  # set up the 'info' Structure----------
  info <- list(census.meds = census.meds,
               census = census,
               census.meds.chosen = census.meds.chosen,
               census.start = census.start,
               census.end = census.end,
               interval.mean = interval.mean$interval.mean,
               intervals = intervals,
               hydro.output = hydro.output,
               hydro = hydro.output.chosen,
               root.param = root.nsam,
               root.param.long = pro.df,
               si.type = current.folder
  )
  # save this variables to a workspace
  save(file = "results/GLUEsetup_part1_BCI.RData", info)

  ##-------------------------------------------------------------
  ## summarise Btran by census interval for each best-fit par.sam -----
  ##--------------------------------------------------------------
  ## creating a list of nsam with a list of n.best with a matrix to store btran by interval
  sdt <- data.table(hydro.output.chosen, key = "depth")
  ncor <- detectCores() - 1
  cl <- parallel::makeForkCluster(ncor)
  doParallel::registerDoParallel(cl)

  Sys.time()
  beg <- Sys.time()
  drop.months.vec <- list("Feb", "Mar", c("Feb", "Mar"))
  btran.nsam.int <- vector(mode = "list", length = length(drop.months.vec))
  for (k in 1: length(drop.months.vec)){
    sdt.k <- sdt[, ':='(month = format(as.Date(date), "%b"))][!month %chin% drop.months.vec[[k]]]
    ## saving btran by nsam would be too big
    #where Btran=i=1i=zRootFracz * gfac
    btran.nsam.int[[k]] <- foreach(ii  = 1 : nsam) %dopar% {
      hydro <- merge(sdt.k, rdt[rf.sam == ii])
      btran.nsam.int[[ii]] <- hydro[, btran := root.frac*gfac][
        , keyby = .(par.sam, date, interval),
        .(btran.sum = sum(btran, na.rm = TRUE))][ # daily sums should range (0-1)
          , keyby = .(par.sam, interval),
          .(btran = mean(btran.sum, na.rm = TRUE))]
    }
    Sys.time()
    end <- Sys.time()
    (end - beg)
    print(k)
  }
  names(btran.nsam.int) <- lapply(drop.months.vec, paste, collapse = "")
  parallel::stopCluster(cl)
  # 4 min for root.nsam = 599
  # 1.15 hrs for nsam 10K n.best 100 1990-2018
  #
  # , interval := cut(date, include.lowest = TRUE, breaks = cut.breaks,
  #                   labels = cut.labels, right = FALSE)]
  btran.nsam.int[[1]][[nsam]]

  rm(sdt)
  Sys.time()
  beg <- Sys.time()
  n.ensembles <- length(btran.nsam.int[[1]])
  btran.matrix <- si.param.rel <- vector(mode = "list", length = length(drop.months.vec))
  for (k in 1: length(drop.months.vec)){
    btran.mat.list <- lapply(btran.nsam.int[[k]], function(x) {
      btran.wide <- x %>% pivot_wider(names_from = interval, values_from = btran) %>% as.matrix()
      btran.wide[, -1]
    })
    btran.matrix[[k]] <- do.call(rbind, btran.mat.list)
    ### record relational info for si to rf.sam and par.sam
    btran.wide.1 <- btran.nsam.int[[k]][[1]] %>% pivot_wider(names_from = interval, values_from = btran) %>%
      as.matrix()
    n.best <- nrow(btran.wide.1); par.sam.vec <- btran.wide.1[, "par.sam"]
    si.param.rel[[k]] <- data.frame(row.num = 1:nrow(btran.matrix[[k]]),
                                    rf.sam = rep(c(1:n.ensembles), each = n.best),
                                    par.sam = rep(par.sam.vec, times = n.ensembles))
    print(k)
  }
  names(btran.matrix) <- names(si.param.rel) <- lapply(drop.months.vec, paste, collapse = "")

  info.3 <- list(btran = btran.nsam.int,
                 btran.matrix = btran.matrix,
                 si.param.rel = si.param.rel,
                 drop.months.vec = drop.months.vec)
  save(file = "results/GLUEsetup_part1.3_BCI.RData", info.3)

  ##*************************************
  ## 3.b  With species level tlp -----
  ##*************************************
  ## IN Kunert's tlp data: swp ranges from -1.13 to -2.42 MPa
  # at -2.42 MPa f(swp) = 0
  # from 0 to 0.5 f(swp) = 1,
  # f(swp) = linearly decreases from 0 to 1 f(swp) =  1/(2.42-0.5)*2.42 - 1/(2.42-0.5)*swp = (2.42 - swp)/(2.42-0.5)
  # f(swp) =  1.260417 - 0.5208333*swp
  ## such that at -2.42 MPa f(swp) = 0
  deci <- read.csv("data-raw/traits/HydraulicTraits_Kunert/deciduous_species_Meakem.csv")
  deci <- deci %>% mutate(sp = as.character(Species.code), Deciduousness = as.character(Deciduousness)) %>%
    select(sp, Deciduousness)
  tlp.deci <- tlp %>% subset(sp %in% unique(deci$sp)) %>% droplevels()

  sp.vec <- list(tlp$sp, tlp.deci$sp, tlp.deci$sp, tlp.deci$sp)
  tlp.vec <- list(tlp$tlp, tlp.deci$tlp, tlp.deci$tlp, tlp.deci$tlp)

  Sys.time()
  beg <- Sys.time()
  sp.hydro.output.chosen <- vector(mode = "list", length = length(drop.months.vec))
  for (k  in 1 : length(drop.months.vec)) {
    psi.k <- psi %>% mutate(psi.positive = -psi) %>%
      filter(date %in% min(census.meds.chosen):max(census.meds.chosen)) %>%
      mutate(month = format(as.Date(date), "%b")) %>%
      filter(!month %in% drop.months.vec[[k]])
    sp.tlp.vec <- tlp.vec[[k]]
    for (jj  in 1 : length(sp.vec[[k]])) {
      sp.tlp <- -sp.tlp.vec[jj]
      sp.hydro.output.chosen[[k]][[jj]] <- psi.k %>%
        mutate(gfac = if_else(psi.positive >= 0 & psi.positive < 0.5, 1,
                              if_else(psi.positive == 0.5, 1,
                                      if_else(psi.positive > sp.tlp, 0,
                                              (sp.tlp - psi.positive)/(sp.tlp - 0.5))))) %>%
        select(date, gfac, par.sam, depth) %>%
        mutate(date = as.Date(date),
               interval = droplevels(cut(date, include.lowest = TRUE, breaks = cut.breaks,
                                         labels = cut.labels, right = FALSE)))

    }
    names(sp.hydro.output.chosen[[k]]) <- sp.vec[[k]]
    print(k)
  }
  Sys.time()
  end <- Sys.time()
  (end - beg)
  rm(psi); rm(psi.k)
  names(sp.hydro.output.chosen) <- lapply(drop.months.vec, paste, collapse = "")
  save(file = "results/sp.hydro.output.chosen.RData", sp.hydro.output.chosen)
  load("results/sp.hydro.output.chosen.RData")
  # 2.2 min
  # 2 min for root.nsam = 599
  # saving sp.swp.gfac first crashes system as ~20 GB # also parallel too slow

  #   # plot(gfac ~ psi.2, data = sp.swp.gfac[[1]][1:500,], ylab = "Growth Factor", xlab = "Soil Water Potential (-MPa)",
  #   #      main = "SWP - Growth factor Relationship")

  ##-------------------------------------------------------------
  ## summarise Btran by census interval for each best-fit par.sam-----
  ##--------------------------------------------------------------
  cl <- parallel::makeForkCluster(ncor)
  doParallel::registerDoParallel(cl)

  Sys.time()
  beg <- Sys.time()

  sp.btran.matrix <- vector(mode = "list", length = length(drop.months.vec))
  names(sp.btran.matrix) <- lapply(drop.months.vec, paste, collapse = "")
  for (k in 1: length(drop.months.vec)){
    for(jj  in 1 : length(sp.vec[[k]])) {
      sdt <- data.table(sp.hydro.output.chosen[[k]][[jj]], key = "depth")
      sp.btran.nsam.int <- list()
      sp.btran.nsam.int <- foreach(ii  = 1 : nsam) %dopar% {
        hydro <- merge(sdt, rdt[rf.sam == ii])
        sp.btran.nsam.int[[ii]] <- hydro[, btran := root.frac*gfac][
          , keyby = .(par.sam, date, interval),
          .(btran.sum = sum(btran, na.rm = TRUE))][ # daily sums should range (0-1)
            , keyby = .(par.sam, interval),
            .(btran = mean(btran.sum, na.rm = TRUE))]
      }
      sp.btran.mat.list <- lapply(sp.btran.nsam.int, function(x) {
        btran.wide <- x %>% pivot_wider(names_from = interval, values_from = btran) %>%
          as.matrix()
        btran.wide[,-1]
      }
      )
      sp.btran.matrix[[k]][[jj]] <- do.call(rbind, sp.btran.mat.list)
      print(jj)
    }
    names(sp.btran.matrix[[k]]) <- sp.vec[[k]]
  }
  parallel::stopCluster(cl)

  Sys.time()
  end <- Sys.time()
  (end - beg)
  ## 18 hrs for 4 dropped.months rf.sam = 599, 51 sp
  # (end - beg) ## 13.4 sec for 1 sp and nsam = 18, so 5.6 hrs for 51 sp & 533 nsam
  duration <- end-beg
  save(file = "results/duration.RData", duration)
  parallel::stopCluster(cl)
  info.2 <- list(sp.si = sp.btran.matrix)
  save(file = "results/GLUEsetup_part1.2_BCI.RData", info.2)
}


# https://github.com/EcoClimLab/HydraulicTraits/blob/master/data/Panama/processed_trait_data/Panama_all_traits_table_indvidual_level.csv
