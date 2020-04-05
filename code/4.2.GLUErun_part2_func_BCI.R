## function to calculate rsq between observed growth vs. simulated SI matrix
growth_by_si.func <- function(fun.version = fun.version,
                              growth_by_si.info = growth_by_si.info,
                              splevel = splevel,
                              goodness.fit = goodness.fit,
                              drop.months = drop.months,
                              statistic = statistic) {

  if (splevel == "on"){
    sp.with.tlp <- names(info.2$sp.si[[drop.months]])
    growth.meta <- growth_by_si.info$growth.meta %>%
      subset(sp %in% sp.with.tlp)
    sp_size.on  <- growth.meta$sp_size #.full %>% subset(size == "large" || size == "small") %>% select(-size, -sp, -interval.2) %>% as.matrix()
    growth <- growth_by_si.info$growth[sp_size.on]
    sp_size <- names(growth)
    ### to create si for each sp_size, sp.si elements should be repeated as sps in growth
    sp.in.sp_size <- stringr::str_split_fixed(names(growth), "_", 2)[,1]
    ## Taking the sp.si list element for the prescribed drop.months
    sp.si.match <- info.2$sp.si[[drop.months]][sp.in.sp_size] ## this repeats sp.si elements (with sp names) in the order of sp.in.sp_size
    ## need to repeat elements of si according to growth[[j]]$interval
    ## since intervals are from 3:7 (8-intervals), and so are corresponding si vec, this would mean removing 3-1 = 2
    location.adj <- as.numeric(strsplit(colnames(sp.si.match[[1]]),".", fixed = TRUE)[[1]][2]) - 1
    sp.si.match <- lapply(sp.si.match, unname)
  } else {
    growth <- growth_by_si.info$growth
    if (drop.months != "None"){
      # only species that are deciduous need to be retained from the growth matrix
      deci <- read.csv("data-raw/traits/HydraulicTraits_Kunert/deciduous_species_Meakem.csv")
      deci <- deci %>% mutate(sp = as.character(Species.code), Deciduousness = as.character(Deciduousness)) %>%
        select(sp, Deciduousness)
      sp.in.sp_size <- stringr::str_split_fixed(names(growth), "_", 2)[,1]
      growth <- growth[sp.in.sp_size %in% unique(deci$sp)]
    }
    sp_size <- names(growth)
    location.adj <- as.numeric(strsplit(colnames(growth_by_si.info$si[[drop.months]]),".", fixed = TRUE)[[1]][2]) - 1
    ## Taking the sp.si list element for the prescribed drop.months
    si <- unname(growth_by_si.info$si[[drop.months]])
  }
  # if (growth_by_si.info$dbh.residuals == "on") {
  #   growth <- lapply(growth, function (x) {
  #     x %>% rename(growth = dbh.residuals)
  #   })
  # }
  n.rf.sam <- length(unique(growth_by_si.info$si.param.rel$rf.sam))
  n.par.sam <- length(unique(growth_by_si.info$si.param.rel$par.sam))
  intervals <- info$intervals

  ncor <- detectCores() - 1
  cl <- parallel::makeForkCluster(ncor)
  doParallel::registerDoParallel(cl)

  GLUE.list <- vector("list", length(sp_size))
  Sys.time()
  beg <- Sys.time()
  if (fun.version == "rsq") {
    GLUE.list <- foreach(jj = 1:length(sp_size)) %dopar% {
      obs <- growth[[jj]]
      if(splevel == "on") {
        si.matrix <- sp.si.match[[jj]]
      } else {
        si.matrix <- si
      }
      GLUE.sub <- apply(si.matrix, 1, function (x) {
        sim <- x[obs$interval - location.adj]
        cor.val <- cor(obs$median, sim, method = "pearson", use = "complete.obs")
        if (!is.na(cor.val) & (abs(cor.val) > (goodness.fit - 0.01))) {
          if (statistic == "rsq") {
            value <- cor.val^2 } else {
              value <- cor.val
            }
        } else {
          value <- NA
        }
        # list(cor = cor, rsq = rsq)
        # model <- lm(y ~ obs)
        # summary(model)$r.squared
        # list(model = model, rsq = rsq)
      }
      )
      GLUE.list[[jj]] <- unlist(GLUE.sub)
    }
  } else if (fun.version == "likelihood") {
    GLUE.list <- foreach(jj = 1:length(sp_size)) %dopar% {
      obs <- growth[[jj]]
      if(splevel == "on") {
        si.matrix <- sp.si.match[[jj]]
      } else {
        si.matrix <- si
      }
      GLUE.sub <- apply(si.matrix, 1, function (x) {
        sim <- x[obs$interval - location.adj]
        model <- lm.fit(matrix(sim), matrix(obs$median))
        within.ci <- data.table::between(model$fitted.values, obs$lwr, obs$upr)
        # return(all(data.table::between(model$fitted.values, obs$lwr, obs$upr), na.rm = FALSE))
        # returns no. of data points(censuses) out of the total data points(censuses) in a time series
        # that fell within the CI
        return(sum(within.ci, na.rm = TRUE))
      })
      GLUE.list[[jj]] <- unlist(GLUE.sub)
    }
  }
  Sys.time()
  end <- Sys.time()
  duration <- end - beg
  parallel::stopCluster(cl)
  GLUE <- data.frame(do.call(rbind, GLUE.list))
  row.names(GLUE) <- sp_size
  # 12 min for splevel = "on" sp-size & ~1500 n.rf.sam, 33 min for splevel = "off"
  return(GLUE)
}
