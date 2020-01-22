## function to calculate rsq between observed growth vs. simulated SI matrix
growth_by_si.func <- function(growth_by_si.info, goodness.fit, statistic, splevel){

  if (splevel == "on"){
    sp.with.tlp <- names(info.2$sp.si)
    growth.meta <- growth_by_si.info$growth.meta %>%
      subset(sp %in% sp.with.tlp)
    sp_size.on  <- growth.meta$sp_size #.full %>% subset(size == "large" || size == "small") %>% select(-size, -sp, -interval.2) %>% as.matrix()
    growth <- growth_by_si.info$growth[sp_size.on]
    sp_size <- names(growth)
    ### to create si for each sp_size, sp.si elements should be repeated as sps in growth
    sp.in.sp_size <- stringr::str_split_fixed(names(growth), "_", 2)[,1]
    sp.si.match <- info.2$sp.si[sp.in.sp_size] ## this repeats sp.si elements (with sp names) in the order of sp.in.sp_size
    ## need to repeat elements of si according to growth[[j]]$interval
    ## since intervals are from 3:7 (8-intervals), and so are corresponding si vec, this would mean removing 3-1 = 2
    location.adj <- as.numeric(strsplit(colnames(sp.si.match[[1]]),".", fixed = TRUE)[[1]][2]) - 1
    sp.si.match <- lapply(sp.si.match, unname)
  } else {
    growth <- growth_by_si.info$growth
    sp_size <- names(growth)
    location.adj <- as.numeric(strsplit(colnames(growth_by_si.info$si),".", fixed = TRUE)[[1]][2]) - 1
    si <- unname(growth_by_si.info$si)
  }
  if (growth_by_si.info$dbh.residuals == "on") {
    growth <- lapply(growth, function (x) {
      x %>% rename(growth = dbh.residuals)
    })
  }
  n.rf.sam <- length(unique(growth_by_si.info$si.param.rel$rf.sam))
  n.par.sam <- length(unique(growth_by_si.info$si.param.rel$par.sam))
  intervals <- info$intervals

  ncor <- detectCores() - 1
  cl <- parallel::makeForkCluster(ncor)
  doParallel::registerDoParallel(cl)

  GLUE <- vector("list", length(sp_size))
  Sys.time()
  beg <- Sys.time()
  GLUE <- foreach(jj = 1:length(sp_size)) %dopar% {
    obs <- growth[[jj]]
    if(splevel == "on") {
      si.matrix <- sp.si.match[[jj]]
    } else {
      si.matrix <- si
    }
    GLUE.sub <- apply(si.matrix, 1, function (x) {
      sim <- x[obs$interval - location.adj]
      cor.val <- cor(obs$growth, sim, method = "pearson", use = "complete.obs")
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
    GLUE[[jj]] <- matrix(unlist(GLUE.sub), nrow = n.par.sam, ncol = n.rf.sam, byrow = FALSE)
  }
  # return(matrix(unlist(GLUE.list, use.names = FALSE), byrow = TRUE, ncol = n.best, nrow = 19))
  Sys.time()
  end <- Sys.time()
  duration <- end - beg
  parallel::stopCluster(cl)
  names(GLUE) <- sp_size
  # 12 min for splevel = "on" sp-size & ~1500 n.rf.sam, 33 min for splevel = "off"
  #GLUE <- data.frame(do.call(rbind, GLUE.list))
  return(GLUE)
}
