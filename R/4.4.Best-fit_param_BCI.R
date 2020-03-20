# ## Not anymore: Finding best fit water-stress/growth model parameters for species growth time series
# ## Just vectorises Rsquares with decreasing order and notes down the order to relate to par.sam and rf.sam from growth_by_si.info$si.param.rel
#
# # for 50 ha obs species groups
# rm(list = ls())
# gc()
# if (!require("pacman")) install.packages("pacman"); library(pacman)
# pacman::p_load(tidyverse, doParallel, foreach, data.table)
#
# bestfit.param <- function(splevel) {
#   # load interval and n.ensembles
#   load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
#   load(file.path("results/4.1GLUEsetup_part2_BCI.RData")) # has n.ensembles and growth and si matrix
#
#   intervals <- info$intervals
#   n.ensembles <- growth_by_si.info$n.ensembles
#   growth.type <- growth_by_si.info$growth.type ## this is not data but type such as "individual"
#   growth.selection <- growth_by_si.info$growth.selection
#   si.type <- growth_by_si.info$si.type
#   dbh.residuals <- growth_by_si.info$dbh.residuals
#   goodness.fit <- 0.3 # rsq0.3
#   if (splevel == "on") {
#     level.folder <- "splevel"
#   } else {
#     level.folder <- "commlevel"
#   }
#   ncor <- detectCores() - 2
#
#   if (splevel == "on") {
#     load(paste0("results/sp.GLUE.rsq_", si.type, "_", n.ensembles, "_",
#               growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_cor", goodness.fit, ".Rdata"))
#     GLUE.rsq <- sp.GLUE.rsq
#     rm(sp.GLUE.rsq)
#   } else {
#     load(paste0("results/GLUE.rsq_", si.type, "_", n.ensembles, "_",
#                 growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, "_cor", goodness.fit, ".Rdata"))
#   }
#   sp_size <- names(GLUE.rsq)
#   # growth_by_si.info$growth only has growth matrix for given info$intervals
#   # to get corresponding species-size, use growth meta
#
#   registerDoParallel(ncor)
#   GLUE.rsq.vector.list <- list()
#   GLUE.rsq.vector.list <-
#     foreach(ii  = 1:length(GLUE.rsq)) %dopar% {
#       mat <- GLUE.rsq[[ii]]
#       GLUE.rsq.vector.list[[ii]] <-
#         as.vector(mat)
#     }
#
#   GLUE.rsq.vector <- data.frame(do.call(rbind, GLUE.rsq.vector.list))
#   row.names(GLUE.rsq.vector) <- sp_size
#   save(GLUE.rsq.vector, file = paste0("results/", level.folder, "/GLUE.rsq.vector_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection, "_", dbh.residuals, "_", intervals, ".Rdata"))
# }
#
#
