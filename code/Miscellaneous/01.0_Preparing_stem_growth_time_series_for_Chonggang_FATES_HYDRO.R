#---------------------------------
# Title: Preparing growth data for inverse modeling
# Author : Rutuja Chitra-Tarak
# Original date: December 11, 2019
#---------------------------------

# for 50 ha obs species groups
rm(list = ls())
gc()
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, purrrlyr, scales, grid, gridExtra, lme4)

# graphics info
theme_set(theme_bw())
theme_update(text = element_text(size = 14),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)
# time series info
years <- seq(1990, 2015, by = 5)
ncensus  <- length(years)
nint <- length(years) - 1
intervals <- nint
# sp4 <- c("ANOL", "TECG", "LAGL", "TERT")
##------------------------------------------------
## Collating tree.full data for all censuses
##------------------------------------------------
# https://repository.si.edu/handle/10088/20925
# load("data-raw/CTFScensuses/BCI.stem1.Rdata") # 1982-1985
load("data-raw/CTFScensuses/BCI.stem2.Rdata") # 1985-1990
load("data-raw/CTFScensuses/BCI.stem3.Rdata")
load("data-raw/CTFScensuses/BCI.stem4.Rdata")
load("data-raw/CTFScensuses/BCI.stem5.Rdata")
load("data-raw/CTFScensuses/BCI.stem6.Rdata")
load("data-raw/CTFScensuses/BCI.stem7.Rdata")
load("data-raw/CTFScensuses/BCI.stem8.Rdata")

length(unique(bci.stem3$sp))
# all censuses in one list
stem.full <- bind_rows(bci.stem2, bci.stem3, bci.stem4, bci.stem5, bci.stem6, bci.stem7, bci.stem8)
# > nrow(stem.full)
# [1] 4740360
# > nrow(stem.full)*6
# [1] 28442160
###------------selection----------
census <- 2:8
###-------------------------------
stem.full$census <- rep(census, each = nrow(bci.stem3))
## remove duplicated tags if any
nrow(stem.full) # 4740360
stem.full <- stem.full %>% mutate(ExactDate = as.Date(ExactDate))
stem.full <- stem.full %>% group_by(census) %>% distinct(stemID, .keep_all= TRUE)
nrow(stem.full)
# 4740360
# so no duplicate stemIDs;
n_stem.full <- stem.full %>% group_by(census) %>% nest()
map_dbl(n_stem.full$data, nrow)
# 790060 790060 790060 790060 790060 790060 790060

##------------------------------------------------
## to remove stems for which hom have changed--------
##------------------------------------------------
# pool together for each stem.full hom in censuses
hom.mat <- map_dfc(n_stem.full$data, "hom")
hom.mat <- data.frame(invisible(lapply(hom.mat, function (x) {as.numeric(x)})))
# check by row (stem.full) if hom value remains the same
hom.same <- apply(hom.mat, 1, function (x) {diff(range(x, na.rm = T)) == 0})
# xx <- dbh.mat %>%
#   rowwise() %>% mutate(diff = function (x) {diff(range(x, na.rm = T)) == 0})
## if hom has changed, resut is FALSE
length(which(hom.same == FALSE))
# 185969
## those with all NA rows are false positively identified as having hom changed, so making those TRUE
hom.same[rowSums(is.na(hom.mat)) == length(census)] <- TRUE
length(which(hom.same == FALSE))
## 4530

## stem.fulls with unchanged hom
keep.stemIDs <- n_stem.full$data[[1]]$stemID[hom.same]
length(keep.stemIDs) #785530
# only selecting those stem.fulls out of
map_dbl(n_stem.full$data, nrow)
# [1] 790060 790060 790060 790060 790060 790060
n_stem <- stem.full %>% filter(stemID %in% keep.stemIDs) %>% group_by(census) %>% nest()
map_dbl(n_stem$data, nrow)
# 785530 785530 785530 785530 785530 785530
# stems removed
map_dbl(n_stem.full$data, nrow) - map_dbl(n_stem$data, nrow)
# 4530 4530 4530 4530 4530 4530
##--------------------End of removing stems with hom change----------------------------


##------------------------------------------------
# Getting a matrix of date by census for these stems
##------------------------------------------------
time.mat <- map_dfc(n_stem$data, "ExactDate")
n_date.mat <- map_dfc(n_stem$data, "date")
colnames(time.mat) <- paste("census", census, sep = ".")
# getting a matrix of dbh by census for these stems
dbh.mat <- map_dfc(n_stem$data, "dbh")
colnames(dbh.mat) <- paste("census", census, sep = ".")

dbh.mat[1:5,]
summary(dbh.mat)
## Removing extreme negative measurement errors in dbh-----
# Condit et al 2017 page 4: if a stem's later dbh measurements fell < 4sd_e (sd_e = 0.006214*dbh + 0.9036) compared to earlier dbh, then an extreme error
dbh.diff.temp <- t(diff(t(dbh.mat))) # in mm
# first note those records and assign them NAs## get growth
dbh.lim <- 4*(0.006214*dbh.mat[,-ncol(dbh.mat)] + 0.903)
## negative of negative growth is larger than allowed error
low.out <- -dbh.diff.temp > dbh.lim
length(which(as.vector(low.out) == TRUE))
# 3338
dbh_low.out <- data.frame(dbh.mat)
## make the identified erroneous second dbh NA
for (i in 1:ncol(low.out)) {dbh_low.out[which(low.out[, i]), i + 1] <- NA}
head(which(low.out[,1] == "TRUE"))
dbh.mat[which(low.out[,1] == "TRUE")[1:20],]
dbh_low.out[which(low.out[,1] == "TRUE")[1:20],]
## This matrix (same dimensions as dbh) needs to be used in growth calculation
##--------------------End-------------------------


##------------------------------------------------
## Now to obtain growth rates-----
##------------------------------------------------
## subtracting successive columns
time.diff <- t(diff(t(n_date.mat)))/365 # from days to yr
colnames(time.diff) <- paste("interval", census[-length(census)], sep = ".")
dbh.diff <- t(diff(t(dbh_low.out))) # in mm
dbh.diff[1:100,]
g <- dbh.diff/time.diff # mm/yr
mean(g, na.rm = T); median(g, na.rm = T)
colnames(g) <- paste("interval", census[-length(census)], sep = ".")

g <- data.frame(g)
summary(g)
# [1] 0.8045554 mean
# [1] 0.394808 median
## not the same

## Removing extreme positive errors in growth-----
# removing positive errors > 75 mm/yr since that's the growth of the fastest growing stems (true growth, from a double blinded test)
## note extreme records and assign them NAs
g_up.out <- apply(g, 2, function (x) {ifelse(x > 75, NA, x)})
mean(g_up.out, na.rm = T); median(g_up.out, na.rm = T)
# [1] 0.8034973
# [1] 0.394808
## This matrix needs to be used in further growth calculation
# saving growth without modulus transformation, after all outliers were removed
# save(g_up.out, file = "results/stem_growth_without_outliers.Rdata")
##------------------End------------------------------
# jpeg("figures/dbh/stem_growth_without_outliers.jpeg")
# plot(x = as.vector(as.matrix(dbh_low.out[, - ncol(dbh_low.out)])), y = as.vector(as.matrix(g_up.out)),
#      ylab = "Growth rate (mm/yr)", xlab = "DBH (mm)")
# graphics.off()

#
# ##--------------------------------------------------
# ## Saving individual stem data for HADAD with following selection criterion:
# ## stems with at least 3 records
# ##--------------------------------------------------
# intervals <- nint; #census.cols = (ncol(gt_back) - intervals + 1):ncol(gt_back)
# rows.3obs <- unlist(apply(gt_back, 1, function(x) {sum(!is.na(x)) >= 3}))
# ###---------------
# sum(rows.3obs); #1,47,913
# sum(rows.3obs)/nrow(gt_back) # 18% of stems have at least 3 records
# gt_back_for_hadad <- gt_back[rows.3obs,]
#
# save(gt_back_for_hadad, file =  paste0("results/growth_individual_stem_intervals_", intervals ,".Rdata"))
# # Now saving stem data for these selected stems, dbh for all censuses are saved
# stem.data <- n_stem$data[[1]][rows.3obs, 1:9]
#
# cut.breaks <- c(10, 50, 100, 300, max(stem.data$dbh, na.rm = T)) # about 45 trees per interval are > 1000 mm in dbh
# cut.labels <- c("tiny", "small", "medium", "large")
# # old cut.breaks <- c(10, 50, 60, 80, 100, 150, 300, 700)
# # cut.labels <- c("tiny", "out", "small", "out", "medium", "out", "large")
# stem.data <- stem.data %>% mutate(size = cut(dbh,
#                                include.lowest = TRUE, breaks = cut.breaks,
#                                labels = cut.labels, right = TRUE))
# ## about 2000 gaps in size, because dbh absent, so taking dbh from other censuses
# #29023
# length(stem.data$size[is.na(stem.data$size)])
# for (i in 1: intervals){
#   dbh.i <- n_stem$data[[1 + i]]$dbh[rows.3obs]
#   size.i <- cut(dbh.i, include.lowest = TRUE, breaks = cut.breaks,
#                 labels = cut.labels, right = TRUE)
#   stem.data <- stem.data %>%
#     mutate(size = if_else(is.na(size), size.i, size))
# }
# #length(stem.data$size[is.na(stem.data$size)])
# stem.data <- stem.data %>% unite("sp_size", c("sp", "size"), remove = FALSE)
# save(stem.data, file =  paste0("results/growth_individual_stem.data_intervals_", intervals ,".Rdata"))
# load(file = paste0("results/growth_individual_stem.data_intervals_", intervals ,".Rdata"))
# ##------------------End------------------------------

cut.breaks <- c(c(0, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100)*10, max(dbh_low.out, na.rm = T)) # about 45 trees per interval are > 1000 mm in dbh
# cut.labels <- c("tiny", "small", "medium", "large")
# old cut.breaks <- c(10, 50, 60, 80, 100, 150, 300, 700)
# cut.labels <- c("tiny", "out", "small", "out", "medium", "out", "large")
size.class <- data.frame(apply(dbh_low.out[, -ncol(dbh_low.out)], 2, cut,
                               include.lowest = TRUE, breaks = cut.breaks*10, right = TRUE))
colnames(size.class) <- paste("interval", census[-length(census)], sep = ".")

## gro and data for one census in n_tree (one list element) have same stem identities by row, so species id is added from n_tree
gro.wide <- g_up.out %>% data.frame() %>% mutate(stemID = n_stem$data[[1]]$stemID) %>% mutate(sp = n_stem$data[[1]]$sp)
## converting growth to long format and adding size class identity to each row/stem

gro.long <- data.frame(gro.wide) %>%
  pivot_longer(cols = starts_with("interval."), names_to = "interval",
               names_prefix = "interval.", values_to = "growth") %>%
##--------------------------------------------------
## a stem's size class can change across censuses
## and thus which size class group it belongs to, whose mean is being taken
##--------------------------------------------------
  mutate(size = as.vector(as.matrix(size.class)),
         dbh = as.vector(as.matrix(dbh_low.out[, -ncol(dbh_low.out)])),
         interval.num = as.numeric(interval),
         interval = recode(interval, `1` = "1982-85", `2` = "1985-90", `3` = "1990-95",
                              `4` = "1995-00", `5` = "2000-05", `6` = "2005-10", `7` = "2010-15")) %>%
  # remove obs outside the defined classes
  filter(!is.na(size) & !is.na(growth)) %>%
  mutate(sp_size = paste(sp, size, sep = "_")) %>%
  arrange(sp_size, interval)

sp_size.stem <- split(gro.long %>% select(sp, size, interval.num, interval, growth), gro.long$sp_size)

sp_size.stats.log <- lapply(sp_size.stem, function(x){
  x %>% group_by(sp, size, interval, interval.num) %>%
    summarise(median = median(log(growth), na.rm = TRUE),
              mean = mean(log(growth), na.rm = TRUE),
              sd.2 = 2*sd(log(growth), na.rm = TRUE)) %>% arrange(interval)
}
)
sp_size.stats <- lapply(sp_size.stem, function(x){
  x %>% group_by(interval, interval.num) %>%
    summarise(median = median(growth, na.rm = TRUE),
              mean = mean(growth, na.rm = TRUE),
              sd.2 = 2*sd(growth, na.rm = TRUE)) %>% arrange(interval)
}
)
## removing sp_size for which record present only for one interval
# sp_size.stem.n <- as.numeric(unlist(lapply(sp_size.stem, function(x){length(x$growth[!is.infinite(x$growth)])})))
# sp_size.stem <- sp_size.stem[-which(sp_size.stem.n %in% c(0,1))]
# sp_size.med.n <-  as.numeric(unlist(lapply(sp_size.med, function(x){length(x$growth[!is.infinite(x$growth)])})))
# sp_size.med <- sp_size.med[-which(sp_size.med.n %in% c(0,1))]
# sp_size.mean.n <-  as.numeric(unlist(lapply(sp_size.mean, function(x){length(x$growth[!is.infinite(x$growth)])})))
# sp_size.mean <- sp_size.mean[-which(sp_size.mean.n %in% c(0,1))]
#
# sp_size.2sd.n <-  as.numeric(unlist(lapply(sp_size.2sd, function(x){length(x$growth[!is.infinite(x$growth)])})))
# sp_size.2sd <- sp_size.2sd[-which(sp_size.2sd.n %in% c(0,1))]
#
# sp_size.stem.names <- names(sp_size.stem) # 889
# sp_size.mean.names <- names(sp_size.mean)
# sp_size.med.names <- names(sp_size.med) # 869
growth.selection <- "growth_refined_size_classes_varying_non_cc_for_Chonggang"
# save(sp_size.stem, file = paste0("results/sp_size.individual_growth_dbh.residuals_off_", intervals, "_", growth.selection, ".Rdata"))
save(sp_size.stats, file = paste0("results/sp_size.stats_growth_dbh.residuals_off_", intervals, "_", growth.selection, ".Rdata"))

growth.selection <- "log-transformed_growth_refined_size_classes_varying_non_cc_for_Chonggang"

save(sp_size.stats.log, file = paste0("results/sp_size.stats.log_growth_dbh.residuals_off_", intervals, "_", growth.selection, ".Rdata"))

