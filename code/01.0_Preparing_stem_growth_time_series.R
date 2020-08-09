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

range01 <- function(x){(x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))}

# sp4 <- c("ANOL", "TECG", "LAGL", "TERT")
##------------------------------------------------
## Collating tree.full data for all censuses
##------------------------------------------------
# https://repository.si.edu/handle/10088/20925
# load("data-raw/CTFScensuses/BCI.stem1.Rdata") # 1982-1985
# load("data-raw/CTFScensuses/BCI.stem2.Rdata") # 1985-1990
load("data-raw/CTFScensuses/BCI.stem3.Rdata")
load("data-raw/CTFScensuses/BCI.stem4.Rdata")
load("data-raw/CTFScensuses/BCI.stem5.Rdata")
load("data-raw/CTFScensuses/BCI.stem6.Rdata")
load("data-raw/CTFScensuses/BCI.stem7.Rdata")
load("data-raw/CTFScensuses/BCI.stem8.Rdata")

length(unique(bci.stem3$sp))
# all censuses in one list
stem.full <- bind_rows(bci.stem3, bci.stem4, bci.stem5, bci.stem6, bci.stem7, bci.stem8)
# > nrow(stem.full)
# [1] 4740360
# > nrow(stem.full)*6
# [1] 28442160
###------------selection----------
census <- 3:8
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
max(g_up.out, na.rm = T) # 69.352
## This matrix needs to be used in further growth calculation
# saving growth without modulus transformation, after all outliers were removed
save(g_up.out, file = "results/stem_growth_without_outliers.Rdata")
##------------------End------------------------------
jpeg("figures/dbh/stem_growth_without_outliers.jpeg")
plot(x = as.vector(as.matrix(dbh_low.out[, - ncol(dbh_low.out)])), y = as.vector(as.matrix(g_up.out)),
     ylab = expression("Growth rate (mmyr"^-1*")"), xlab = "DBH (mm)")
graphics.off()

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

cut.breaks <- c(10, 50, 100, 300, max(dbh_low.out, na.rm = T)) # about 45 trees per interval are > 1000 mm in dbh
cut.labels <- c("tiny", "small", "medium", "large")
# old cut.breaks <- c(10, 50, 60, 80, 100, 150, 300, 700)
# cut.labels <- c("tiny", "out", "small", "out", "medium", "out", "large")
size.class <- data.frame(apply(dbh_low.out[, -ncol(dbh_low.out)], 2, cut,
                               include.lowest = TRUE, breaks = cut.breaks,
                               labels = cut.labels, right = TRUE))
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
         interval = as.numeric(interval)) %>%
  # remove obs outside the defined classes
  filter(!is.na(size) & !is.na(growth)) %>%
  mutate(sp_size = paste(sp, size, sep = "_")) %>%
  arrange(sp_size, interval)

sp_size.stem <- split(gro.long %>% select(interval, growth), gro.long$sp_size)
## taking summary stat for each sp-size by census interval (stem's size identity varying with interval)
sp_size.med <- lapply(sp_size.stem, function(x){
  x %>% group_by(interval) %>%
    summarise(growth = median(growth, na.rm = TRUE), .groups = "drop_last") %>% arrange(interval)
  }
  )
sp_size.mean <- lapply(sp_size.stem, function(x){
  x %>% group_by(interval) %>%
    summarise(growth = mean(growth, na.rm = TRUE), .groups = "drop_last") %>% arrange(interval)
}
)

## removing sp_size for which record present only for one interval
sp_size.stem.n <- as.numeric(unlist(lapply(sp_size.stem, function(x){nrow(x)})))
sp_size.stem <- sp_size.stem[-which(sp_size.stem.n == 1)]
sp_size.med.n <-  as.numeric(unlist(lapply(sp_size.med, function(x){nrow(x)})))
sp_size.med <- sp_size.med[-which(sp_size.med.n == 1)]
sp_size.mean.n <-  as.numeric(unlist(lapply(sp_size.mean, function(x){nrow(x)})))
sp_size.mean <- sp_size.mean[-which(sp_size.mean.n == 1)]

sp_size.stem.names <- names(sp_size.stem) # 889
sp_size.mean.names <- names(sp_size.mean)
sp_size.med.names <- names(sp_size.med) # 869
growth.selection <- "size_class_varying_non_cc"
save(sp_size.stem, file = paste0("results/sp_size.individual_growth_dbh.residuals_off_", intervals, "_", growth.selection, ".Rdata"))
save(sp_size.med, file = paste0("results/sp_size.med_growth_dbh.residuals_off_", intervals, "_", growth.selection, ".Rdata"))
save(sp_size.mean, file = paste0("results/sp_size.mean_growth_dbh.residuals_off_", intervals, "_", growth.selection, ".Rdata"))
save(sp_size.stem.names, file = paste0("results/sp_size.individual.names_", intervals,"_", growth.selection, ".Rdata"))
save(sp_size.med.names, file = paste0("results/sp_size.med.names_", intervals, "_", growth.selection, ".Rdata"))
save(sp_size.mean.names, file = paste0("results/sp_size.mean.names_", intervals, "_", growth.selection, ".Rdata"))
save(gro.long, file = paste0("results/gro.long_", intervals, "_", growth.selection, ".Rdata"))

## summary stats across all census intervals
sp_size.stem.all.cen <- split(gro.long %>% select(growth), gro.long$sp_size)
sp_size.stats.all.cen <- lapply(sp_size.stem.all.cen, function(x){
  x %>% droplevels() %>%
    summarise(median = median(growth, na.rm = TRUE),
                  mean = mean(growth, na.rm = TRUE),
                  n = length(growth),
                  sd = sd(growth, na.rm = TRUE),
                  se = sd/sqrt(n),
                  trees = length(unique(x$stemID)),
                  .groups = "drop_last")
}
)
save(sp_size.stats.all.cen, file = paste0("results/sp_size.stats_growth_dbh.residuals_off_", intervals, "_", growth.selection, ".Rdata"))

sp_size.stats.all.cen.df <- gro.long %>% group_by(sp_size) %>%
  summarise(median = median(growth, na.rm = TRUE),
            mean = mean(growth, na.rm = TRUE),
            n = length(growth),
            sd = sd(growth, na.rm = TRUE),
            se = sd/sqrt(n),
            trees = length(unique(stemID)),
            .groups = "drop_last")
large.stats.all.cen.df <- gro.long %>%
  subset(size %in% c("large")) %>%
  group_by(sp) %>%
  summarise(median = median(growth, na.rm = TRUE),
            mean = mean(growth, na.rm = TRUE),
            n = length(growth),
            sd = sd(growth, na.rm = TRUE),
            se = sd/sqrt(n),
            trees = length(unique(stemID)),
            .groups = "drop_last")
adult.stats.all.cen.df <- gro.long %>%
  subset(size %in% c("medium", "large")) %>%
  group_by(sp) %>%
  summarise(median = median(growth, na.rm = TRUE),
            mean = mean(growth, na.rm = TRUE),
            n = length(growth),
            sd = sd(growth, na.rm = TRUE),
            se = sd/sqrt(n),
            trees = length(unique(stemID)),
            .groups = "drop_last")
juvenile.stats.all.cen.df <- gro.long %>%
  subset(size %in% c("tiny", "small")) %>%
  group_by(sp) %>%
  summarise(median = median(growth, na.rm = TRUE),
            mean = mean(growth, na.rm = TRUE),
            n = length(growth),
            sd = sd(growth, na.rm = TRUE),
            se = sd/sqrt(n),
            trees = length(unique(stemID)),
            .groups = "drop_last")
save(sp_size.stats.all.cen.df, file = paste0("results/sp_size.stats_growth_dbh.residuals_off_", intervals, "_", growth.selection, ".df.Rdata"))
save(large.stats.all.cen.df, file = paste0("results/large.stats_growth_dbh.residuals_off_", intervals, "_", growth.selection, ".df.Rdata"))
save(adult.stats.all.cen.df, file = paste0("results/adult.stats_growth_dbh.residuals_off_", intervals, "_", growth.selection, ".df.Rdata"))
save(juvenile.stats.all.cen.df, file = paste0("results/juvenile.stats_growth_dbh.residuals_off_", intervals, "_", growth.selection, ".df.Rdata"))

##-----------------------
### getting dbh residuals
##-----------------------
### power law (growth = r*dbh^d) as in literature/Ecol Lett 2006 Muller-Landau.pdf
## hardly any data by species and it's wild (not enough to get mean by size class bind either)
## So using a mixed effects model

## growth = r*dbh^d
## log(growth) = log(r*dbh^d)) = log(r) + d*log(dbh) ## the linear model
## growth = exp(log(growth))

#### doing the same thing with size bins
cut.breaks.1 <- 10*c(1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.1, 4.4, 4.7, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 41, 44, 47, 50, 55, 60, 65, 70, 75, 80, 85 ,90, 95, 100, 1000)
gro.long.mod <- gro.long %>%
  mutate(size.bins.1 = cut(dbh, include.lowest = TRUE, breaks = cut.breaks.1,
                           right = FALSE))
level_key <- c( `3` = "1990-1995", `4` = "1995-2000", `5` = "2000-2005",  `6` = "2005-2010",  `7` = "2010-2015")
gro.long.mod.med.allsp <- gro.long.mod %>%
  group_by(interval, size.bins.1) %>%
  summarise(growth.mean.allsp = mean(growth, na.rm = TRUE),
            growth.median.allsp = median(growth, na.rm = TRUE),
         dbh.mean.allsp = mean(dbh, na.rm = TRUE),
         .groups = "drop_last") %>%
  ungroup(interval, size.bins.1) %>%
  arrange(size.bins.1, interval) %>%
  mutate(interval = recode_factor(as.factor(interval), !!!level_key))
View(gro.long.mod.med.allsp)

####

lm.model <- lm(log(growth.mean.allsp) ~ log(dbh.mean.allsp) , data = gro.long.mod.med.allsp %>% subset(growth.mean.allsp != 0))
glm.model <- glm(growth.mean.allsp ~ splines::bs(dbh.mean.allsp, 5), data = gro.long.mod.med.allsp)
summary(lm.model)
### using log(growth.mean.allsp) --------
lm.model$coefficients
# (Intercept) log(dbh.mean.allsp)
# -0.6528035           0.1397581
gro.long.mod.med.allsp <- gro.long.mod.med.allsp %>%
  subset(growth.mean.allsp != 0) %>%
  mutate(lm.predict.growth.allsp = exp(predict(lm.model, newdata = gro.long.mod.med.allsp)),
           # exp(as.numeric(lm.model.1$coefficients[1]))*
           # dbh.mean.allsp^as.numeric(lm.model.1$coefficients[2]),
         glm.predict.growth.allsp = predict(glm.model, newdata = gro.long.mod.med.allsp))

pred.data <- data.frame(dbh.mean.allsp = seq(10, 1300, by = 1))
pred.data <- pred.data %>%
  mutate(lm.predict.growth.allsp = exp(predict(lm.model, newdata = pred.data)),
         glm.predict.growth.allsp = predict(glm.model, newdata = pred.data))
g1 <- ggplot(gro.long.mod.med.allsp,
             aes(x = dbh.mean.allsp/10, y = growth.mean.allsp/10)) +
  geom_point(aes(color = as.factor(interval)), alpha = 0.9) +
  scale_color_discrete(name = "Interval") + theme(legend.position = c(0.2, 0.25)) +
  scale_y_log10() + scale_x_log10(breaks = c(1, 2, 5, 10, 20, 30, 50, 100)) +
  # ggtitle("All data") +
  ylab(expression("Growth rate (cmyr"^-1*")"))  + xlab("DBH (cm)")
g2 <- g1 +
  geom_line(data = pred.data, aes(y = glm.predict.growth.allsp/10), color = "black")
ggsave("growthrate_dbh_predicted_growth_with_spline_degree_5_dbh_cutoff_100.jpeg", plot = g2, path =
         file.path("figures/dbh/"), height = 4.5, width = 4.5, units='in')
g3 <- g1 +  geom_line(data = pred.data, aes(y = lm.predict.growth.allsp/10), color = "black")
ggsave("growthrate_dbh_predicted_growth_with_lm.jpeg", plot = g3, path =
         file.path("figures/dbh/"), height = 4.5, width = 4.5, units='in')
# to show with equation
formula.glm <- y ~ poly(x, 5, raw = TRUE)
g4 <- g1 +
  stat_poly_eq(aes(label = paste(stat(eq.label))),
               npcx = 0.05, npcy = 0.97, rr.digits = 2,
               formula = formula.glm, parse = TRUE, size = 3) +
  stat_poly_eq(aes(label = paste(stat(adj.rr.label))),
               npcx = 0.05, npcy = 0.90, rr.digits = 2,
               formula = formula.glm, parse = TRUE, size = 3) +
  geom_smooth(method = glm, formula = y ~ splines::bs(x, 5), color = "black", size = 0.7)
ggsave("growthrate_dbh_predicted_growth_with_spline_degree_5_dbh_cutoff_100_model.jpeg", plot = g4, path =
         file.path("figures/dbh/"), height = 4.5, width = 4.5, units='in')

## not enough species-wise data to fit a trend...often giving rise to negative trends
## So using community wide fit
gro.long.mod <- gro.long.mod %>%
  mutate(dbh.predict.growth = predict(glm.model, newdata = data.frame(dbh.mean.allsp = gro.long.mod$dbh)),
         dbh.residuals = growth - dbh.predict.growth)
####
## plotting the predicted growth
sp.stem.dbh.full <- split(gro.long.mod %>% select(interval, growth, dbh.predict.growth, dbh, dbh.residuals), gro.long.mod$sp)

## finding out sample size per size class sample size
sample.size <- 3 # at least those number of stems to calculate the power law or mean growth/residuals
sp.n.list <- lapply(sp.stem.dbh.full, function(x){
  data.frame(n = nrow(x))})
sp.n  <- rbindlist(sp.n.list, idcol = "sp_size")
select.sp <- sp.n$sp[sp.n$n >= sample.size]
sp.stem.dbh <- sp.stem.dbh.full[select.sp]

pdf("figures/dbh/growth_dbh_splines_plot_by_species.pdf", onefile=TRUE)
for (i in 1:length(sp.stem.dbh)) {
  df <- sp.stem.dbh[[i]]
  plot(x = df$dbh*0.1, y = df$growth*0.1, log = "xy", main = names(sp.stem.dbh[i]),
       ylab = "Growth rate (cm/yr)", xlab = "DBH (cm)")
  #ylim = c(0, max(gro.long.mod.join$growth*0.1)), xlim = c(1, max(gro.long.mod.join$dbh*0.1)))
  points(x = df$dbh*0.1, y = df$dbh.predict.growth*0.1, col = "red")
}
graphics.off()

pdf("figures/dbh/residuals.growth_dbh_splines_plot_by_species.pdf", onefile=TRUE)
for (i in 1:length(sp.stem.dbh)) {
  df <- sp.stem.dbh[[i]]
  plot(x = df$dbh*0.1, y = df$dbh.residuals*0.1, main = names(sp.stem.dbh[i]),
       ylab = "DBH Residual Growth rate (cm/yr)", xlab = "DBH (cm)")
  points(x = df$dbh*0.1, y = df$dbh.predict.growth*0.1, col = "red")
}
graphics.off()

## not enough species-wise data to fit a trend...often giving rise to negative trends--------
#
# mixed.model.1 <- lmer(log(growth.med.sp) ~ log(dbh.mean.sp) + (1 | sp), data = gro.long.mod %>% subset(growth.med.sp != 0))
# summary(mixed.model.1)
# ### using raw data log(growth)
# ### using log(growth.mean.sp)
# fixef(mixed.model.1)
# # (Intercept) log(dbh.mean.sp)
# # -0.4331291        0.1016471
#
# ranef.sp.1 <- data.frame(ranef(mixed.model.1)) %>%
#   select(grp, condval) %>%
#   rename(sp = grp, intercept = condval) %>%
#   mutate(sp = as.character(sp))
# require(lattice)
# dotplot(ranef(mixed.model.1))  ## default
#
# gro.long.mod.2 <- gro.long.mod %>%
#   left_join(ranef.sp.1, by = "sp") %>%
#   subset(growth.med.sp != 0) %>%
#   mutate(predict.growth.allsp = exp(predict(lm.model.1,
#                                             gro.long.mod %>% subset(growth.med.sp != 0) %>% select(dbh) %>% rename(dbh.mean.allsp = dbh))),
#          predict.growth.sp = exp(predict(mixed.model.1,
#                                          gro.long.mod %>% subset(growth.med.sp != 0) %>% select(dbh, sp) %>% rename(dbh.mean.sp = dbh))))
#            # exp(as.numeric(fixef(mixed.model.half1)[1]) +
#            #       intercept)*dbh^as.numeric(fixef(mixed.model.half1)[2]))
# jpeg("figures/dbh/powerlaw.jpeg")
# plot(x = gro.long.mod.med.allsp$dbh.mean.allsp*0.1, y = gro.long.mod.med.allsp$growth.med.allsp*0.1, log = "xy",
#      main = "All data", ylab = "Growth rate (cm/yr)", xlab = "DBH (cm)")
# points(x = gro.long.mod.2$dbh.mean.allsp*0.1, y = gro.long.mod.2$predict.growth.allsp*0.1,  col = "red")
# graphics.off()
#####---------------

##--------------------------------------------------
## Defining size at the beginning of selected censuses
## (ie. here defnied at census 3, not for each census after)
## And selecting only complete cases stems with growth records in all selected census intervals
##--------------------------------------------------
row.names(gro.wide) <- gro.wide$stemID
gro.wide.cc <- gro.wide %>% subset(complete.cases(gro.wide %>% select(-stemID, -sp))) ## 87922 stem IDs

dbh.long.cc <- dbh_low.out %>% select(-census.8) %>% mutate(stemID = gro.wide$stemID) %>%
  subset(stemID %in% row.names(gro.wide.cc)) %>%
  pivot_longer(cols = starts_with("census."), names_to = "interval",
                           names_prefix = "census.", values_to = "dbh") %>%
  arrange(stemID, interval)

size.class.long.cc <- size.class %>%
  mutate(stemID = gro.wide$stemID) %>%
  subset(stemID %in% row.names(gro.wide.cc)) %>%
  # Defining size at the beginning of census 3
  select(stemID, interval.3) %>%
  pivot_longer(cols = starts_with("interval."), names_to = "interval",
               names_prefix = "interval.", values_to = "size") %>%
  arrange(stemID, interval)

gro.long.cc <- data.frame(gro.wide) %>%
  pivot_longer(cols = starts_with("interval."), names_to = "interval",
               names_prefix = "interval.", values_to = "growth") %>%
  subset(stemID %in% row.names(gro.wide.cc)) %>%
  arrange(stemID, interval) %>%
  ## dbh in the beginning of each census interval
  mutate(dbh.predict.growth = predict(glm.model,
                                      newdata = data.frame(dbh.mean.allsp = dbh.long.cc$dbh)),
                  dbh.residuals = growth - dbh.predict.growth,
        dbh = dbh.long.cc$dbh,
        ## size class in the initial census repeated for all censuses
        size = rep(size.class.long.cc$size, each = intervals)) %>%
  arrange(sp, stemID, interval) %>%
  unite("sp_size", sp, size, remove = FALSE) %>%
  mutate(interval = as.numeric(interval))

gro.long.cc.sub <- gro.long.cc %>% subset(interval == 3 & dbh < 50)
g10 <- ggplot(gro.long.cc %>% subset(stemID %in% gro.long.cc.sub$stemID),
             aes(x = interval)) +
  scale_y_log10() +
  ylab("DBH (cm)")  + xlab("Interval")
#geom_line(data = gro.long.mod.med.allsp %>% subset(interval == 7), aes(y = predict.growth.allsp/10), color = "red")
g10 + geom_line(aes(y = dbh/10, group = stemID, color = as.factor(stemID)), show.legend = FALSE)
g10 + geom_line(aes(y = dbh.residuals/10, group = stemID, color = as.factor(stemID)), show.legend = FALSE)

# this takes a long time
gro.long.cc.norm.stem <- gro.long.cc %>% group_by(stemID) %>%
  mutate(dbh.resid.range = range01(dbh.residuals),
         growth.range = range01(growth),
         dbh.resid.scale = scale(dbh.residuals),
         growth.scale = scale(growth),
         dbh.resid.center = scale(dbh.residuals, scale = FALSE),
         growth.center = scale(growth, scale = FALSE)) %>%
  ungroup(stemID)
## subsetting based on sample size
sample.size <- 3 # at least those number of stems to calculate the mean

gro.long.cc.norm.med <- gro.long.cc.norm.stem %>%
  group_by(sp, size, sp_size, interval) %>%
  summarise(n = n()/intervals,
            med.dbh.resid.range = median(dbh.resid.range, na.rm = TRUE),
            med.growth.range = median(growth.range, na.rm = TRUE),
            med.dbh.resid.scale = median(dbh.resid.scale, na.rm = TRUE),
            med.growth.scale = median(growth.scale, na.rm = TRUE),
            med.dbh.resid.center = median(dbh.resid.center, na.rm = TRUE),
            med.growth.center = median(growth.center, na.rm = TRUE)) %>%
  subset(n >= sample.size) %>%
  ungroup(sp_size, interval)

gro.long.cc.med <- gro.long.cc %>%
  group_by(sp, size, sp_size, interval) %>%
  summarise(med.dbh.resid = median(dbh.residuals, na.rm = TRUE),
            med.growth = median(growth, na.rm = TRUE)) %>%
  subset(sp_size %in% gro.long.cc.norm.med$sp_size) %>%
  ungroup(sp_size, interval)

g11 <- ggplot(gro.long.cc.norm %>% subset(sp_size == "alsebl_large"),
              aes(x = interval)) +
  geom_line(aes(y = dbh.resid.range, group = stemID, color = as.factor(stemID)), show.legend = FALSE) +
  geom_line(data = gro.long.cc.norm.med %>% subset(sp_size == "alsebl_large"), aes(y = med.dbh.resid.range), color = "black", lwd = 1) +
  geom_line(data = gro.long.cc.med %>% subset(sp_size == "alsebl_large"), aes(y = med.dbh.resid), color = "red", lwd = 1)
g11
ggsave("dbh.residuals_stem_range01 and scaled_median range01 in black and median raw overlaid_in red__alsebl_large.jpeg", plot = g11, path =
         file.path("figures/dbh"), height = 5, width = 5, units='in')

g12 <- ggplot(gro.long.cc.norm %>% subset(sp_size == "alsebl_large"),
              aes(x = interval)) +
  geom_line(aes(y = dbh.resid.scale, group = stemID, color = as.factor(stemID)), show.legend = FALSE) +
  geom_line(data = gro.long.cc.norm.med %>% subset(sp_size == "alsebl_large"), aes(y = med.dbh.resid.scale), color = "black", lwd = 1) +
  geom_line(data = gro.long.cc.med %>% subset(sp_size == "alsebl_large"), aes(y = med.dbh.resid), color = "red", lwd = 1)
g12
ggsave("dbh.residuals_stem_center and scaled_median center in black and scaled and median raw overlaid_in red_alsebl_large.jpeg", plot = g12, path =
         file.path("figures/dbh"), height = 5, width = 5, units='in')
g13 <- ggplot(gro.long.cc.norm %>% subset(sp_size == "alsebl_large"),
              aes(x = interval)) +
  geom_line(aes(y = dbh.resid.center, group = stemID, color = as.factor(stemID)), show.legend = FALSE) +
  geom_line(data = gro.long.cc.norm.med %>% subset(sp_size == "alsebl_large"), aes(y = med.dbh.resid.center), color = "black", lwd = 1) +
  geom_line(data = gro.long.cc.med %>% subset(sp_size == "alsebl_large"), aes(y = med.dbh.resid), color = "red", lwd = 1)
g13
ggsave("dbh.residuals_stem_centered_median centered in black and median raw overlaid_in red__alsebl_large.jpeg", plot = g13, path =
         file.path("figures/dbh"), height = 15, width = 5, units='in')
### centering and scaling growth for each stem

growth.wide.cc.1 <- gro.long.cc %>% select(sp_size, stemID, interval, growth) %>%
  pivot_wider(names_from = c("interval"), values_from = "growth")
growth.wide.cc <- growth.wide.cc.1 %>% select(-sp_size, -stemID)
growth.wide.cc.norm <- data.frame(t(apply(growth.wide.cc, 1, scale, center = TRUE, scale = TRUE)))
colnames(growth.wide.cc.norm) <- colnames(growth.wide.cc)
sp_size.stem.cc.growth.list <- split(growth.wide.cc.norm, growth.wide.cc.1$sp_size)
sp_size.stem.cc.growth <- lapply(sp_size.stem.cc.growth.list, function(x){
  x %>% pivot_longer(everything(), names_to = c("interval"), values_to = "growth") %>%
    mutate(interval = as.numeric(interval))})

### centering and scaling the residuals for each stem
dbh.res.wide.cc.1 <- gro.long.cc %>% select(sp_size, stemID, interval, dbh.residuals) %>%
  pivot_wider(names_from = c("interval"), values_from = "dbh.residuals")
dbh.res.wide.cc <- dbh.res.wide.cc.1 %>% select(-sp_size, -stemID)
dbh.res.wide.cc.norm <- data.frame(t(apply(dbh.res.wide.cc, 1, scale, center = TRUE, scale = TRUE)))
colnames(dbh.res.wide.cc.norm) <- colnames(dbh.res.wide.cc)
sp_size.stem.cc.dbh.res.list <- split(dbh.res.wide.cc.norm, dbh.res.wide.cc.1$sp_size)

sp_size.stem.cc.dbh.res <- lapply(sp_size.stem.cc.dbh.res.list, function(x){
  x %>% pivot_longer(everything(), names_to = c("interval"), values_to = "dbh.residuals") %>%
    mutate(interval = as.numeric(interval))})

sp_size.stem.cc.n.dbh.res.list <- lapply(sp_size.stem.cc.dbh.res, function(x){
  data.frame(n = nrow(x)/intervals)})
sp_size.n.cc.dbh.res  <- rbindlist(sp_size.stem.cc.n.dbh.res.list, idcol = "sp_size")
select.sp_size <- sp_size.n.cc.dbh.res$sp_size[sp_size.n.cc.dbh.res$n >= sample.size]


sp_size.stats.cc.dbh.res <- lapply(sp_size.stem.cc.dbh.res[select.sp_size], function(x){
  x %>% group_by(interval) %>%
    summarise(median = median(dbh.residuals, na.rm = TRUE),
              mean = median(dbh.residuals, na.rm = TRUE),
              sd = sd(dbh.residuals, na.rm = TRUE),
              trees = n(),
              se = sd/sqrt(trees),
              .groups = "drop_last") %>%
    arrange(interval)
})
sp_size.stats.cc.growth <- lapply(sp_size.stem.cc.growth[select.sp_size], function(x){
  x %>% group_by(interval) %>%
    summarise(median = median(growth, na.rm = TRUE),
              mean = mean(growth, na.rm = TRUE),
              sd = sd(growth, na.rm = TRUE),
              trees = n(),
              se = sd/sqrt(trees),
              .groups = "drop_last") %>%
    arrange(interval)
})
## 493 sp_size has at least stems = sample.size

sp_size.stem.cc.names <- names(sp_size.stem.cc.growth) # 660
sp_size.stats.cc.names <- names(sp_size.stats.cc.growth) # 520
growth.selection <- "size_class_predefined_cc_scaled"

save(gro.long.cc, file = paste0("results/gro.long.cc_", intervals, "_", growth.selection, ".Rdata"))
save(gro.long.cc.norm.stem , file = paste0("results/gro.long.cc.norm.stem _", intervals, "_", growth.selection, ".Rdata"))
save(gro.long.cc.norm.med, file = paste0("results/gro.long.cc.norm.med_", intervals, "_", growth.selection, ".Rdata"))
save(gro.long.cc.med, file = paste0("results/gro.long.cc.med_", intervals, "_", growth.selection, ".Rdata"))

save(sp_size.stem.cc.dbh.res, file = paste0("results/sp_size.individual_growth_dbh.residuals_on_", intervals, "_", growth.selection, ".Rdata"))
save(sp_size.stats.cc.dbh.res, file = paste0("results/sp_size.stats_growth_dbh.residuals_on_", intervals, "_", growth.selection, ".Rdata"))
save(sp_size.stem.cc.growth, file = paste0("results/sp_size.individual_growth_dbh.residuals_off_", intervals, "_", growth.selection, ".Rdata"))
save(sp_size.stats.cc.growth, file = paste0("results/sp_size.stats_growth_dbh.residuals_off_", intervals, "_", growth.selection, ".Rdata"))
save(sp_size.stem.cc.names, file = paste0("results/sp_size.individual.names_", intervals, "_", growth.selection, ".Rdata"))
save(sp_size.stats.cc.names, file = paste0("results/sp_size.stats.names_", intervals, "_", growth.selection, ".Rdata"))


##--------------------------------------------------

##--------------------------------------------------
## Plotting stem growth time series
##--------------------------------------------------

# graphics info
theme_set(theme_bw())
theme_update(text = element_text(size = 14),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)

load(file.path("data-raw/swp.gfac.rda"))
load(file.path("data-raw/psi.rda"))
census.meds <- readr::read_rds("results/census.mediandates.rds")
census.beg <- census.meds[3: length(census.meds)]
cut.breaks <- census.beg
cut.labels.2 <- paste0(seq(1990, 2010, by = 5), "-", seq(1995, 2015, by = 5))

swp.gfac <- swp.gfac %>%
  mutate(interval.yrs = cut(date, include.lowest = TRUE, breaks = cut.breaks,
                            labels = cut.labels.2, right = TRUE))

swp.gfac.stat.2 <- swp.gfac %>%
  subset(!is.na(interval.yrs)) %>%
  group_by(interval.yrs) %>%
  summarise(mean = mean(gfac, na.rm = TRUE)) %>%
  droplevels()
plot.swp.gfac.stat.3 <- ggplot(swp.gfac.stat.2, aes(x = interval.yrs, y = mean)) +
  geom_point() + ylab("Mean Growth factor (0-1)") + xlab("Census Interval") +
  scale_fill_viridis_c("Growth\nFactor\n[unitless, 0-1]", trans = "reverse", option = "plasma") +
  ggtitle("Average Growth Factor by interval and depth across best-fit ensembles")
ggsave("swp.gfac_mean_across_params.top.few_full_interval_no-depth.jpeg", plot = plot.swp.gfac.stat.3, path =
         file.path("figures"), height = 5, width = 5, units='in')

psi <- psi %>%
  mutate(interval.yrs = cut(date, include.lowest = TRUE, breaks = cut.breaks,
                            labels = cut.labels.2, right = TRUE))

psi.stat.2 <- psi %>%
  subset(!is.na(interval.yrs)) %>%
  group_by(interval.yrs) %>%
  summarise(mean = mean(psi, na.rm = TRUE)) %>%
  droplevels()
plot.psi.stat.3 <- ggplot(psi.stat.2, aes(x = interval.yrs, y = mean)) +
  geom_point() + ylab("PSI (MPa)") + xlab("Census Interval") +
  scale_fill_viridis_c("Soil\nWater\nPotential\n[MPa]", trans = "reverse", option = "plasma") +
  ggtitle("Mean Soil Water Potential by interval and depth across best-fit ensembles")
ggsave("psi_mean_across_params.top.few_full_interval_no-depth.jpeg", plot = plot.psi.stat.3, path =
         file.path("figures"), height = 5, width = 5, units='in')

####------------------------------------------------------
### Also detrending for dbh and solar radiation in the same loop
####------------------------------------------------------
### Boris's data - 2018 substituted by rutuja
####------------------------------------------------------

## data starts only from 1985-01-01
clim.data <- read.csv("data-raw/BCI_1985_2018c_mod_2018substituted.csv")
str(clim.data)
census.meds <- readr::read_rds("results/census.mediandates.rds")
cut.breaks <- census.meds
cut.labels <- 1:(length(census.meds) -1)

clim.data <- clim.data %>%
  # converting to system time zone
  mutate(datetime = as.POSIXct(clim.data$DateTime, format = "%m/%d/%y %H:%M", tz = "America/Panama"),
         date = as.Date(datetime),
         ## adding interval label to each date
         interval = cut(date, include.lowest = TRUE, breaks = cut.breaks,
                        labels = cut.labels, right = FALSE)) %>% subset(date %in% cut.breaks[1]: cut.breaks[length(cut.breaks)])
sola <- clim.data %>%
  group_by(interval) %>%
  summarise(solar = mean(SR_W_m2, na.rm = TRUE)) %>%
  mutate(interval = as.numeric(interval)) %>% data.table()
####------------------------------------------------------

####------------Loop starts-----------------------------

growth.types <- c("individual", "mean")
growth.selections <- c("size_class_predefined_cc", "size_class_varying_non_cc")
for (k in 1: length(growth.selection)){
  for (i in 1: length(growth.types)) {
    growth.type <- growth.types[i]
    growth.selection <- growth.selections[k]

    g.list.dfname <- load(file =  paste0("results/sp_size.", growth.type, "_back_transformed_lambda_growth_", intervals, "_", growth.selection, ".Rdata"))
    g.names.dfname <- load(file =  paste0("results/sp_size.", growth.type, ".names_", intervals, "_", growth.selection, ".Rdata"))
    g.list <- get(g.list.dfname)
    g.names <- get(g.names.dfname)
    level_key <- c(`3` = "1990-1995", `4` = "1995-2000", `5` = "2000-2005", `6` = "2005-2010", `7` = "2010-2015")

    for (j in 1 : length(g.list)) {
      g.list[[j]]$sp_size <- g.names[j]
    }
    ###------------------------------------
    ## Plotting stem growth time series----
    ###------------------------------------
    gdf <- do.call(rbind, g.list)
    gdf <- gdf %>%
      separate(sp_size, into = c("sp", "size"), sep = "_", remove = FALSE) %>%
      mutate(census = as.factor(recode(interval, !!!level_key)))

    sp.n <- read.csv(file.path(paste("results/sp.n_med_growth_sp_by_size_6.csv")), row.names = 1)
    select.sp.df <- sp.n %>% subset(size == "large") %>% arrange(desc(n))
    gdf <- gdf %>%
      left_join(swp.gfac.stat.2 %>%
                  rename(gfac = mean, census = interval.yrs), by = "census")
    gdf.plot.1 <- ggplot(gdf %>% subset(sp %in% as.character(select.sp.df$sp[1:20])),
                         aes(x = gfac, y = growth, color = size)) +
      facet_wrap(size ~ sp, scale = "free_y") +
      geom_point(alpha = 0.3) +
      geom_smooth(method = "lm", lty = "dashed") +
      xlab("Mean Growth factor (0-1)") + ylab("Growth rate (cm/yr)")
    ggsave(file.path(paste0("figures/", growth.type, "_growth_by_gfac_by_size_", growth.selection, ".jpeg")), plot = gdf.plot.1,
           height = 15, width = 30, units ='in')
    ###---------------------------
    ### Detrending for dbh--------
    ###---------------------------
  }
}


### detrending for solar radiation
###-------------------
### Boris's data - 2018 substituted by rutuja
###-------------------
## data starts only from 1985-01-01
clim.data <- read.csv("data-raw/BCI_1985_2018c_mod_2018substituted.csv")
str(clim.data)
census.meds <- readr::read_rds("results/census.mediandates.rds")
cut.breaks <- census.meds
cut.labels <- 1:(length(census.meds) -1)

clim.data <- clim.data %>%
  # converting to system time zone
  mutate(datetime = as.POSIXct(clim.data$DateTime, format = "%m/%d/%y %H:%M", tz = "America/Panama"),
         date = as.Date(datetime),
  ## adding interval label to each date
         interval = cut(date, include.lowest = TRUE, breaks = cut.breaks,
                 labels = cut.labels, right = FALSE)) %>% subset(date %in% cut.breaks[1]: cut.breaks[length(cut.breaks)])
sola <- clim.data %>%
  group_by(interval) %>%
  summarise(solar = mean(SR_W_m2, na.rm = TRUE)) %>%
  mutate(interval = as.numeric(interval)) %>% data.table()

growth.types <- c("individual", "mean")
intervals <- 5

for (i in 1 : 2) {
  growth.type <- growth.types[i]
  load(file =  paste0("results/sp_size.", growth.type, "_back_transformed_lambda_growth_", intervals, ".Rdata"))
  growth.sola <- lapply(sp_size.mean, function(x) {
    x <- data.table(x)[sola, on = "interval", `:=`(solar = i.solar)]
    lm.model <- lm(growth ~ solar, data = x)
    x[, residuals := lm.model$residuals]
  })
  save(growth.sola, file = paste0("results/sp_size.", growth.type, "_back_transformed_lambda_growth_&_residuals_", intervals, ".Rdata"))

  residuals <- lapply(growth.sola, function(x) {
    x[, c("growth", "solar") := NULL]
    })
  save(residuals, file = paste0("results/sp_size.", growth.type, "_back_transformed_lambda_growth_residuals_", intervals, ".Rdata"))
}


growth.sola.plot <- lapply(sp_size.mean[10:20], function(xx) {
  xx <- data.table(xx)[sola, on = "interval", `:=`(solar = i.solar)]
  lm.model <- lm(growth ~ solar, data = xx)
  xx[, residuals := lm.model$residuals]
  ggplot(xx, aes(y = growth, x = solar)) +
    geom_point() +
    geom_smooth(method = "lm", title)
})

for (i in 1: length(growth.types)) {
  growth.type <- growth.types[i]
  growth.selection <- growth.selections[k]

  g.list.dfname <- load(file =  paste0("results/sp_size.", growth.type, "_back_transformed_lambda_growth_", intervals, "_", growth.selection, ".Rdata"))
  g.names.dfname <- load(file =  paste0("results/sp_size.", growth.type, ".names_", intervals, "_", growth.selection, ".Rdata"))
  g.list <- get(g.list.dfname)
  g.names <- get(g.names.dfname)

  level_key <- c(`3` = "1990-1995", `4` = "1995-2000", `5` = "2000-2005", `6` = "2005-2010", `7` = "2010-2015")

  for (j in 1 : length(g.list)) {
    g.list[[j]]$sp_size <- g.names[j]
  }

  gdf <- do.call(rbind, g.list)
  gdf <- gdf %>%
    separate(sp_size, into = c("sp", "size"), sep = "_", remove = FALSE) %>%
    mutate(census = as.factor(recode(interval, !!!level_key)))

  sp.n <- read.csv(file.path(paste("results/sp.n_med_growth_sp_by_size_6.csv")), row.names = 1)
  select.sp.df <- sp.n %>% subset(size == "large") %>% arrange(desc(n))

  gdf.plot.1 <- ggplot(gdf %>% subset(sp %in% as.character(select.sp.df$sp[1:20])),
                       aes(x = solar, y = growth, color = size)) +
    facet_wrap(size ~ sp, scale = "free_y") +
    geom_smooth(method = "lm", lty = "dashed") +
    geom_point() +
    xlab("Mean Solar Radiation (W/m2)") + ylab("Growth rate (cm/yr)")
  ggsave(file.path(paste0("figures/", growth.types[i], "_growth_by_solar_radiation_by_size.jpeg")), plot = gdf.plot.1,
         height = 15, width = 30, units ='in')
}
