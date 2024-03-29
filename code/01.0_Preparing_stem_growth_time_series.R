#________________________________
# Title: Preparing growth data for inverse modeling
# Author : Rutuja Chitra-Tarak
# Original date: December 11, 2019
#________________________________

# for 50 ha obs species groups
rm(list = ls())
gc()

#*******************************************
####   Load Libraries, Prep for graphics, folders  ####
#*******************************************
#### Written with R version 4 ###
#*******************************************
if (!require("groundhog")) install.packages("groundhog"); library(groundhog)
groundhog.folder <- paste0("groundhog.library")
if(!dir.exists(file.path(groundhog.folder))) {dir.create(file.path(groundhog.folder))}
set.groundhog.folder(groundhog.folder)
groundhog.day = "2021-01-01" # "2020-05-01"
pkgs=c('tidyverse', 'purrrlyr', 'scales', 'grid', 'gridExtra', 'lme4')
groundhog.library(pkgs, groundhog.day)


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
#________________________________
## Collating tree.full data for all censuses
#________________________________
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
stem.full$census <- rep(census, each = nrow(bci.stem3))
#________________________________
## adding crown exposure---------
#________________________________
# crown1 <- read.csv("data-raw/traits/CrownExposure_MatteoDetto/bci1.csv") %>% mutate(census = 1)
# crown2 <- read.csv("data-raw/traits/CrownExposure_MatteoDetto/bci2.csv") %>% mutate(census = 2)
crown3 <- read.csv("data-raw/traits/CrownExposure_MatteoDetto/bci3.csv") %>% mutate(census = 3)
crown4 <- read.csv("data-raw/traits/CrownExposure_MatteoDetto/bci4.csv") %>% mutate(census = 4)
crown5 <- read.csv("data-raw/traits/CrownExposure_MatteoDetto/bci5.csv") %>% mutate(census = 5)
crown6 <- read.csv("data-raw/traits/CrownExposure_MatteoDetto/bci6.csv") %>% mutate(census = 6)
crown7 <- read.csv("data-raw/traits/CrownExposure_MatteoDetto/bci7.csv") %>% mutate(census = 7)
crown8 <- read.csv("data-raw/traits/CrownExposure_MatteoDetto/bci8.csv") %>% mutate(census = 8)

crown.full <- bind_rows(crown3, crown4, crown5, crown6, crown7, crown8) %>%
  rename(crown = rank)
stem.full <- stem.full %>% left_join(crown.full %>% select(-dbh, -sp), by = c("treeID", "census"))
#________________________________
## remove duplicated tags if any
nrow(stem.full) # 4740360
stem.full <- stem.full %>% mutate(ExactDate = as.Date(ExactDate))
stem.full <- stem.full %>% group_by(census) %>% distinct(stemID, .keep_all= TRUE)
save(stem.full, file = "results/stem.full_crown.rda")

nrow(stem.full)
# 4740360
# so no duplicate stemIDs;
n_stem.full <- stem.full %>% group_by(census) %>% nest()
map_dbl(n_stem.full$data, nrow)
# 790060 790060 790060 790060 790060 790060 790060

#________________________________
## to remove stems for which hom have changed--------
#________________________________
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
#_____End of removing stems with hom change#____________


#________________________________
# Getting a matrix of date by census for these stems
#________________________________
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
#________________________________End________________________________


#________________________________
## Now to obtain growth rates-----
#________________________________
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
#________________End_____________________________________________________
jpeg("figures/dbh/stem_growth_without_outliers.jpeg")
plot(x = as.vector(as.matrix(dbh_low.out[, - ncol(dbh_low.out)])), y = as.vector(as.matrix(g_up.out)),
     ylab = expression("Growth rate (mmyr"^-1*")"), xlab = "DBH (mm)")
graphics.off()

#
# ##________________________________
# ## Saving individual stem data for HADAD with following selection criterion:
# ## stems with at least 3 records
# #________________________________
# intervals <- nint; #census.cols = (ncol(gt_back) - intervals + 1):ncol(gt_back)
# rows.3obs <- unlist(apply(gt_back, 1, function(x) {sum(!is.na(x)) >= 3}))
# #________________________________
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
# #____________________End___________________________________________

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
#________________________________
## a stem's size class can change across censuses
## and thus which size class group it belongs to, whose mean is being taken
#________________________________
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

#________________________________
### getting dbh residuals--------
#________________________________
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
# Call:
#   glm(formula = growth.mean.allsp ~ splines::bs(dbh.mean.allsp,
#                                                 5), data = gro.long.mod.med.allsp)
#
# Deviance Residuals:
#   Min       1Q   Median       3Q      Max
# -1.3968  -0.1374  -0.0169   0.1150   3.1867
#
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                      0.74826    0.08408   8.899  < 2e-16 ***
#   splines::bs(dbh.mean.allsp, 5)1  0.18948    0.12911   1.468  0.14319
#   splines::bs(dbh.mean.allsp, 5)2 -0.11721    0.10447  -1.122  0.26268
#   splines::bs(dbh.mean.allsp, 5)3  2.34282    0.24441   9.586  < 2e-16 ***
#   splines::bs(dbh.mean.allsp, 5)4 -0.83232    0.29359  -2.835  0.00487 **
#   splines::bs(dbh.mean.allsp, 5)5  1.53714    0.25768   5.965 6.31e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# (Dispersion parameter for gaussian family taken to be 0.165665)
#
# Null deviance: 90.606  on 334  degrees of freedom
# Residual deviance: 54.504  on 329  degrees of freedom
# AIC: 356.38
#
# Number of Fisher Scoring iterations: 2
summary(lm.model)
### using log(growth.mean.allsp)
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
# without color
g5 <- ggplot(gro.long.mod.med.allsp,
             aes(x = dbh.mean.allsp/10, y = growth.mean.allsp/10)) +
  geom_point(alpha = 0.7) +
  scale_y_log10() + scale_x_log10(breaks = c(1, 2, 5, 10, 20, 30, 50, 100)) +
  # ggtitle("All data") +
  ylab(expression("Growth rate (cmyr"^-1*")"))  + xlab("DBH (cm)") +
  stat_poly_eq(aes(label = paste(stat(eq.label))),
               npcx = 0.05, npcy = 0.97, rr.digits = 2,
               formula = formula.glm, parse = TRUE, size = 3) +
  stat_poly_eq(aes(label = paste(stat(adj.rr.label))),
               npcx = 0.05, npcy = 0.90, rr.digits = 2,
               formula = formula.glm, parse = TRUE, size = 3) +
  geom_smooth(method = glm, formula = y ~ splines::bs(x, 5), size = 0.7)
ggsave("growthrate_dbh_predicted_growth_with_spline_degree_5_dbh_cutoff_100_model_no_color.jpeg",
       plot = g5, path = file.path("figures/dbh/"), height = 4, width = 4.5, units='in')
# > 30 cm
g6 <- g4 + coord_cartesian(xlim = c(30, max(gro.long.mod.med.allsp$dbh.mean.allsp/10, na.rm = TRUE)))
ggsave("growthrate_dbh_predicted_growth_with_spline_degree_5_dbh_cutoff_100_model_above30cm.jpeg",
       plot = g6, path = file.path("figures/dbh/"), height = 4.5, width = 4.5, units='in')

## not enough species-wise data to fit a trend...often giving rise to negative trends
## So using community wide fit
gro.long.mod <- gro.long.mod %>%
  mutate(dbh.predict.growth = predict(glm.model, newdata = data.frame(dbh.mean.allsp = gro.long.mod$dbh)),
         dbh.residuals = growth - dbh.predict.growth)
## throws a warning
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

## not enough species-wise data to fit a trend...often giving rise to negative trends
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
#________________________________

#_____________________________________________________
## Defining size at the beginning of selected censuses
## (ie. here defnied at census 3, not for each census after)
## And selecting only complete cases stems with growth records in all selected census intervals----
#_____________________________________________________
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
         # throws a warning
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
            med.growth.center = median(growth.center, na.rm = TRUE), .groups = "drop_last") %>%
  subset(n >= sample.size) %>%
  ungroup(sp_size, interval)

gro.long.cc.med <- gro.long.cc %>%
  group_by(sp, size, sp_size, interval) %>%
  summarise(med.dbh.resid = median(dbh.residuals, na.rm = TRUE),
            med.growth = median(growth, na.rm = TRUE), .groups = "drop_last") %>%
  subset(sp_size %in% gro.long.cc.norm.med$sp_size) %>%
  ungroup(sp_size, interval)

g11 <- ggplot(gro.long.cc.norm.stem %>% subset(sp_size == "alsebl_large"),
              aes(x = interval)) +
  geom_line(aes(y = dbh.resid.range, group = stemID, color = as.factor(stemID)), show.legend = FALSE) +
  geom_line(data = gro.long.cc.norm.med %>% subset(sp_size == "alsebl_large"), aes(y = med.dbh.resid.range), color = "black", lwd = 1) +
  geom_line(data = gro.long.cc.med %>% subset(sp_size == "alsebl_large"), aes(y = med.dbh.resid), color = "red", lwd = 1)
g11
ggsave("dbh.residuals_stem_range01 and scaled_median range01 in black and median raw overlaid_in red__alsebl_large.jpeg", plot = g11, path =
         file.path("figures/dbh"), height = 5, width = 5, units='in')

g12 <- ggplot(gro.long.cc.norm.stem %>% subset(sp_size == "alsebl_large"),
              aes(x = interval)) +
  geom_line(aes(y = dbh.resid.scale, group = stemID, color = as.factor(stemID)), show.legend = FALSE) +
  geom_line(data = gro.long.cc.norm.med %>% subset(sp_size == "alsebl_large"), aes(y = med.dbh.resid.scale), color = "black", lwd = 1) +
  geom_line(data = gro.long.cc.med %>% subset(sp_size == "alsebl_large"), aes(y = med.dbh.resid), color = "red", lwd = 1)
g12
ggsave("dbh.residuals_stem_center and scaled_median center in black and scaled and median raw overlaid_in red_alsebl_large.jpeg", plot = g12, path =
         file.path("figures/dbh"), height = 5, width = 5, units='in')
g13 <- ggplot(gro.long.cc.norm.stem %>% subset(sp_size == "alsebl_large"),
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
dbh.res.wide.cc.1 <- gro.long.cc %>%
  # retaining stems that are limited to canopy exposure > 80 removes most of the trees, so not employing
  # subset(stemID %in% crown.wide.select$stemID) %>%
  select(sp_size, stemID, interval, dbh.residuals) %>%
  pivot_wider(names_from = c("interval"), values_from = "dbh.residuals")
dbh.res.wide.cc <- dbh.res.wide.cc.1 %>% select(-sp_size, -stemID)
dbh.res.wide.cc.norm <- data.frame(t(apply(dbh.res.wide.cc, 1, scale, center = TRUE, scale = TRUE)))
colnames(dbh.res.wide.cc.norm) <- colnames(dbh.res.wide.cc)
table(dbh.res.wide.cc.1$sp_size)

sp_size.stem.cc.dbh.res.list <- split(dbh.res.wide.cc.norm, dbh.res.wide.cc.1$sp_size)

sp_size.stem.cc.dbh.res <- lapply(sp_size.stem.cc.dbh.res.list, function(x){
  x %>% pivot_longer(everything(), names_to = c("interval"), values_to = "dbh.residuals") %>%
    mutate(interval = as.numeric(interval))})

sp_size.stem.cc.n.dbh.res.list <- lapply(sp_size.stem.cc.dbh.res, function(x){
  data.frame(n = nrow(x)/intervals)})
sp_size.n.cc.dbh.res  <- rbindlist(sp_size.stem.cc.n.dbh.res.list, idcol = "sp_size")
select.sp_size <- sp_size.n.cc.dbh.res$sp_size[sp_size.n.cc.dbh.res$n >= sample.size]

nrow(select.sp_size)

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

