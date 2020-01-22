# for 50 ha obs species groups
rm(list = ls())
source(file = "R/utilities.R")
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, purrrlyr, scales, grid, gridExtra)
# graphics info
tex <- element_text(size = 16, face = "plain") # , family = "gara"
tex <- element_text(size = 5, face = "plain")
my.theme <-  theme(axis.text = tex, axis.title = tex,
                   title = tex, legend.title = tex, legend.text = tex, strip.text.y = tex, strip.text.x = tex)
my.bg <- theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "black", size = 0.25),
                            panel.grid.minor = element_blank())
my.adjust <- theme(axis.title.y = element_text(vjust = 1), axis.title.x = element_text(vjust = -0.6), title = element_text(vjust = 2))

# time series info
years <- c(1982, seq(1985, 2015, by = 5))
ncensus  <- length(years)
nint <- length(years) - 1
# sp4 <- c("ANOL", "TECG", "LAGL", "TERT")
##------------------------------------------------
## Collating tree.full data for all censuses
##------------------------------------------------
## All data frames have same number of rows. Each row represent a tree,
## and a given row number across censuses represent the same tree
load("data-raw/CTFScensuses/BCI.tree1.Rdata")
load("data-raw/CTFScensuses/BCI.tree2.Rdata")
load("data-raw/CTFScensuses/BCI.tree3.Rdata")
load("data-raw/CTFScensuses/BCI.tree4.Rdata")
load("data-raw/CTFScensuses/BCI.tree5.Rdata")
load("data-raw/CTFScensuses/BCI.tree6.Rdata")
load("data-raw/CTFScensuses/BCI.tree7.Rdata")
load("data-raw/CTFScensuses/BCI.tree8.Rdata")

length(unique(bci.tree1$sp))
# all censuses in one list
tree.full <- bind_rows(bci.tree2, bci.tree3, bci.tree4, bci.tree5, bci.tree6, bci.tree7, bci.tree8)
# > nrow(bci.tree.full1)
# [1] 423617
# > nrow(bci.tree.full1)*8
# [1] 3388936
###------------selection----------
census <- 2:8 # keeping 1985-1990 census for testing
###-------------------------------
tree.full$census <- rep(census, each = nrow(bci.tree2))
## remove duplicated tags if any
nrow(tree.full)
tree.full <- tree.full %>% mutate(ExactDate = as.Date(ExactDate))
tree.full <- tree.full %>% group_by(census) %>% distinct(tag, .keep_all= TRUE)
nrow(tree.full)
# so no duplicate tags; may be the largest stem
n_tree.full <- tree.full %>% group_by(census) %>% nest()
map_dbl(n_tree.full$data, nrow)
# 423617

## to remove tree.fulls for which hom have changed--------
##------------------------------------------------
# pool together for each tree.full hom in censuses
hom.mat <- map_dfc(n_tree.full$data, "hom")
hom.mat <- data.frame(invisible(lapply(hom.mat, function (x) {as.numeric(x)})))
# check by row (tree.full) if hom value remains the same
hom.same <- apply(hom.mat, 1, function (x) {diff(range(x, na.rm = T)) == 0})
# xx <- dbh.mat %>%
#   rowwise() %>% mutate(diff = function (x) {diff(range(x, na.rm = T)) == 0})
## if hom has changed, resut is FALSE
length(which(hom.same == FALSE))
# 28959
## those with all NA rows are false positively identified as having hom changed, so making those TRUE
hom.same[rowSums(is.na(hom.mat)) == length(censuses)] <- TRUE
length(which(hom.same == FALSE))
## 5515
# length(which(hom.same == FALSE))
# [1] 6027 ## this is the number when census 1 is included; hom.same <- apply(hom.mat, 1, function (x) {diff(range(x, na.rm = T)) == 0})
## but most pom change has happened from census one to census 2, and we are not interested in including census 1

## tree.fulls with unchanged hom
keep.tags <- n_tree.full$data[[1]]$tag[hom.same]
length(keep.tags) #394658
# only selecting those tree.fulls out of
map_dbl(n_tree.full$data, nrow)
# [1] 423617 423617 423617 423617 423617 423617 423617
tree <- tree.full %>% group_by(census) %>% filter(tag %in% keep.tags)
n_tree <- tree %>% group_by(census) %>% nest()
map_dbl(n_tree$data, nrow)
# [1] 394658 394658 394658 394658 394658 394658 394658
# 423617 - 394658
# [1] 28959
##--------------------End----------------------------

##---
# mod_fun <- function(df) lm(Sepal.Length ~ ., data = df)
# m_iris <- n_iris %>%
#   mutate(model = map(data, mod_fun))
##-------

##------------------------------------------------
# Getting a matrix of date by census for these trees
##------------------------------------------------
time.mat <- map_dfc(n_tree$data, "ExactDate")
n_date.mat <- map_dfc(n_tree$data, "date")
colnames(time.mat) <- paste("census", census, sep = ".")
# getting a matrix of dbh by census for these trees
dbh.mat <- map_dfc(n_tree$data, "dbh")
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
# 20048
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
# 1] 0.8798
# [1] 0.4037611
## not the same

## Removing extreme positive errors in growth-----
# removing positive errors > 75 mm/yr since that's the growth of the fastest growing trees (true growth, from a double blinded test)
## note extreme records and assign them NAs
g_up.out <- apply(g, 2, function (x) {ifelse(x > 75, NA, x)})
mean(g_up.out, na.rm = T); median(g_up.out, na.rm = T)
# [1] 0.8767508
# [1] 0.4037611
## This matrix needs to be used in further growth calculation
# saving growth without modulus transformation, after all outliers were removed
save(g_up.out, file = "results/tree_growth_without_outliers.Rdata")
##------------------End---------------------------

##------------------------------------------------
## Modulus transformation (growth_lambda = growth^lambda) to rein in negative and positive outliers-------
##------------------------------------------------
## ~From Condit et al 2017~
# "We found that in the range λ ∈ (0.3, 0.6), transformed growth rates have
# low skewness, and median and mean are close. For any
# sample of growth increments, a λ can be located that minimizes
# skewness, but we sought one value that worked reasonably for all species and dbh categories.
# The main purpose is to reign in the big outliers that can cause peculiar model results,
# and λ = 0.4 was satisfactory for both saplings and large trees."
## ~.~
##------------------------------------------------

lambda <- 0.4
g_sign <- sign(g_up.out)
# g_sign
# [1] -1 1 1 -1
g_t <- g_sign*(g_sign*g_up.out)^lambda
mean(g_t, na.rm = T); median(g_t, na.rm = T)
# [1] 0.6727
# [1] 0.6957445
## mean, median similar-------
gt_back <- g_t^(1/lambda)
# g_t
# [1] -1.000000 1.319508 1.551846 -1.741101
g_up.out[1:30,]; gt_back[1:30,]
mean(gt_back, na.rm = T); median(gt_back, na.rm = T)
# [1] 0.9637
# [1] 0.4291
mean(g_up.out, na.rm = T); median(g_up.out, na.rm = T)
# [1] 0.8791
# [1] 0.4037
## medians are supposed to be exact for raw and back transformed growth, but they arent? They are close, yes
save(gt_back, file = "results/tree_back_transformed_lambda_growth.Rdata")
load(file = "results/tree_back_transformed_lambda_growth.Rdata")
##------------------End------------------------------


##--------------------------------------------------
## Saving individual tree data for HADAD with following selection criterion:
## excluding census 1985-1990; that is selecting last 5 intervals
## trees with at least 3 records
##--------------------------------------------------
intervals <- 5; census.cols = (ncol(gt_back) - intervals + 1):ncol(gt_back)
rows.3obs <- rowSums(!is.na(gt_back[, census.cols])) >= 3
###---------------
sum(rows.3obs); #1,58,107 # with 1985-90 interval, 1,78,973
sum(rows.3obs)/nrow(gt_back) # 40% of trees have at least 3 records
gt_back_for_hadad <- gt_back[rows.3obs,]
##---
## that is this data has 1985-1990 census growth (6 census intervals), but the tree selection with at least 3 records was performed on last 5 census intervals (excluding census 1985-1990)
##---
save(gt_back_for_hadad, file =  paste0("results/growth_individual_tree_intervals_", intervals ,".Rdata"))
# Now saving tree data for these selected trees, dbh for all censuses are saved
tree.data <- n_tree$data[[census.cols[1]]][rows.3obs, 1:8]
save(tree.data, file =  paste0("results/growth_individual_tree.data_intervals_", intervals ,".Rdata"))
load(file = paste0("results/growth_individual_tree_intervals_", intervals ,".Rdata"))
##------------------End------------------------------


##--------------------------------------------------
## Selecting size class based on each census, so a tree's size class identity can change across censuses -------
##--------------------------------------------------

# large(300, 500], medium(100, 150], small (60, 80], very small (10, 50]
## Scaled growth
## this has all the trees data with census intervals included in (variable "census" - 1)
gro <- data.frame(t(scale(t(gt_back), center = TRUE))) #, scale = TRUE # substract mean across time series and divide by sd
save(gro, file = paste0("results/tree_back_transformed_growth_clean_intervals_", length(census) - 1 ,".Rdata"))
load(file = paste0("results/tree_back_transformed_growth_clean_intervals_", length(census) - 1 ,".Rdata"))

nrow(gro)
# gro <- data.frame(gt_back)
cut.breaks <- c(10, 50, 100, 300, max(dbh_low.out, na.rm = T)) # about 45 trees per interval are > 1000 mm in dbh
cut.labels <- c("tiny", "small", "medium", "large")
# old cut.breaks <- c(10, 50, 60, 80, 100, 150, 300, 700)
# cut.labels <- c("tiny", "out", "small", "out", "medium", "out", "large")
size.class <- data.frame(apply(dbh_low.out[, -ncol(dbh_low.out)], 2, cut,
                               include.lowest = TRUE, breaks = cut.breaks,
                               labels = cut.labels, right = TRUE))
colnames(size.class) <- paste("interval", census[-length(census)], sep = ".")
## gro and data for one census in n_tree (one list element) have same tree identities by row, so species id is added from n_tree
gro.wide <- gro %>% mutate(tree = n_tree$data[[1]]$treeID) %>% mutate(sp = n_tree$data[[1]]$sp)
## converting growth to long format and adding size class identity to each row/tree
gro.long <- gather(data.frame(gro.wide), key = interval, value = growth, -tree, -sp)
size.long <- gather(data.frame(size.class), key = interval, value = size)
## Since the census columns match exactly in wide format, they also do in long format.
## Therefore, the row ID is the same between gro.long and size.long. So no need to left_join, which takes a long time
gro.long$size <- size.long$size
# remove obs outside the defined classes
gro.long <- gro.long %>% filter(size != "out")
##--------------------------------------------------
## organised into a tibble by size
##--------------------------------------------------

by_size <- select(gro.long, -sp) %>% group_by(size) %>% nest() %>% rename(data.long = data)
names(by_size$data.long) <- by_size$size
## converting to wide format, but needs to store a row id/key (in this case tree) to gather the data later
by_size <- by_size %>%
  mutate(data.wide = map(data.long, ~ spread(.x, key = interval, value = growth)))
## but tree column is not needed for complete.cases
by_size <- by_size %>%
  mutate(data = map(data.wide, select, -tree))
## rows to keep with complete cases
by_size <- by_size %>% mutate(rows.to.keep = map(data.wide, ~select(.x, -tree) %>% complete.cases))
## subsetting data with only complete cases
by_size <- by_size %>% mutate(data.cc = map2(data, rows.to.keep, ~ .x[.y,]))
## transforming subsetted complete cases data to long format
by_size <- by_size %>% mutate(data.long.cc = map2(data.wide, rows.to.keep,
  ~ gather(.x[.y,], key = interval, value = growth, -tree)))
## no. of trees present across all intervals with at least one record
by_size <- by_size %>% mutate(n = map(data, nrow))
## no. of trees with complete cases == nrow of data.cc
by_size <- by_size %>% mutate(n.cc = map(data.cc, nrow))
## no. of tree records per interval that are non-na
count.non.na <- function (vec) {sum(!is.na(vec))}
by_size <- by_size %>% mutate(n.interval = map(data, summarise_all, count.non.na))
# To be able to select those intervals that have 5 samples per interval for a size class for at least 3 intervals
## Does each size class have at least 5 samples per interval for at least 3 intervals?
by_size <- by_size %>% mutate(enough.sample.size = map_lgl(n.interval, ~sum(.x >= 5) >= 3))
## Does each size class have at least 5 samples per interval?
by_size <- by_size %>% mutate(enough.sample.size.interval = map(n.interval, ~{.x >= 5}))
## Median growth per interval on data with complete cases
by_size <- by_size %>% mutate(med = map(data.cc, summarise_all, median, na.rm = T))
## Quantiles per interval on data with complete cases
# by_size <- by_size %>% mutate(q = map(data.c, summarise_all, quantile, probs = c(0, 0.05, 0.5, 0.95, 1), na.rm = T))
# does not work
## Function for bootstrapping medians with replacement; not used anymore###----------
# resamp.func <- function(df, iter) {
#   dfx <- setNames(data.frame(matrix(ncol = ncol(df), nrow = iter)), colnames(df))
#     for (i in 1:iter) {
#       # each census is resampled separately with replacement and taken a median of
#       dfx[i, ] <-  apply(df, 2, function (vec) {
#         ## because sample() does not work when vec[!is.na(vec)] returns 0, when all vec are NAs;
#         if (sum(!is.na(vec)) == 0) {x = vec} else {x = vec[!is.na(vec)]}
#         median(sample(x, size = floor(0.9*length(vec)), replace = T), na.rm = T)
#       }
#       )
#       }
#   return(dfx)
# }

# xx <- resamp.func(by_size$data$large, 2)
# str(xx)
# Sys.time()
# iterations = 1000
# r.gro <- by_size %>% mutate(boot.med = map(data, resamp.func, iter = iterations))
# Sys.time() #takes about 13 min
# r.gro <- r.gro %>% mutate(data.norm = map(data, apply, 1, rescale, to = c(-1, 1)))
# should I be rescaling on flattened vector?
# r.gro <- r.gro %>% mutate(boot.med.norm = map(boot.med, function(x) {xx <- data.frame(t(apply(data.frame(x), 1, rescale, to = c(-1, 1)))); return (xx)}))
# boot.med.long <- mapply(cbind, r.gro$boot.med,
#                    "iter" = rep(list("iter" = 1:nrow(r.gro$boot.med[[1]])), times = length(r.gro$boot.med)), SIMPLIFY = F) %>%
#   lapply(gather, key = census, value = growth, -iter)
# boot.med.norm.long <- mapply(cbind, r.gro$boot.med.norm,
#                         "iter" = rep(list("iter" = 1:nrow(r.gro$boot.med[[1]])), times = length(r.gro$boot.med)), SIMPLIFY = F) %>%
#   lapply(gather, key = census, value = growth, -iter)
# for (i in 1:length(gro.long)) {
#   gro.long[[i]]$interval <- recode(gro.long[[i]]$interval, !!!level_key)
#   # boot.med.norm.long[[i]]$interval <- recode(boot.med.norm.long[[i]]$interval, !!!level_key)
#
# }
# str(boot.med.long$large)
# r.gro <- r.gro %>%
#   mutate(boot.med.long = boot.med.long,
#          boot.med.norm.long = boot.med.norm.long)#-------

## To plot data, long format is used, but instead of the interval, corresponding median of the interval is needed
## for mapping time on x axis
## obtaining median time interval
time.med <- time.mat %>% summarise_all(median, na.rm = T)
interval.med <- time.med[,-ncol(time.med)] + apply(time.diff*365/2, 2, median, na.rm = T); colnames(interval.med) <- colnames(time.diff)
level_key <- list(interval.2 = interval.med$interval.2, interval.3 = interval.med$interval.3, interval.4 = interval.med$interval.4,
                  interval.5 = interval.med$interval.5, interval.6 = interval.med$interval.6, interval.7 = interval.med$interval.7)
level_key <- lapply(level_key, as.Date)
## creating new long format data in which corresponding date is substituted in place of interval
# gro.long.dates <- gro.long %>% mutate(interval = recode(interval, !!!level_key))
by_size <- by_size %>% mutate(data.long.cc.dates = map(data.long.cc,
  ~.x %>% mutate(interval = recode(interval, !!!level_key))))

by_size <- by_size %>%
  mutate(plot = imap(data.long.cc.dates,
                     ~ggplot(.x, aes(y = growth, x = interval)) +
                     ggtitle(paste0("size class = ", .y)) +
                     geom_violin(aes(group = interval), draw_quantiles = c(0.25, 0.5, 0.75), color = "gray", alpha = 0.3) +
                     geom_line(data = bind_cols(.x %>% group_by(interval) %>%
                                                  summarise_at("growth", median, na.rm = T)),
                               color = "red") +
                     ylab("Median Growth [cm/year]") +
                     xlab("Census Interval Date")
                   )
         )
print(by_size$plot$med)
grid.arrange(rectGrob(), rectGrob())
# gridExtra::grid.arrange()
m1 <- marrangeGrob(by_size$plot, nrow = 2, ncol = 2)
ggsave("figures/Growth/community_sc_growth_ts.pdf", m1, width = 12, height = 9)
# by_size$n.interval$large
# interval.2 interval.3 interval.4 interval.5 interval.6 interval.7
# 1306       1264       1405       1363       1390       1371
# $med
# interval.2 interval.3 interval.4 interval.5 interval.6 interval.7
# 6364       6531       7475       7309       7462       7283
# $small
# interval.2 interval.3 interval.4 interval.5 interval.6 interval.7
# 8595       8152       9642       9119       9551       9401
# $tiny
# interval.2 interval.3 interval.4 interval.5 interval.6 interval.7
# 151669     132252     138534     123202     126314     121473
##--------------------End---------------------------


##--------------------------------------------------
## by species and size
##--------------------------------------------------

# gro.sp <- data.frame(gro) %>% mutate(sp = n_tree$data[[1]]$sp) %>% group_by(sp) %>% nest()
# gro.sp <- gro.sp %>% mutate(med = map(data, med.func))

sp_by_size <- gro.long %>% group_by(sp, size) %>% nest() %>% rename(data.long = data)
names(sp_by_size$data.long) <- paste(sp_by_size$sp, sp_by_size$size, sep = ".")

## converting to wide format, but needs to store a row id/key (here, tree) to gather the data later
sp_by_size <- sp_by_size %>%
  mutate(data.wide = map(data.long, ~ spread(.x, key = interval, value = growth))) %>%
## but tree column is not needed for complete.cases
  mutate(data = map(data.wide, select, -tree)) %>%
## rows to keep with complete cases
  mutate(rows.to.keep = map(data.wide, ~select(.x, - tree) %>% complete.cases)) %>%
## subsetting data with only complete cases
  mutate(data.cc = map2(data, rows.to.keep, ~ .x[.y,])) %>%
## transforming subsetted complete cases data to long format for plotting
  mutate(data.long.cc = map2(data.wide, rows.to.keep,
                             ~ gather(.x[.y,], key = interval, value = growth, -tree)))

## no. of trees present across all intervals with at least one record
sp_by_size <- sp_by_size %>%
  mutate(n = map(data, nrow)) %>%
## no. of trees with complete cases == no. of observations per interval
  mutate(n.cc = map(data.cc, nrow)) %>%
## no. of tree records per interval that are non-na
  mutate(n.interval = map(data, summarise_all, function (vec) {sum(!is.na(vec))})) %>%
## to be able to select those species-size class that have at least 5 trees/interval for at least 3 intervals
## Does each sp-size class have at least 5 samples per interval for at least 3 intervals?
  mutate(enough.sample.size.cc = map_lgl(n.cc, ~{.x >= 5})) %>%
  mutate(enough.sample.size = map_lgl(n.interval, ~sum(.x >= 5) >= 3)) %>%
## Not sure how to code teh following
# ##Does each sp-size class have at least 5 samples per interval in at least three of the (5) chosen "intervals"
#   mutate(enough.sample.size.chosen.intervals = map_lgl(n.interval, ~sum(.x[, 2:6] >= 5) >= 3))  %>%
## Does each sp-size class have at least 5 samples per interval?
  mutate(enough.sample.size.interval = map(n.interval, ~{.x >= 5}))

## Median growth per interval on data without complete cases
sp_by_size <- sp_by_size %>%
  mutate(med = map(data, summarise_all, median, na.rm = T)) %>%
## Median growth per interval on data with complete cases
  mutate(med.cc = map(data.cc, summarise_all, median, na.rm = T)) %>%
## making those census intervals with insufficient sample sizes as NAs
  mutate(med.enough.sample = map2(med, enough.sample.size.interval,
                                  ~{.x[!as.vector(.y)] <- NA; return (.x)}))


## To plot data, long format is used, but instead of interval index corresponding median date of the interval is needed
## for mapping time on x axis
## obtaining median time interval
time.med <- time.mat %>% summarise_all(median, na.rm = T)
interval.med <- time.med[,-ncol(time.med)] + apply(time.diff*365/2, 2, median, na.rm = T); colnames(interval.med) <- colnames(time.diff)
level_key <- list(interval.2 = interval.med$interval.2, interval.3 = interval.med$interval.3, interval.4 = interval.med$interval.4,
                  interval.5 = interval.med$interval.5, interval.6 = interval.med$interval.6, interval.7 = interval.med$interval.7)
level_key <- lapply(level_key, as.Date)
## creating new long format data in which corresponding date is substituted in place of interval
# gro.long.dates <- gro.long %>% mutate(interval = recode(interval, !!!level_key))
sp_by_size <- sp_by_size %>% mutate(data.long.cc.dates = map(data.long.cc,
                                                       ~.x %>% mutate(interval = recode(interval, !!!level_key))))

sp_by_size <- sp_by_size %>%
  mutate(plot = map2(data.long.cc.dates, paste0("sp = ", sp_by_size$sp, ", size = ", sp_by_size$size),
                     ~ggplot(.x, aes(y = growth, x = interval)) +
                       ggtitle(.y) +
                       geom_violin(aes(group = interval), draw_quantiles = c(0.25, 0.5, 0.75), color = "gray", alpha = 0.3) +
                       geom_line(data = bind_cols(.x %>% group_by(interval) %>%
                                                    summarise_at("growth", median, na.rm = T)),
                                 color = "red") +
                       ylab("Median Growth [cm/year]") +
                       xlab("Census Interval Date")
  )
  )
## only plot sp.size class with enough sample size
print(sp_by_size$plot[[which(sp_by_size$enough.sample.size.cc)[2]]])
grid.arrange(rectGrob(), rectGrob())
large <- sp_by_size %>% filter(size == "large")
sum(large$enough.sample.size)
# 51
m1 <- marrangeGrob(large$plot[which(large$enough.sample.size.cc)], nrow = 2, ncol = 2)
ggsave("figures/Growth/sp-specific_growth_ts_large.pdf", m1, width = 12, height = 9)

small <- sp_by_size %>% filter(size == "small")
sum(small$enough.sample.size)
# 173
m2 <- marrangeGrob(small$plot[which(small$enough.sample.size.cc)], nrow = 2, ncol = 2)
ggsave("figures/Growth/sp-specific_growth_ts_small.pdf", m2, width = 12, height = 9)
##--------------------End---------------------------

##--------------------------------------------------
## Saving growth data by size or by size-sp and medians for HADAD given enough sample sizes
##--------------------------------------------------
## sp_by_size & by_size have each tree growth time series centered and scaled on 6 census intervals (including 1985-1990)
## "intervals" in the saved file indicate that only those trees were retained in the dataset that had data for at least 3 intervals in the last (here 5) "intervals"
save(sp_by_size, file = paste0("results/growth_sp_by_size_intervals_", length(census) - 1 ,".Rdata"))
save(by_size, file = paste0("results/growth_by_size_intervals_", length(census) - 1 ,".Rdata"))
load(paste0("results/growth_sp_by_size_intervals_", length(census) - 1 ,".Rdata"))
load(paste0("results/growth_by_size_intervals_", length(census) - 1 ,".Rdata"))
##--------------------------------------------------
## taking a median growth rate/interval for those species-size class that have at least 5 trees/interval for at least 3 intervals
##--------------------------------------------------
# Only retain medians for those intervals for which 5 tree/intervals, make the rest of intervals NAs
med.df <- sp_by_size %>% #filter(enough.sample.size == TRUE) %>%
  select(size, sp, med.enough.sample) %>% unnest()
## but this does not ensure that there are 5 trees/interval in the chosen "intervals" (last 5 or excluding 1985-1990)
## so excluding such species-sizes. Those will have NAs
# load("results/GLUEsetup_BCI.RData") # has model info and data on obs
# intervals = info$intervals
rows.with.sample.sizes <- apply(med.df[, (ncol(med.df) - intervals + 1): ncol(med.df)], 1, function (x) {sum(!is.na(x)) >= 3})
med_growth_sp_by_size <- med.df[rows.with.sample.sizes,]
save(med_growth_sp_by_size, file = paste0("results/tree_med_growth_sp_by_size_intervals_", intervals ,".Rdata"))
#saving trees per species
sp_by_size <- unite_(sp_by_size, "sp_size", c("sp","size"), remove = FALSE)
med_growth_sp_by_size <- unite_(med_growth_sp_by_size, "sp_size", c("sp","size"), remove = FALSE)
n <- as.numeric(sp_by_size$n[match(med_growth_sp_by_size$sp_size, sp_by_size$sp_size)])
sp.n <-  select(med_growth_sp_by_size, sp, size, sp_size) %>% mutate(n = n)
write.csv(sp.n, file.path(paste("results/tree_sp.n_med_growth_sp_by_size_", intervals ,".csv", sep = "")), row.names = TRUE)

## all sizes have at least 5 trees/interval for all intervals, so no need to filter or use med.enough.sample
med_growth_by_size <- by_size %>% select(size, med) %>% unnest()
# by_size.select <- by_size %>% select(size, med)
# med_growth_by_size <- bind_cols(size = by_size$size, map_dfr(by_size.select$med, bind_rows))
save(med_growth_by_size, file = paste0("results/tree_med_growth_by_size_intervals_", intervals ,".Rdata"))
load(paste0("results/tree_med_growth_sp_by_size_intervals_", intervals ,".Rdata"))
load(paste0("results/tree_med_growth_by_size_intervals_", intervals ,".Rdata"))
##--------------------End---------------------------

#####----------------------------------------------------------------------------------------------------
##### This is enough for now; move onto 4.Best-fit_param_BCI.R; older code below-------------
#####----------------------------------------------------------------------------------------------------

apply(dbh.wide, 2, range, na.rm = T)
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
# [1,]  300  306  309  307  312  311  313  313
# [2,]  498  571  638  732  791  852  941 1007

library(reshape2) ##--------
dbh.long <- melt(dbh.wide)
range(dbh.long$value, na.rm = T)
# [1]  300 1007
str(dbh.long)
unique(dbh.long$Var2)
ggplot(dbh.long, aes(value, colour = as.factor(Var2))) +
  geom_density() +
  # geom_jitter(aes(x = Var2, y = value), color = "gray", size = 0.5, width = 0.05) +
  my.theme + my.bg + my.adjust + theme(panel.grid.major.x = element_line(colour = "grey", size = 0.01)) +
  scale_x_continuous(limits = c(300, 800)) +
  guides(colour = guide_legend(title = "Interval")) +
  ylab(expression("Density")) + xlab("DBH(mm)")
ggsave(file.path("figures/Growth/DBH_distribution_dynamics_selection_first dbh_cc_7intervals.tiff"), height = 5, width = 9, units='in', compression = "lzw")


gl.wide_sp <- split(gl.wide, gl.wide$sp, drop = TRUE)
# Obtain median per species per census-----
sp.gl.wide_list <- lapply(gl.wide_sp, function (x) apply(x[,-1], 2, median, na.rm = T))
str(sp.gl.wide_list)
sp.gl.wide <- data.frame(do.call(rbind, sp.gl.wide_list))
str(sp.gl.wide)
head(sp.gl.wide)

## trees per species
## numbers per species, so that one can choose minimum criterion
library(dplyr)
sp.N <- summarise(group_by(gl.wide, sp), N = length(sp))
head(sp.N)
View(sp.N)
nrow(sp.N)
colnames(sp.N)[1] <- "Species"
# 48
write.csv(sp.N, file.path("results/sp.N_cc_7intervals.csv"), row.names = F)
read.csv(file.path("results/sp.N_cc_7intervals.csv"))
# ## Choosing top 10 most sbundant species within size class

top10.sp.sc <- read.csv(file.path("results/top10.sp.sc.csv"))
str(top10.sp.sc)
#top10.sp.sc <- sp.N[order(sp.N$N, decreasing = T),]$sp[1:10]
sp.5N <- subset(sp.N, N >=5)$Species
length(sp.5N)
# 19
# so there are at least 19 species with 5 trees with complete cases from 1982-2015 in size class (300,500]

sp.gl <- reshape(sp.gl.wide, direction = "long", varying = list(1:ncol(sp.gl.wide)),
                 timevar = "interval", v.names = "gt_back", ids = row.names(sp.gl.wide), idvar = "Species")
head(sp.gl)
# ## for select species
# sp.gl <-subset(sp.gl.all, Species %in% sp.5N)
# head(sp.gl)
## Obtain community median per census
dat.comm <- data.frame(gt_back = apply(gl.wide[, -1], 2, median, na.rm = T))
head(dat.comm)
str(dat.comm)
dat.comm$interval <- seq(1:nrow(dat.comm))
dat.comm$growth <- dat.comm$gt_back
dat.comm
##-------
## Back transform median gt_back into growth (g_hat)-----
sp.gl$growth <- sp.gl$gt_back
head(sp.gl)
dat.sp <- sp.gl
# dat.sp$growth <- sp.gl$growth # in mm/yr
dat.sp <- subset(dat.sp)
head(dat.sp)
str(dat.sp)
dat.sp <-left_join(dat.sp, sp.N, "Species")

## plotting at the mean of difference between all tree-wise consecutive census dates-----
# date.mat <- sapply(tree.sub1, function(x) cbind(x[, "ExactDate"]))
# head(date.mat)
# [,1]         [,2]         [,3]         [,4] [,5] [,6] [,7] [,8]
# [1,] "1981-06-12" "1985-06-28" "1990-10-20" NA   NA   NA   NA   NA
# [2,] "1982-07-01" "1985-03-06" NA           NA   NA   NA   NA   NA
# as.Date(time.mat[1,1], origin = as.Date("1960-01-01"))
# [1] "1981-06-12"
# so 1960-01-01 indeed seems to be the origin-----

t(diff(t(time.mat)))/365
head(time.diff)
jul.interval <- time.mat[,-ncol(time.mat)] + time.diff*365/2
jul.interval.mean <- apply(jul.interval, 2, mean, na.rm = T)
date.interval.mean <-  as.Date(jul.interval.mean, origin = as.Date("1960-01-01"))
library(lubridate)
interval.mean <- decimal_date(date.interval.mean)
interval.mean
write.csv(interval.mean, file.path("results/interval.mean.csv"), row.names = T)

interval.mean <- read.csv(file.path("results/interval.mean.csv"))$x

library(plyr);
##-older forms-----
# years.strp <- strptime(years, "%Y") ; years.short <- format(years.strp, "%y")
# dat.sp$censusint <- mapvalues(dat.sp$interval, from = 1:length(years[-1]), to =  paste(years[-length(years)],years[-1], sep = "-"))
# dat.sp$censusint.s <- mapvalues(dat.sp$interval, from = 1:length(years[-1]), to =  paste(years.short[-length(years)],years.short[-1], sep = "-"))
# dat.sp$censusint.m <- mapvalues(dat.sp$interval, from = 1:length(years[-1]), to =  paste(years[-length(years)],years.short[-1], sep = "-"))
##-------
dat.sp$censusint.y <- mapvalues(dat.sp$interval, from = 1:length(years[-1]), to =  interval.mean)

dat.comm$censusint.y <- mapvalues(dat.comm$interval, from = 1:length(years[-1]), to =  interval.mean)
detach("package:plyr", unload = TRUE)
head(dat.comm)

dat.sp$norm.growth <- NA
for (i in 1:length(unique(dat.sp$Species))){
  grouprows <- which(dat.sp$Species == unique(dat.sp$Species)[i])
  dat.sp$norm.growth[grouprows] <- rescale(dat.sp$growth[grouprows], to = c(-1,1))
}
head(dat.sp)
##------normalising growth of only last two censuses
dat.sp$norm.growth <- NA
for (i in 1:length(unique(dat.sp$Species))){
  grouprows <- which(dat.sp$Species == unique(dat.sp$Species)[i])
  dat.sp$norm.growth[grouprows] <- rescale(dat.sp$growth[grouprows], to = c(-1,1))
}
head(dat.sp)
##----------------------------------------------------
dat.comm$comm <- "comm"
dat.comm$norm.growth <- rescale(dat.comm$growth, to = c(-1,1))

## top10 species
top10.sp.dat <- subset(dat.sp, Species %in% top10.sp.sc$sp)
ggplot(dat.sp, aes(x = censusint.y, y = growth)) +
  geom_jitter(aes(x = censusint.y, y = growth), color = "gray", size = 0.5, width = 0.15) +
  geom_line(data = top10.sp.dat, aes(x = censusint.y, y = growth, group = Species), color = "red") +
  geom_point(data = top10.sp.dat, aes(x = censusint.y, y = growth), color = "red") +
  geom_line(data = dat.comm, aes(x = censusint.y, y = growth, group = "comm"), color = "black", size = 0.7) +
  geom_point(data = dat.comm, aes(x = censusint.y, y = growth), color = "black") +
  # geom_errorbar(aes(ymin = growth + se, ymax = growth - se), colour = "black", width = 0.01) +
  my.theme + my.bg + my.adjust + theme(panel.grid.major.x = element_line(colour = "grey", size = 0.01)) +
  scale_y_sqrt(breaks = c(0.2, 0.5, seq(1, 20)), limits = c(0, 20)) +
  scale_x_continuous(breaks = seq(1980, 2015, by = 5)) +
  ylab(expression("Growth Rate(" *mm.yr^-1*")")) + xlab("Census Interval")
ggsave(file.path("figures/Growth/Species raw growth_same_panel_selection_first dbh_cc_7intervals.tiff"), height = 5, width = 9, units='in', compression = "lzw")

ggplot(dat.sp, aes(x = censusint.y, y = norm.growth)) +
  geom_jitter(aes(x = censusint.y, y = norm.growth), color = "gray", size = 0.5, width = 0.15) +
  geom_line(data = top10.sp.dat, aes(x = censusint.y, y = norm.growth, group = Species), color = "red") +
  geom_point(data = top10.sp.dat, aes(x = censusint.y, y = norm.growth), color = "red") +
  geom_line(data = dat.comm, aes(x = censusint.y, y = norm.growth, group = "comm"), color = "black", size = 0.7) +
  geom_point(data = dat.comm, aes(x = censusint.y, y = norm.growth), color = "black") +
  scale_x_continuous(breaks = seq(1980, 2015, by = 5)) +
  # scale_y_sqrt(breaks = c(0.2, 0.5, seq(1, 20))) +
  # geom_errorbar(aes(ymin = growth + se, ymax = growth - se), colour = "black", width = 0.01) +
  my.theme + my.bg + my.adjust + theme(panel.grid.major.x = element_line(colour = "grey", size = 0.01)) +
  ylab(expression("Normalised Growth Rate (unitless)")) + xlab("Census Interval")
ggsave(file.path("figures/Growth/Species norm growth_same_panel_selection_first dbh_cc_7intervals.tiff"), height = 5, width = 9, units='in', compression = "lzw")

ggplot(top10.sp.dat, aes(x = censusint.y, y = growth)) +
  geom_line(aes(x = censusint.y, y = growth, group = Species), color = "red") +
  facet_wrap(~ Species) +
  my.theme + my.bg + my.adjust + theme(panel.grid.major.x = element_line(colour = "grey", size = 0.01)) +
  ylab(expression("Growth Rate (mm/yr)")) + xlab("Census Interval")
ggsave(file.path("figures/Growth/Top-10 Species raw growth_same_panel_selection_first dbh_cc_7intervals.tiff"), height = 5, width = 9, units='in', compression = "lzw")

ggplot(top10.sp.dat, aes(x = censusint.y, y = norm.growth)) +
  geom_line(aes(x = censusint.y, y = norm.growth, group = Species), color = "red") +
  facet_wrap(~ Species) +
  my.theme + my.bg + my.adjust + theme(panel.grid.major.x = element_line(colour = "grey", size = 0.01)) +
  ylab(expression("Normalised Growth Rate (unitless)")) + xlab("Census Interval")
ggsave(file.path("figures/Growth/Top-10 Species norm growth_same_panel_selection_first dbh_cc_7intervals.tiff"), height = 5, width = 9, units='in', compression = "lzw")

sub.dat.sp <- as.data.frame(subset(dat.sp, select = c("Species", "interval", "norm.growth")))

norm.growth <- reshape(sub.dat.sp, direction = "wide",
                       idvar = "Species", timevar = "interval", v.names = "norm.growth")
head(norm.growth)
write.csv(norm.growth, file.path("results/norm.growth_BCI_cc_7intervals.csv"), row.names = T)

sub.dat.sp2 <- as.data.frame(subset(dat.sp, select = c("Species", "interval", "growth")))

growth <- reshape(sub.dat.sp2, direction = "wide",
                  idvar = "Species", timevar = "interval", v.names = "growth")
head(growth)
write.csv(growth, file.path("results/growth_BCI_selection_cc_7intervals.csv"), row.names = T)

###-----------
## Choosing complete.cases only in last 5 censuses
## tried that and its OK (comparable to previous)
## so check
large.rows <- which(tree.sub1[[3]]$dbh >= 300 & tree.sub1[[3]]$dbh < 500) # medium(150, 100), small (80, 60)
# med.rows <- which(tree.sub1[[1]]$dbh >= 100 | tree.sub1[[1]]$dbh < 150) # medium(150, 100), small (80, 60)
# small.rows <- which(tree.sub1[[1]]$dbh >= 60 | tree.sub1[[1]]$dbh < 80) # medium(150, 100), small (80, 60)
## working on only large trees
length(large.rows)
# 1242
gl.all.wide <- gl_tree[large.rows,] ## medium(150, 100), small (80, 60)
head(gl.all.wide)

gl.all.sub <- gl.all.wide[,-c(2,3)]
head(gl.all.sub)
###-----------
gl.wide <-  gl.all.wide[complete.cases(gl.all.sub),]
nrow(gl.wide); nrow(gl.all.wide)
# [1] 451
# [1] 1242
head(gl.wide)
## how much does dbh distribution changes with time
dbh.all.wide <- dbh_low.out[large.rows,]
dbh.wide <- dbh.all.wide[complete.cases(gl.all.wide),]
apply(dbh.all.wide, 2, range, na.rm = T)
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
# [1,]   35  145  300  260   33   17   17   53
# [2,]  766  642  498  596  670  736  833  895
nrow(dbh.wide)
apply(dbh.wide, 2, median, na.rm = T)
# 303.0 327.0 353.0 373.0 391.5 407.0 417.5 431.5
nrow(dbh_low.out)
apply(dbh.wide, 2, range, na.rm = T)
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
# [1,]   82  176  300  299  300  303  305  305
# [2,]  491  485  498  568  670  736  833  895
library(reshape2)
dbh.long <- melt(dbh.wide)
range(dbh.long$value, na.rm = T)
dbh.long$Var2 <- as.factor(dbh.long$Var2)
str(dbh.long)
unique(dbh.long$Var2)

ggplot(dbh.long, aes(value, colour = as.factor(Var2))) +
  geom_density() +
  # geom_jitter(aes(x = Var2, y = value), color = "gray", size = 0.5, width = 0.05) +
  my.theme + my.bg + my.adjust + theme(panel.grid.major.x = element_line(colour = "grey", size = 0.01)) +
  scale_x_continuous(limits = c(300, 800)) +
  guides(colour = guide_legend(title = "Interval")) +
  ylab(expression("Density")) + xlab("DBH(mm)")
ggsave(file.path("figures/Growth/DBH_distribution_dynamics_selection_first dbh_cc_5intervals.tiff"), height = 5, width = 9, units='in', compression = "lzw")

gl.wide_sp <- split(gl.wide, gl.wide$sp, drop = TRUE)
head(gl.wide_sp)
str(gl.wide_sp)
# Obtain median per species per census-----
sp.gl.wide_list <- lapply(gl.wide_sp, function (x) apply(x[,-1], 2, median, na.rm = T))
str(sp.gl.wide_list)
sp.gl.wide <- data.frame(do.call(rbind, sp.gl.wide_list))
str(sp.gl.wide)
head(sp.gl.wide)

## trees per species
## numbers per species, so that one can choose minimum criterion
library(dplyr)
sp.N <- summarise(group_by(gl.wide, sp), N = n())
colnames(sp.N)[1] <- "Species"
# dont understand why sp.N.x is not a dataframe
top10.sp.sc <- read.csv(file.path("results/top10.sp.sc.csv"))
str(sp.N)
nrow(sp.N)
#57
# write.csv(sp.N, file.path("results/sp.N.csv"), row.names = T)
# ## Choosing top 10 most sbundant species within size class
# top10.sp.sc <- sp.N[order(sp.N$N, decreasing = T),]$sp[1:10]
sp.5N <- subset(sp.N, N >=5)$Species
length(sp.5N)
#29
# so there are at least 192 species with 5 trees with complete cases from 1990-2015 in size class (300,500]
## compared to 109 with complete cases from 1982-2015

sp.gl <- reshape(sp.gl.wide, direction = "long", varying = list(1:ncol(sp.gl.wide)),
                 timevar = "interval", v.names = "gt_back", ids = row.names(sp.gl.wide), idvar = "Species")
head(sp.gl)
# ## for select species
# sp.gl <-subset(sp.gl.all, Species %in% sp.5N)
# head(sp.gl)
## Obtain community median per census
dat.comm <- data.frame(gt_back = apply(gl.wide[, -1], 2, median, na.rm = T))
head(dat.comm)
str(dat.comm)
dat.comm$interval <- seq(1:nrow(dat.comm))
dat.comm$growth <- (dat.comm$gt_back)^(1/lambda)
dat.comm
##-------
## Back transform median gt_back into growth (g_hat)-----
sp.gl$growth <- (sp.gl$gt_back)^(1/lambda)
head(sp.gl)
dat.sp <- sp.gl
# dat.sp$growth <- sp.gl$growth # in mm/yr
dat.sp <- left_join(dat.sp, sp.N, by = "Species")
head(dat.sp)
str(dat.sp)
library(plyr)
interval.mean <- read.csv(file.path("results/interval.mean.csv"))$x

dat.sp$censusint.y <- mapvalues(dat.sp$interval, from = 1:length(years[-1]), to =  interval.mean)

dat.comm$censusint.y <- mapvalues(dat.comm$interval, from = 1:length(years[-1]), to =  interval.mean)
detach("package:plyr", unload = TRUE)
head(dat.comm)

choose.intervals <- unique(dat.sp$interval)[-c(1,2)]
dat.sp$norm.growth <- NA
for (i in 1:length(unique(dat.sp$Species))){
  grouprows <- which(dat.sp$Species == unique(dat.sp$Species)[i] & dat.sp$interval %in% choose.intervals)
  dat.sp$norm.growth[grouprows] <- rescale(dat.sp$growth[grouprows], to = c(-1,1))
}
head(dat.sp)

##----------------------------------------------------
dat.comm$comm <- "comm"
dat.comm$norm.growth <- NA
dat.comm$norm.growth[choose.intervals] <-   dat.sp$norm.growth[grouprows] <- rescale(dat.comm$growth[choose.intervals], to = c(-1,1))

## top10 species
top10.sp.dat <- subset(dat.sp, Species %in% top10.sp.sc$sp)
ggplot(dat.sp, aes(x = censusint.y, y = growth)) +
  geom_jitter(aes(x = censusint.y, y = growth), color = "gray", size = 0.5, width = 0.15) +
  geom_line(data = top10.sp.dat, aes(x = censusint.y, y = growth, group = Species), color = "red") +
  geom_point(data = top10.sp.dat, aes(x = censusint.y, y = growth), color = "red") +
  geom_line(data = dat.comm, aes(x = censusint.y, y = growth, group = "comm"), color = "black", size = 0.7) +
  geom_point(data = dat.comm, aes(x = censusint.y, y = growth), color = "black") +
  # geom_errorbar(aes(ymin = growth + se, ymax = growth - se), colour = "black", width = 0.01) +
  my.theme + my.bg + my.adjust + theme(panel.grid.major.x = element_line(colour = "grey", size = 0.01)) +
  scale_y_sqrt(breaks = c(0.2, 0.5, seq(1, 20)), limits = c(0, 20)) +
  scale_x_continuous(breaks = seq(1980, 2015, by = 5)) +
  ylab(expression("Growth Rate(" *mm.yr^-1*")")) + xlab("Census Interval")
ggsave(file.path("figures/Growth/Species raw growth_same_panel_selection_first dbh_cc_5intervals.tiff"), height = 5, width = 9, units='in', compression = "lzw")

ggplot(dat.sp, aes(x = censusint.y, y = norm.growth)) +
  geom_jitter(aes(x = censusint.y, y = norm.growth), color = "gray", size = 0.5, width = 0.15) +
  geom_line(data = top10.sp.dat, aes(x = censusint.y, y = norm.growth, group = Species), color = "red") +
  geom_point(data = top10.sp.dat, aes(x = censusint.y, y = norm.growth), color = "red") +
  geom_line(data = dat.comm, aes(x = censusint.y, y = norm.growth, group = "comm"), color = "black", size = 0.7) +
  geom_point(data = dat.comm, aes(x = censusint.y, y = norm.growth), color = "black") +
  scale_x_continuous(breaks = seq(1980, 2015, by = 5)) +
  # scale_y_sqrt(breaks = c(0.2, 0.5, seq(1, 20))) +
  # geom_errorbar(aes(ymin = growth + se, ymax = growth - se), colour = "black", width = 0.01) +
  my.theme + my.bg + my.adjust + theme(panel.grid.major.x = element_line(colour = "grey", size = 0.01)) +
  ylab(expression("Normalised Growth Rate (unitless)")) + xlab("Census Interval")
ggsave(file.path("figures/Growth/Species norm growth_same_panel_selection_first dbh_cc_5intervals.tiff"), height = 5, width = 9, units='in', compression = "lzw")

ggplot(top10.sp.dat, aes(x = censusint.y, y = growth)) +
  geom_line(aes(x = censusint.y, y = growth, group = Species), color = "red") +
  facet_wrap(~ Species) +
  my.theme + my.bg + my.adjust + theme(panel.grid.major.x = element_line(colour = "grey", size = 0.01)) +
  ylab(expression("Growth Rate (mm/yr)")) + xlab("Census Interval")
ggsave(file.path("figures/Growth/Top-10 Species raw growth_same_panel_selection_first dbh_cc_5intervals.tiff"), height = 5, width = 9, units='in', compression = "lzw")

ggplot(top10.sp.dat, aes(x = censusint.y, y = norm.growth)) +
  geom_line(aes(x = censusint.y, y = norm.growth, group = Species), color = "red") +
  facet_wrap(~ Species) +
  my.theme + my.bg + my.adjust + theme(panel.grid.major.x = element_line(colour = "grey", size = 0.01)) +
  ylab(expression("Normalised Growth Rate (unitless)")) + xlab("Census Interval")
ggsave(file.path("figures/Growth/Top-10 Species norm growth_same_panel_selection_first dbh_cc_5intervals.tiff"), height = 5, width = 9, units='in', compression = "lzw")

sub.dat.sp <- as.data.frame(subset(dat.sp, select = c("Species", "interval", "norm.growth")))

norm.growth <- reshape(sub.dat.sp, direction = "wide",
                       idvar = "Species", timevar = "interval", v.names = "norm.growth")
head(norm.growth)
write.csv(norm.growth, file.path("results/norm.growth_BCI_cc_5intervals.csv"), row.names = T)

sub.dat.sp2 <- as.data.frame(subset(dat.sp, select = c("Species", "interval", "growth")))

growth <- reshape(sub.dat.sp2, direction = "wide",
                  idvar = "Species", timevar = "interval", v.names = "growth")
head(growth)
write.csv(growth, file.path("results/growth_BCI_cc_5intervals.csv"), row.names = T)


