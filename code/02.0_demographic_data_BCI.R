#-----------------------------------------------------
# Title: Preapring demographic data by sp & sp_size
# Author : Rutuja Chitra-Tarak
# Original date: Mar 20, 2020
#-----------------------------------------------------

rm(list=ls())
gc()
# load("/Library/Frameworks/R.framework/Versions/3.4/Resources/library/CTFSRPackage/CTFSRPackage.Rdata")
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, hms, ggpmisc)

# graphics info
theme_set(theme_bw())
theme_update(text = element_text(size=14),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)

census.years <- c(1982, 1985, 1990, 1995, 2000, 2005, 2010, 2015)

require(scales)
rev_sqrt_trans <- function() {
  scales::trans_new(
    name = "rev_sqrt",
    transform = function(x) -sqrt(abs(x)),
    inverse = function(x) x^2);
}
####**********************************************************************************************
##### ---- If you dont need to recalculate demographic rates, jump to "Load demographic data" ----
####**********************************************************************************************

## load full tree data
## only single stem present for a tree, usually the largest in dbh, status refers to the entire tree
load("data-raw/CTFScensuses/bci.stem1.Rdata")
load("data-raw/CTFScensuses/bci.tree1.Rdata")
load("data-raw/CTFScensuses/bci.tree2.Rdata")
load("data-raw/CTFScensuses/bci.tree3.Rdata")
load("data-raw/CTFScensuses/bci.tree4.Rdata")
load("data-raw/CTFScensuses/bci.tree5.Rdata")
load("data-raw/CTFScensuses/bci.tree6.Rdata")
load("data-raw/CTFScensuses/bci.tree7.Rdata")
load("data-raw/CTFScensuses/bci.tree8.Rdata")

head(bci.tree1)

bci.tree1$census <- census.years[1]
bci.tree2$census <- census.years[2]
bci.tree3$census <- census.years[3]
bci.tree4$census <- census.years[4]
bci.tree5$census <- census.years[5]
bci.tree6$census <- census.years[6]
bci.tree7$census <- census.years[7]
bci.tree8$census <- census.years[8]

bci.tree <- rbind.data.frame(bci.tree1, bci.tree2, bci.tree3, bci.tree4, bci.tree5, bci.tree6, bci.tree7, bci.tree8)
range(bci.tree$dbh, na.rm = T)

## trees that are in status dead have no size.
# for the purposes of classifying trees in size classes, they need to have dbh (size class) they originally had
## in tree data files, each row represents a tree, so only the larget stem is represented.
# When the main stem dies the dbh of the largest stem among the other stems, if present,
# is represented. So dbh can decrease over time or size class reduce.
# Need to maintain the size class of the original.

# getting a matrix of dbh by census for these stems
n_tree <- bci.tree %>% group_by(census) %>% nest()

dbh.mat <- map_dfc(n_tree$data, "dbh")
dbh.mod <- t(apply(as.matrix(dbh.mat), 1, function (x) {
  x <- as.numeric(x)
  x[is.na(x)] <- max(x, na.rm = TRUE) # it doesnt matter if NA before the tree was A (was P) turns non-NA and gets a size. It won't be counted in alive trees.
  return(x)
}))
bci.tree$dbh.mod <- as.vector(dbh.mod)
status.mat <- map_dfc(n_tree$data, "status")
cut.breaks <- c(10, 50, 100, 300, max(bci.tree$dbh.mod, na.rm = T)) # about 45 trees per interval are > 1000 mm in dbh
cut.labels <- c("tiny", "small", "medium", "large")
bci.tree <- bci.tree %>% mutate(size = cut(dbh.mod,
                           include.lowest = TRUE, breaks = cut.breaks,
                           labels = cut.labels, right = TRUE),
                           sp_size = paste(sp, size, sep = "_"))

unique(bci.tree$status)
sp_size_by_status <- summarise(bci.tree %>%
                            group_by(census, sp_size, status), n = n())
head(sp_size_by_status)
unique(sp_size_by_status$status)
### Documentation: data-raw/CTFScensuses/RoutputFull_documentation
# Status. Alive (A) and dead (D) refer to the entire tree, so if any stem is alive,
# the tree is alive, and a tree is only dead when every stem is dead. Status = ’lost stem’
# indicates that the stem had the associated code; it usually means the was broken in the
# given census, while the tree had no other stem. Status = ’missing’ (M) are cases where
# dbh and codes for a tree were not given, so it is not certain whether the tree was alive
# or dead. Status = ’prior’ (P) indicates a tree had not yet recruited at this census.
# The lost stem status is now deprecated, since it should always be safer to check stemID
# to determine whether a tree’s measurement changed stems between censuses.

# pooling all uncertain status as alive
A <- subset(sp_size_by_status, subset = status %in% c("A", "AD", "AM", "AR", "M"), select = -status)
D <- subset(sp_size_by_status, subset = status %in% c("D"), select = -status)
head(D)
library(reshape2)
A.wide <- dcast(A, sp_size ~ census, value.var = "n", fun.aggregate = sum, na.rm = TRUE)
head(A.wide)
D.wide <- dcast(D, sp_size ~ census, value.var = "n")
head(D.wide)
nrow(A.wide); nrow(D.wide)
# [1] 1032
# [1] 947
all.sp_size <- data.frame(sp_size = as.character(unique(A.wide$sp_size,
                                              D.wide$sp_size)))
D.wide <- left_join(all.sp_size, D.wide, by = "sp_size")
# ignore warning
nrow(D.wide)
# [1] 1032
head(D.wide)
d.census2 <- D.wide[, c( -1, - 2)]
head(d.census2)
d.census1 <- D.wide[, c(-1, -length(D.wide))]
head(d.census1)
#since trees that died in the last census would be still called dead in the next census,
# to count new dead, need to remove the carry-overs:
new.dead2 <- data.frame(d.census2 - d.census1)
str(new.dead2)
new.dead <- data.frame(sp_size = all.sp_size$sp_size, `1985` = D.wide$`1985`, new.dead2)
head(new.dead)
head(A.wide)
abund <- A.wide[, c(-1, -ncol(A.wide))]
dead <- new.dead[, -1]
# surv <-  abund - dead
# # changing columnnames the the year that the trees have survived:
# colnames(surv) <- colnames(A.wide[, c(-1, -2)])
# head(surv)
meanDate <- summarise(bci.tree %>% group_by(census),
                      t = mean(date, na.rm = T))
duration <- (meanDate$t[-1] - meanDate$t[-nrow(meanDate)])/365 # in years
## for sp_size:
mrate <- dead/abund*100/duration # log(abund) - log(surv))/duration
colnames(mrate) <- colnames(D.wide)[-1]
head(mrate)
# finite <- apply(apply(mrate, 2, is.finite), 1, all)
# mrate.select <- mrate[finite,]
# head(mrate.select)
## default can be 0
mrate$avg.abund <- round(rowMeans(abund), 0)
mrate$sp_size <- all.sp_size$sp_size
abund$sp_size <- all.sp_size$sp_size

abund.long <- pivot_longer(abund, cols = 1:7, names_to = "census",
                           values_to = "trees")
mrate.long <- pivot_longer(mrate, cols = 1:7, names_to = "census",
             values_to = "mrate") %>%
  left_join(abund.long, by = c("census", "sp_size")) %>%
  mutate(censusint.m = recode(census, `1985` = "1982-85", `1990` = "1985-90", `1995` = "1990-95",
                              `2000` = "1995-00", `2005` = "2000-05", `2010` = "2005-10", `2015` = "2010-15"),
         interval.num = as.numeric(recode(census, `1985` = "1",
                                          `1990` = "2", `1995` = "3",
                                          `2000` = "4", `2005` = "5",
                                          `2010` = "6", `2015` = "7")))
head(mrate.long)
save(mrate.long, file = ("results/mrate.long.RData"))

colnames(dead) <- colnames(D.wide)[-1]
dead$sp_size <- all.sp_size$sp_size
dead$abund.1982 <- abund$`1982`

dead.long <- pivot_longer(dead, cols = 1:7, names_to = "census",
                           values_to = "dead") %>%
  left_join(abund.long, by = c("census", "sp_size")) %>%
  mutate(censusint.m = recode(census, `1985` = "1982-85", `1990` = "1985-90", `1995` = "1990-95",
                              `2000` = "1995-00", `2005` = "2000-05", `2010` = "2005-10", `2015` = "2010-15"),
         interval.num = as.numeric(recode(census, `1985` = "1",
                                          `1990` = "2", `1995` = "3",
                                          `2000` = "4", `2005` = "5",
                                          `2010` = "6", `2015` = "7"))) %>%
  separate(sp_size, c("sp", "size"), remove = FALSE)
head(dead.long)
save(dead.long, file = ("results/dead.long.RData"))

sp_size.mrate.mean <- mrate.long %>%
  group_by(sp_size) %>%
  summarize_at(vars(mrate, avg.abund), mean, na.rm = TRUE) %>%
  separate(sp_size, into = c("sp", "size"), remove = FALSE) %>%
  mutate(size = factor(size, levels = c("tiny", "small", "medium", "large"))) %>%
  mutate(mrate = ifelse(!is.finite(mrate),
                        rep(NA, length(mrate)), mrate))
## for sp:
sp.abund <- A.wide %>% separate(sp_size, into = c("sp", "size"), "_") %>% select(-size) %>%
  group_by(sp) %>%
  summarise_all(list(~sum(., na.rm = TRUE)))
all.sp <- sp.abund$sp
sp.abund <- sp.abund %>% select(-sp, -`2015`)
sp.dead <- new.dead %>% separate(sp_size, into = c("sp", "size"), "_") %>% select(-size) %>%
  group_by(sp) %>%
  summarise_all(list(~sum(., na.rm = TRUE))) %>% select(-sp)
sp.mrate <- sp.dead/sp.abund*100/duration # % per year
colnames(sp.mrate) <- colnames(D.wide)[-1]
head(sp.mrate)
sp.mrate$avg.abund <- round(rowMeans(sp.abund), 0)
sp.mrate$sp <- all.sp
sp.abund$sp <- all.sp
sp.abund.long <- pivot_longer(sp.abund, cols = 1:7, names_to = "census",
                           values_to = "trees")
sp.mrate.long <- pivot_longer(sp.mrate, cols = 1:7, names_to = "census",
                           values_to = "mrate") %>%
  left_join(sp.abund.long, by = c("census", "sp")) %>%
  mutate(censusint.m = recode(census, `1985` = "1982-85", `1990` = "1985-90", `1995` = "1990-95",
                           `2000` = "1995-00", `2005` = "2000-05", `2010` = "2005-10", `2015` = "2010-15"),
         interval.num = as.numeric(recode(census, `1985` = "1",
                                          `1990` = "2", `1995` = "3",
                                          `2000` = "4", `2005` = "5",
                                          `2010` = "6", `2015` = "7")))

str(sp.mrate.long)
save(sp.mrate.long, file = ("results/sp.mrate.long.RData"))

sp.mrate.mean <- sp.mrate.long %>%
  group_by(sp) %>%
  summarize_at(vars(mrate, avg.abund), mean, na.rm = TRUE) %>%
  mutate(mrate = ifelse(!is.finite(mrate),
                        rep(NA, length(mrate)), mrate))
###------- for large
## for sp:
large.abund <- A.wide %>% separate(sp_size, into = c("sp", "size"), "_") %>%
  subset(size %in% c("medium", "large")) %>% select(-size) %>%
  group_by(sp) %>%
  summarise_all(list(~sum(., na.rm = TRUE)))
all.sp <- data.frame(sp = large.abund$sp)
large.abund <- large.abund %>% select(-sp, -`2015`)
large.dead <- new.dead %>% separate(sp_size, into = c("sp", "size"), "_") %>%
  subset(size %in% c("medium", "large")) %>% select(-size) %>%
  group_by(sp) %>%
  summarise_all(list(~sum(., na.rm = TRUE))) %>% select(-sp)
large.mrate <- large.dead/large.abund*100/duration
colnames(large.mrate) <- colnames(D.wide)[-1]
head(large.mrate)
large.mrate$avg.abund <- round(rowMeans(large.abund), 0)
large.mrate$sp <- all.sp$sp
large.abund$sp <- all.sp$sp
large.abund.long <- pivot_longer(large.abund, cols = 1:7, names_to = "census",
                                 values_to = "trees")
large.mrate.long <- pivot_longer(large.mrate, cols = 1:7, names_to = "census",
                                 values_to = "mrate") %>%
  left_join(large.abund.long, by = c("census", "sp")) %>%
  mutate(censusint.m = recode(census, `1985` = "1982-85", `1990` = "1985-90", `1995` = "1990-95",
                              `2000` = "1995-00", `2005` = "2000-05", `2010` = "2005-10", `2015` = "2010-15"),
         interval.num = as.numeric(recode(census, `1985` = "1",
                                          `1990` = "2", `1995` = "3",
                                          `2000` = "4", `2005` = "5",
                                          `2010` = "6", `2015` = "7")))

large.mrate.mean <- large.mrate.long %>%
  group_by(sp) %>%
  summarize(sd = sd(mrate, na.rm = TRUE),
            se = sd/sqrt(n()),
            mrate = mean(mrate, na.rm = TRUE),
            avg.abund = mean(mrate, na.rm = TRUE), .groups = "drop_last") %>%
  mutate(mrate = ifelse(!is.finite(mrate),
                        rep(NA, length(mrate)), mrate))

save(large.mrate.long, file = ("results/large.mrate.long.RData"))
###------- for adult
## for sp:
adult.abund <- A.wide %>% separate(sp_size, into = c("sp", "size"), "_") %>%
  subset(size %in% c("medium", "large")) %>% select(-size) %>%
  group_by(sp) %>%
  summarise_all(list(~sum(., na.rm = TRUE)))
all.sp <- data.frame(sp = adult.abund$sp)
adult.abund <- adult.abund %>% select(-sp, -`2015`)
adult.dead <- new.dead %>% separate(sp_size, into = c("sp", "size"), "_") %>%
  subset(size %in% c("medium", "large")) %>% select(-size) %>%
  group_by(sp) %>%
  summarise_all(list(~sum(., na.rm = TRUE))) %>% select(-sp)
adult.mrate <- adult.dead/adult.abund*100/duration
colnames(adult.mrate) <- colnames(D.wide)[-1]
head(adult.mrate)
adult.mrate$avg.abund <- round(rowMeans(adult.abund), 0)
adult.mrate$sp <- all.sp$sp
adult.abund$sp <- all.sp$sp
adult.abund.long <- pivot_longer(adult.abund, cols = 1:7, names_to = "census",
                              values_to = "trees")
adult.mrate.long <- pivot_longer(adult.mrate, cols = 1:7, names_to = "census",
                              values_to = "mrate") %>%
  left_join(adult.abund.long, by = c("census", "sp")) %>%
  mutate(censusint.m = recode(census, `1985` = "1982-85", `1990` = "1985-90", `1995` = "1990-95",
                              `2000` = "1995-00", `2005` = "2000-05", `2010` = "2005-10", `2015` = "2010-15"),
         interval.num = as.numeric(recode(census, `1985` = "1",
                                          `1990` = "2", `1995` = "3",
                                          `2000` = "4", `2005` = "5",
                                          `2010` = "6", `2015` = "7")))

adult.mrate.mean <- adult.mrate.long %>%
  group_by(sp) %>%
  summarize(sd = sd(mrate, na.rm = TRUE),
            se = sd/sqrt(n()),
            mrate = mean(mrate, na.rm = TRUE),
            avg.abund = mean(mrate, na.rm = TRUE), .groups = "drop_last") %>%
  mutate(mrate = ifelse(!is.finite(mrate),
                        rep(NA, length(mrate)), mrate))

save(adult.mrate.long, file = ("results/adult.mrate.long.RData"))

###----for juvenile
juve.abund <- A.wide %>% separate(sp_size, into = c("sp", "size"), "_") %>%
  subset(size %in% c("tiny", "small")) %>% select(-size) %>%
  group_by(sp) %>%
  summarise_all(list(~sum(., na.rm = TRUE)))
all.sp <- data.frame(sp = juve.abund$sp)
juve.abund <- juve.abund %>% select(-sp, -`2015`)
juve.dead <- new.dead %>% separate(sp_size, into = c("sp", "size"), "_") %>%
  subset(size %in% c("tiny", "small")) %>% select(-size) %>%
  group_by(sp) %>%
  summarise_all(list(~sum(., na.rm = TRUE))) %>% select(-sp)
juve.mrate <- juve.dead/juve.abund*100/duration
colnames(juve.mrate) <- colnames(D.wide)[-1]
head(juve.mrate)
juve.mrate$avg.abund <- round(rowMeans(juve.abund), 0)
juve.mrate$sp <- all.sp$sp
juve.abund$sp <- all.sp$sp
juve.abund.long <- pivot_longer(juve.abund, cols = 1:7, names_to = "census",
                                 values_to = "trees")
juve.mrate.long <- pivot_longer(juve.mrate, cols = 1:7, names_to = "census",
                                 values_to = "mrate") %>%
  left_join(juve.abund.long, by = c("census", "sp")) %>%
  mutate(censusint.m = recode(census, `1985` = "1982-85", `1990` = "1985-90", `1995` = "1990-95",
                              `2000` = "1995-00", `2005` = "2000-05", `2010` = "2005-10", `2015` = "2010-15"),
         interval.num = as.numeric(recode(census, `1985` = "1",
                                          `1990` = "2", `1995` = "3",
                                          `2000` = "4", `2005` = "5",
                                          `2010` = "6", `2015` = "7")))

juve.mrate.mean <- juve.mrate.long %>%
  group_by(sp) %>%
  summarize(sd = sd(mrate, na.rm = TRUE),
            se = sd/sqrt(n()),
            mrate = mean(mrate, na.rm = TRUE),
            avg.abund = mean(mrate, na.rm = TRUE), .groups = "drop_last") %>%
  mutate(mrate = ifelse(!is.finite(mrate),
                        rep(NA, length(mrate)), mrate))

save(juve.mrate.long, file = ("results/juve.mrate.long.RData"))

##---------Combinging growth rates mean across all size and at juvenile and adult level
## Only those for which avg.abundance greater than 10

sp.mrate.adult.juve <- sp.mrate.mean %>%
  full_join(juve.mrate.mean %>%
              rename(mrate.juve = mrate, m.juve.se = se, m.juve.avg.abund = avg.abund) %>% select(-sd), by = "sp") %>%
  full_join(adult.mrate.mean %>%
              rename(mrate.adult = mrate, m.adult.se = se, m.adult.avg.abund = avg.abund) %>% select(-sd), by = "sp") %>%
  full_join(large.mrate.mean %>%
              rename(mrate.large = mrate, m.large.se = se, m.large.avg.abund = avg.abund) %>% select(-sd), by = "sp")


###***********************
## growth-----------------
###***********************

growth.name <- load(file = paste0("results/sp_size.stats_growth_dbh.residuals_off_5_size_class_varying_non_cc.df.Rdata"))
growth <- get(growth.name); rm(growth.name)
large.growth.name <- load(file = paste0("results/large.stats_growth_dbh.residuals_off_5_size_class_varying_non_cc.df.Rdata"))
large.growth <- get(large.growth.name); rm(large.growth.name)
adult.growth.name <- load(file = paste0("results/adult.stats_growth_dbh.residuals_off_5_size_class_varying_non_cc.df.Rdata"))
adult.growth <- get(adult.growth.name); rm(adult.growth.name)
juve.growth.name <- load(file = paste0("results/juvenile.stats_growth_dbh.residuals_off_5_size_class_varying_non_cc.df.Rdata"))
juve.growth <- get(juve.growth.name); rm(juve.growth.name)

sp_size.growth <- growth %>%
  separate(sp_size, into = c("sp", "size"), remove = FALSE) %>%
  mutate(size = factor(size, levels = c("tiny", "small", "medium", "large")),
         grate = median,
         grate.trees = trees,
         grate.se = se,
         grate.sd = sd) %>%
  select(-median,  -mean, -n, -trees, -se, -sd)

sp.growth.large <- large.growth %>%
  mutate(grate.large = median,
         g.large.se = se,
         g.trees.large = trees) %>%
  select(sp, grate.large, g.large.se, g.trees.large)
sp.growth.adult <- adult.growth %>%
  mutate(grate.adult = median,
         g.adult.se = se,
         g.trees.adult = trees) %>%
  select(sp, grate.adult, g.adult.se, g.trees.adult)
sp.growth.juve <- juve.growth %>%
  mutate(grate.juve = median,
         g.juve.se = se,
         g.trees.juve = trees) %>%
  select(sp, grate.juve, g.juve.se, g.trees.juve)

demo.sp_size <- sp_size.mrate.mean %>%
  full_join(sp_size.growth %>%
              select(-sp, -size), by = "sp_size") %>%
  mutate(size = factor(size, levels = c("tiny", "small", "medium", "large")))
save(demo.sp_size, file = ("results/demo.sp_size.RData"))
demo.sp <- sp.growth.juve  %>%
  full_join(full_join(sp.growth.adult, sp.growth.large, by = "sp"), by = "sp") %>%
  full_join(sp.mrate.adult.juve, by = "sp")
# View(sp.demo)
save(demo.sp, file = ("results/demo.sp.RData"))
