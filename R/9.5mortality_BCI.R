rm(list=ls())
gc()
# load("/Library/Frameworks/R.framework/Versions/3.4/Resources/library/CTFSRPackage/CTFSRPackage.Rdata")
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, hms, ggpmisc)

# graphics info
tex <- element_text(size = 16, face = "plain") # , family = "gara"
tex <- element_text(size = 5, face = "plain")
my.theme <-  theme(axis.text = tex, axis.title = tex,
                   title = tex, legend.title = tex, legend.text = tex, strip.text.y = tex, strip.text.x = tex)
my.bg <- theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "black", size = 0.25),
                            panel.grid.minor = element_blank())
my.adjust <- theme(axis.title.y = element_text(vjust = 1), axis.title.x = element_text(vjust = -0.6), title = element_text(vjust = 2))
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

mrate.long <- pivot_longer(mrate, cols = 1:7, names_to = "census",
             values_to = "mrate")
head(mrate.long)
save(mrate.long, file = ("results/mrate.long.RData"))

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
sp.mrate.long <- pivot_longer(sp.mrate, cols = 1:7, names_to = "census",
                           values_to = "mrate")
str(sp.mrate.long)
save(sp.mrate.long, file = ("results/sp.mrate.long.RData"))

sp.mrate.mean <- sp.mrate.long %>%
  group_by(sp) %>%
  summarize_at(vars(mrate, avg.abund), mean, na.rm = TRUE) %>%
  mutate(mrate = ifelse(!is.finite(mrate),
                        rep(NA, length(mrate)), mrate))
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
adult.mrate.long <- pivot_longer(adult.mrate, cols = 1:7, names_to = "census",
                              values_to = "mrate")
head(adult.mrate.long)
adult.mrate.mean <- adult.mrate.long %>%
  group_by(sp) %>%
  summarize_at(vars(mrate, avg.abund), mean, na.rm = TRUE) %>%
  mutate(mrate = ifelse(!is.finite(mrate),
                        rep(NA, length(mrate)), mrate))
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
juve.mrate.long <- pivot_longer(juve.mrate, cols = 1:7, names_to = "census",
                                 values_to = "mrate")
head(juve.mrate.long)
juve.mrate.mean <- juve.mrate.long %>%
  group_by(sp) %>%
  summarize_at(vars(mrate, avg.abund), mean, na.rm = TRUE) %>%
  mutate(mrate = ifelse(!is.finite(mrate),
                        rep(NA, length(mrate)), mrate))
##---------Combinging growth rates mean acros all size and at juvenile and adult level
## Only those for which avg.abundance greater than 10
sp.mrate.adult.juve <- sp.mrate.mean %>% subset(avg.abund >= 10) %>%
  full_join(juve.mrate.mean %>% subset(avg.abund >= 10) %>% rename(mrate.juve = mrate) %>% select(-avg.abund), by = "sp") %>%
  full_join(adult.mrate.mean %>% subset(avg.abund >= 10) %>% rename(mrate.adult = mrate)%>% select(-avg.abund), by = "sp")

###***********************
## growth-----------------
###***********************
growth.name <-load(file = paste0("results/sp_size.med_growth_dbh.residuals_off_5_size_class_varying_non_cc.Rdata"))
growth <- get(growth.name); rm(growth.name)
sp_size  <- names(growth) #.full %>% subset(size == "large" || size == "small") %>% select(-size, -sp, -interval.2) %>% as.matrix()
sp_size.growth.mean.named <- growth %>% lapply(function(x) {mean(x$growth, na.rm = TRUE)}) %>%
  unlist() %>% data.frame()
sp_size.growth.mean <- data.frame(sp_size = row.names(sp_size.growth.mean.named),
                                  grate = sp_size.growth.mean.named$.) %>%
  separate(sp_size, into = c("sp", "size"), remove = FALSE) %>%
  mutate(size = factor(size, levels = c("tiny", "small", "medium", "large")))
sp.growth <- sp_size.growth.mean %>% group_by(sp) %>%
  summarise(grate = mean(grate, na.rm = TRUE))
sp.growth.adult <- sp_size.growth.mean %>% subset(size %in% c("medium", "large")) %>%
  group_by(sp) %>%
  summarise(grate.adult = mean(grate, na.rm = TRUE))
sp.growth.juvenile <- sp_size.growth.mean %>% subset(size %in% c("tiny", "small")) %>%
  group_by(sp) %>%
  summarise(grate.juve = mean(grate, na.rm = TRUE))
demo.sp_size <- sp_size.mrate.mean %>%
  full_join(sp_size.growth.mean %>% select(-sp, -size), by = "sp_size") %>%
  mutate(size = factor(size, levels = c("tiny", "small", "medium", "large")))
save(demo.sp_size, file = ("results/demo.sp_size.RData"))
demo.sp <- sp.growth %>% full_join(full_join(sp.growth.adult, sp.growth.juvenile, by = "sp"), by = "sp") %>%
  full_join(sp.mrate.adult.juve, by = "sp")
# View(sp.demo)
save(demo.sp, file = ("results/demo.sp.RData"))

####********************************************************************
##### Load demographic data ----------------------------------------
##### If you dont need to update the above, start directly here:
####********************************************************************

load("results/demo.sp.RData")
load("results/demo.sp_size.RData")
load("results/mrate.long.RData")
load("results/sp.mrate.long.RData")

##plotting mortality rates by sp
census.years.short <- format(strptime(census.years, "%Y"), "%y")
# mrate.long$censusint <- mapvalues(mrate.long$census, from = census.years[-1], to =  paste(census.years[-length(census.years)],census.years[-1], sep = "-"))
# mrate.long$censusint.s <- mapvalues(mrate.long$census, from = census.years[-1], to =  paste(census.years.short[-length(census.years)],census.years.short[-1], sep = "-"))
mrate.long$censusint.m <- recode(mrate.long$census, `1985` = "1982-85", `1990` = "1985-90", `1995` = "1990-95", `2000` = "1995-00", `2005` = "2000-05", `2010` = "2005-10", `2015` = "2010-15")
sp.mrate.long$censusint.m <- recode(sp.mrate.long$census, `1985` = "1982-85", `1990` = "1985-90", `1995` = "1990-95", `2000` = "1995-00", `2005` = "2000-05", `2010` = "2005-10", `2015` = "2010-15")

## adding uptake depth index
load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
load(file.path("results/4.1GLUEsetup_part2_BCI.RData")) # has n.ensembles and growth and si matrix

intervals <- info$intervals
n.ensembles <- growth_by_si.info$n.ensembles
growth.type <- growth_by_si.info$growth.type
growth.selection <- growth_by_si.info$growth.selection
si.type <- growth_by_si.info$si.type
goodness.fit <- 0.3
soil.depths <- unique(info$root.param.long$depth)
##
##
dbh.residuals <- "on"#growth_by_si.info$dbh.residuals
dryseason <- "1992"
root.selection <- "on"
iso.subset <- "off"
##
file.extension.base4 <- paste0(goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection,
                               "_", dbh.residuals, "_", intervals, "_dryseason_", dryseason, "_iso.subset_",
                               iso.subset, "_root.selection_", root.selection)

level.folder = "splevel"

if(!dir.exists(paste0("figures/mortality/", growth.type, "/", level.folder))) {dir.create(paste0("figures/mortality/", growth.type, "/", level.folder))}
load(file = paste0("results/", level.folder, "/ds.bestfit_cor", file.extension.base4, ".Rdata"))
ds <- ds.bestfit
# load(file = paste("results/commlevel/ds.bestfit_cor", file.extension.base4, ".Rdata", sep = ""))
# ds <- rbind(ds, ds.bestfit)
ds <- ds %>% mutate(tlplevel = as.factor(tlplevel)) %>% subset(!is.na(udi)) %>% droplevels()

head(ds)
udi <- subset(ds, select = c("sp", "sp_size", "size", "udi.best", "tlplevel")) # best.type.rsq > 0.7

head(udi)
nrow(udi)
head(sp.mrate.long)
sp.mrate.long <- left_join(sp.mrate.long, subset(udi, size == "large") %>% select(sp, udi.best), by = "sp") %>%
  mutate(sp = as.character(sp)) %>%
  transform(sp = reorder(sp, -udi.best))
mrate.long <- left_join(mrate.long, select(udi, sp_size, udi.best), by = "sp_size") %>%
  mutate(sp_size = as.character(sp_size)) %>%
  transform(sp_size = reorder(sp_size, -udi.best)) %>%
  separate(sp_size, c("sp", "size", sep = "_"), remove = FALSE) %>% select(-"_")

# ggplot(mrate.long, aes(x = censusint.m, y = mrate)) +
#   geom_line(aes(x = censusint.m, y = mrate, group = sp, color = sp), size = 1, show.legend = F) +
#   my.theme + my.bg + my.adjust + theme(panel.grid.major.x = element_line(colour = "grey", size = 0.01)) +
#   ylab(expression("Mortality Rate")) + xlab("Census Interval")
# ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_Mortality_rate_above", dbhthreshold/10, "cmDBH_aboveN50.jpeg")), height = 5, width = 9, units='in')
#
# ggplot(mrate.long, aes(x = censusint.m, y = mrate)) +
#   geom_line(aes(x = censusint.m, y = mrate, group = sp, color = sp), size = 1, show.legend = F) +
#   geom_point(aes(x = censusint.m, y = mrate, color = sp), size = 3, show.legend = F) +
#   facet_wrap( ~ sp, scales = "free_y") +
#   theme(axis.text.y = element_text(size = 8),
#         axis.text.x = element_text(size = 6, face = "plain", angle = 45, vjust = 1, hjust = 1)) +
#   ylab(expression("Mortality Rate")) + xlab("Census Interval")
# ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_Mortality_rate_by_period_above", dbhthreshold/10, "cmDBH_aboveN50.jpeg")), height = 12, width = 15, units='in')

## selecting sp with avg abundance to be greater than a threshold, say 50
n.threshold <- 10
# ggplot(mrate.long %>% subset(tlplevel == "comm" &
#                              avg.abund >= n.threshold),
#        aes(x = udi.best, y = mrate, color = avg.abund)) +
#   scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
#                        low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
#   facet_grid(size ~ censusint.m, scales = "free_y") +
#   geom_point() +
#   geom_smooth(method = "loess", color = "black") +
#   ylab(expression("Mean Mortality Rate (per year)")) + xlab("Water Uptake Depth Index (m)")
# ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_Mortality_rate_by_period_udi_by_size_aboveN", n.threshold, "_tlpcommn.jpeg")), height = 7, width = 12, units='in')
## restricting analysis to only canopy and large canopy species
bci.traits <- read.csv("data-raw/traits/BCITRAITS_20101220.csv") %>%
  rename(form1 = GRWFRM1., form2 = GRWFRM2., sp = SP.) %>% mutate(sp = tolower(sp))

# load(file = "results/sp.mode.canopy.understorey")
canopy.bci.traits <- bci.traits %>% subset(form1 %in% c("T")) %>% droplevels()
sp.select <- unique(canopy.bci.traits$sp)
ggplot(mrate.long %>% subset(!is.na(size) & !is.na(udi.best)),
       aes(x = udi.best, y = mrate, color = avg.abund)) +
  scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
                       low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
  facet_grid(size ~ censusint.m, scales = "free_y") +
  geom_point() +
  geom_smooth(method = "loess", color = "black") +
  scale_x_continuous(trans="sqrt", breaks = c(soil.depths[1], signif(soil.depths[c(8, 10, 12)], 1))) +
  ylab(expression("Mean Mortality Rate (% per year)")) + xlab("Water Uptake Depth Index (m)")
ggsave(file.path(paste0("figures/mortality/", growth.type, "/", level.folder,
                        "/sp_Mortality_rate_by_period_udi_by_size_aboveN",
                        n.threshold, "_tlpsp.jpeg")), height = 7, width = 13, units='in')
ggplot(sp.mrate.long %>% subset(!is.na(udi.best)),
       aes(x = udi.best, y = mrate, color = udi.best)) +
  scale_color_gradient(name = "UDI", trans = "rev_sqrt",
                       low = "red", high = "blue") +
  # scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
  #                      low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
  facet_grid(. ~ censusint.m) +
  geom_point() +
  geom_smooth(method = "loess", color = "black") +
  scale_x_continuous(trans="sqrt", breaks = c(soil.depths[1], signif(soil.depths[c(8, 10, 12)], 1))) +
  ylab(expression("Mean Mortality Rate (% per year)")) + xlab("Water Uptake Depth Index (m)")
ggsave(file.path(paste0("figures/mortality/", growth.type, "/", level.folder,
                        "/sp_sp.Mortality_rate_by_period_udi_aboveN",
                        n.threshold, "_tlpsp.jpeg")), height = 3, width = 13, units='in')

ggplot(sp.mrate.long %>% subset(!is.na(udi.best) & avg.abund >= n.threshold),
       aes(y = udi.best, x = mrate, color = udi.best)) +
  scale_color_gradient(name = "UDI", trans = "rev_sqrt",
                       low = "red", high = "blue", breaks = c(0, 0.5, 2.5, 5, 7.5)) +
  facet_grid(censusint.m ~ .) +
  geom_point() +
  geom_smooth(method = "loess", color = "black") +
  scale_y_continuous(trans="rev_sqrt",
                     breaks = c(round(soil.depths[c(8, 10, 12)], 0))) +
  xlab(expression("Mean Mortality Rate (% per year)")) + ylab("Water Uptake Depth Index (m)")
ggsave(file.path(paste0("figures/mortality/", growth.type, "/", level.folder,
                        "/sp_sp.Mortality_rate_by_period_udi_aboveN",
                        n.threshold, "_tlpsp_vertical.jpeg")), height = 13, width = 5, units='in')
formula = y ~ x
ggplot(sp.mrate.long %>% subset(!is.na(udi.best)),
       aes(y = udi.best, x = mrate, color = udi.best)) +
  # scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
  #                      low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
  scale_color_gradient(name = "UDI", trans = "rev_sqrt",
                       low = "red", high = "blue", breaks = c(0.5, 2.5, 5, 7.5, 10)) +
  facet_grid(censusint.m ~ .) +
  geom_point() +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.1, npcy = 0.1, rr.digits = 2,
               formula = formula, parse = TRUE, size = 3) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                  npcx = 0.75, npcy = 0.1, size = 3) +
  # scale_x_log10() +
  scale_y_continuous(trans="rev_sqrt", breaks = c(soil.depths[1], signif(soil.depths[c(8, 10, 12)], 1))) +
  xlab(expression("Mean Mortality Rate (% per year)")) + ylab("Water Uptake Depth Index (m)")
ggsave(file.path(paste0("figures/mortality/", growth.type,  "/", level.folder,
                        "/sp_sp.Mortality_rate_by_period_udi_aboveN",
                        n.threshold, "_tlpsp_vertical_lm.jpeg")), height = 13, width = 5, units='in')

sp.mrate.long.class <- sp.mrate.long %>% subset(avg.abund >= n.threshold) %>%
  mutate(udi.best.class = forcats::fct_explicit_na(cut(udi.best, breaks = c(0, 2.5, 5, 11)))) %>%
  group_by(udi.best.class, censusint.m) %>%
  summarise(mean.mrate = mean(mrate, na.rm = TRUE),
            se = sd(mrate, na.rm = TRUE)/sqrt(n()))
# udi.best.class = factor(udi.best.class, levels = c("(5,11]", "(2.5,5]", "(0,2.5]", "(Missing)"
ggplot(sp.mrate.long.class %>% subset(!is.na(udi.best.class)),
       aes(y = udi.best.class, x = mean.mrate)) +
  # scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
  #                      low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
  # scale_color_gradient(name = "UDI", trans = "rev_sqrt",
  #                      low = "red", high = "blue") +
  facet_grid(censusint.m ~ .) +
  geom_point() +
  geom_errorbarh(aes(xmax = mean.mrate + se, xmin = mean.mrate - se), size = 0.5) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  # scale_x_log10() +
  # scale_y_continuous(trans="rev_sqrt", breaks = c(soil.depths[1], signif(soil.depths[c(8, 10, 12)], 1))) +
  xlab(expression("Mean Mortality Rate (% per year)")) + ylab("Water Uptake Depth Index (m)")
ggsave(file.path(paste0("figures/mortality/", growth.type,  "/", level.folder,
                        "/sp_sp.Mortality_rate_by_period_udi_aboveN",
                        n.threshold, "_tlpsp_UDI.class.jpeg")), height = 13, width = 5, units='in')

## plotting diff drought-induced - background:
sp.mrate.long.bck <- sp.mrate.long %>%
  subset(censusint.m %in% c("1990-95", "1995-00") & sp %in% sp.select & !is.na(udi.best)) %>%
  mutate(sp = as.character(sp)) %>%
  group_by(avg.abund, sp, udi.best) %>%
  summarise(mrate = mean(mrate, na.rm = TRUE)) %>%
  mutate(period = "1990-2000")
sp.mrate.long.drought <- sp.mrate.long %>%
  subset(censusint.m %in% c("2000-05", "2005-10", "2010-15") & sp %in% sp.select & !is.na(udi.best)) %>%
  group_by(avg.abund, sp, udi.best) %>%
  summarise(mrate = mean(mrate, na.rm = TRUE)) %>%
  mutate(period = "2000-2015")
sp.mrate.long.diff <- bind_rows(sp.mrate.long.drought, sp.mrate.long.bck) %>% arrange(sp, period) %>%
  group_by(sp) %>%
  mutate(mrate.diff = mrate - lag(mrate)) %>% ungroup(sp)
sp.mrate.long.diff
ggplot(sp.mrate.long.diff, #  %>% subset(udi.best < 0.8)
       aes(y = udi.best, x = mrate.diff, color = udi.best)) +
  scale_color_gradient(name = "UDI", trans = "rev_sqrt",
                       low = "red", high = "blue", breaks = c(0, 0.5, 2.5, 5, 7.5)) +
  geom_point() +
  geom_vline(xintercept = 0) +
  geom_smooth(method = "lm", color = "black") +
  scale_y_continuous(trans="rev_sqrt",
                     breaks = c(0.01, signif(soil.depths[c(8, 10, 12)], 1))) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.1, npcy = 0.1, rr.digits = 2,
               formula = formula, parse = TRUE, size = 3) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                  npcx = 0.65, npcy = 0.1, size = 3) +
  theme(axis.title.x = element_text(vjust = -2)) +
  xlab(expression("2000-2015 - 1990-2000\nDifference Mean Mortality Rate (% per year)")) + ylab("Water Uptake Depth Index (m)")
ggsave(file.path(paste0("figures/mortality/", growth.type,  "/", level.folder,
                        "/sp_sp.Mortality_rate_by_period_udi_aboveN",
                        n.threshold, "_tlpsp_diff_mean_2000-2015&1990-2000.jpeg")), height = 5, width = 5.5, units='in')
#ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_mean_Mortality_vs_udi_with_outliers.jpeg")), height = 5, width = 5, units='in')


##### Is this because UDI has a relationship with mean growth rate?----

demo.sp_size.udi <- left_join(demo.sp_size, select(udi, sp_size, udi.best), by = "sp_size") %>%
  transform(sp_size = reorder(sp_size, -udi.best)) %>%
  separate(sp_size, c("sp", "size", sep = "_"), remove = FALSE) %>%
  select(-"_") %>% subset(!size == "NA") %>% droplevels() %>%
  mutate(size = factor(size, levels = c("tiny", "small", "medium", "large")))

ggplot(demo.sp_size.udi, aes(x = udi.best, y = grate)) +
  geom_point() +
  facet_wrap(. ~ size, scales = "free_y") +
  geom_smooth(method = "lm") +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.1, npcy = 0.8, rr.digits = 2,
               formula = formula, parse = TRUE, size = 3) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                  npcx = 0.65, npcy = 0.79, size = 3) +
  ylab(expression("Mean Growth Rate (mm/yr)")) +  xlab("Water Uptake Depth Index (m)")
ggsave(file.path(paste0("figures/mortality/", growth.type, "/", level.folder, "/sp_mean_Growth_vs_udi_size.jpeg")), height = 6, width = 6, units='in')

ggplot(demo.sp_size.udi %>% subset(avg.abund >= n.threshold), aes(x = udi.best, y = mrate)) +
  geom_point() +
  facet_wrap(. ~ size, scales = "free_y") +
  geom_smooth(method = "loess", lty = "dashed") +
  scale_x_continuous(trans="sqrt", breaks = soil.depths[-c(2,3, 4, 6, 7, 9)]) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.1, npcy = 0.8, rr.digits = 2,
               formula = formula, parse = TRUE, size = 3) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                  npcx = 0.65, npcy = 0.78, size = 3) +
  ylab(expression("Mean Mortality Rate (% per year)")) + xlab("Water Uptake Depth Index (m)")
ggsave(file.path(paste0("figures/mortality/", growth.type, "/", level.folder, "/sp_mean_Mortality_vs_udi_with_outliers_avg.abund_above", n.threshold, ".jpeg")), height = 5, width = 6, units='in')

ggplot(demo.sp_size %>% subset(size != "NA"), aes(y = mrate, x = grate)) +
  geom_point() +
  facet_wrap(. ~ size, scales = "free_y") +
  geom_smooth(method = "lm") +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.1, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 3) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                  npcx = 0.65, npcy = 0.88, size = 3) +
  xlab(expression("Mean Growth Rate (mm/yr)")) +  ylab(expression("Mean Mortality Rate (% per year)"))
ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_mean_Growth_vs_mrate.jpeg")), height = 6, width = 6, units='in')

