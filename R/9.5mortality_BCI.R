rm(list=ls())
gc()
# load("/Library/Frameworks/R.framework/Versions/3.4/Resources/library/CTFSRPackage/CTFSRPackage.Rdata")
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, hms)

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
## load full tree data
## only single stem present for a tree, usually the largest in dbh, status refers to the entire tree
load("data-raw/CTFScensuses/bci.tree1.Rdata")
load("data-raw/CTFScensuses/bci.tree2.Rdata")
load("data-raw/CTFScensuses/bci.tree3.Rdata")
load("data-raw/CTFScensuses/bci.tree4.Rdata")
load("data-raw/CTFScensuses/bci.tree5.Rdata")
load("data-raw/CTFScensuses/bci.tree6.Rdata")
load("data-raw/CTFScensuses/bci.tree7.Rdata")
load("data-raw/CTFScensuses/bci.tree8.Rdata")

head(bci.tree1)

census.years <- c(1982, 1985, 1990, 1995, 2000, 2005, 2010, 2015)
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
D.wide[is.na(D.wide)] <- 0
d.census2 <- D.wide[, c( -1, - 2)]
head(d.census2)
d.census1 <- D.wide[, c(-1, -length(D.wide))]
head(d.census1)
#since trees that died in the last census would be still called dead in the next census,
# to count new dead, need to remove the carry overs:
new.dead2 <- data.frame(d.census2 - d.census1)
str(new.dead2)
new.dead <- data.frame(sp = all.sp_size$sp_size, `1985` = D.wide$`1985`, new.dead2)
head(new.dead)


head(A.wide)
abund <- A.wide[, c(-1, -ncol(A.wide))]
dead <- new.dead[, -1]
surv <-  abund - dead
# changing columnnames the the year that the trees have survived:
colnames(surv) <- colnames(A.wide[, c(-1, -2)])
head(surv)
meanDate <- summarise(bci.tree %>% group_by(census),
                      t = mean(date, na.rm = T))
duration <- (meanDate$t[-1] - meanDate$t[-nrow(meanDate)])/365 # in years

mrate <- (log(abund) - log(surv))/duration
colnames(mrate) <- colnames(surv)
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

##plotting mortality rates by sp
census.years.short <- format(strptime(census.years, "%Y"), "%y")
# mrate.long$censusint <- mapvalues(mrate.long$census, from = census.years[-1], to =  paste(census.years[-length(census.years)],census.years[-1], sep = "-"))
# mrate.long$censusint.s <- mapvalues(mrate.long$census, from = census.years[-1], to =  paste(census.years.short[-length(census.years)],census.years.short[-1], sep = "-"))
mrate.long$censusint.m <- recode(mrate.long$census, `1985` = "1982-85", `1990` = "1985-90", `1995` = "1990-95", `2000` = "1995-00", `2005` = "2000-05", `2010` = "2005-10", `2015` = "2010-15")

## adding uptake depth index
load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
load(file.path("results/4.1GLUEsetup_part2_BCI.RData")) # has n.ensembles and growth and si matrix

intervals <- info$intervals
n.ensembles <- growth_by_si.info$n.ensembles
growth.type <- growth_by_si.info$growth.type
si.type <- growth_by_si.info$si.type
goodness.fit <- 0.3
##
load(file = paste("results/splevel/ds.bestfit_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", intervals, "_id.Rdata", sep = ""))
ds <- ds.bestfit
load(file = paste("results/commlevel/ds.bestfit_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", intervals, "_id.Rdata", sep = ""))
ds <- rbind(ds, ds.bestfit) %>% mutate(tlplevel = as.factor(tlplevel)) %>% subset(!is.na(udi)) %>% droplevels()

head(ds)
udi <- subset(ds, select = c("sp", "sp_size", "size", "udi", "tlplevel")) # best.type.rsq > 0.7

head(udi)
nrow(udi)
mrate.long <- left_join(mrate.long, udi, by = "sp_size") %>%
  transform(sp = reorder(sp_size, -udi))

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
require(scales)
rev_sqrt_trans <- function() {
  scales::trans_new(
    name = "rev_sqrt",
    transform = function(x) -sqrt(abs(x)),
    inverse = function(x) x^2);
}
## selecting sp with avg abundance to be greater than a threshold, say 50
n.threshold <- 10
ggplot(mrate.long %>% subset(tlplevel == "comm" &
                             avg.abund >= n.threshold),
       aes(x = udi, y = mrate, color = avg.abund)) +
  scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
                       low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
  facet_grid(size ~ censusint.m, scales = "free_y") +
  geom_point() +
  geom_smooth(method = "loess", color = "black") +
  ylab(expression("Mean Mortality Rate (per year)")) + xlab("Water Uptake Depth Index (0-1)")
ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_Mortality_rate_by_period_udi_by_size_aboveN", n.threshold, "_tlpcommn.jpeg")), height = 7, width = 12, units='in')
ggplot(mrate.long %>% subset(tlplevel == "sp"),
       aes(x = udi, y = mrate, color = avg.abund)) +
  scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
                       low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
  facet_grid(size ~ censusint.m, scales = "free_y") +
  geom_point() +
  geom_smooth(method = "loess", color = "black") +
  ylab(expression("Mean Mortality Rate (per year)")) + xlab("Water Uptake Depth Index (0-1)")
ggsave(file.path(paste0("figures/mortality/", growth.type,
                        "/sp_Mortality_rate_by_period_udi_by_size_aboveN",
                        n.threshold, "_tlpsp.jpeg")), height = 7, width = 12, units='in')

#ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_mean_Mortality_vs_udi_with_outliers.jpeg")), height = 5, width = 5, units='in')

mrate.mean <- mrate.long %>% subset(tlplevel == "sp") %>%
  group_by(sp_size) %>%
  summarize_at(vars(mrate, udi, avg.abund), mean, na.rm = TRUE) %>%
  separate(sp_size, into = c("sp", "size"), remove = FALSE) %>%
  mutate(size = factor(size, levels = c("tiny", "small", "medium", "large"))) %>%
  mutate(mrate = ifelse(!is.finite(mrate),
                         rep(NA, length(mrate)), mrate))
summary(mrate.mean)
mrate.m.1 <- lm(mrate ~ udi, data = mrate.mean %>% subset(avg.abund >= n.threshold))
summary(mrate.m.1)
summ.mrate.m.1 <- summary(mrate.m.1)
mrate.m.label.1 = paste0("y = ", round(mrate.m.1$coefficients[1], 2), " + ",
                        round(mrate.m.1$coefficients[2], 2),
                        " * x\nR-squared = ", round(summ.mrate.m.1$r.squared, 4),
                        "\np-val = ", round(summ.mrate.m.1$coefficients[2, 4], 2))

ggplot(mrate.mean %>% subset(avg.abund >= n.threshold), aes(x = udi, y = mrate)) +
  geom_point() +
  facet_wrap(. ~ size, scales = "free_y") +
  geom_smooth(method = "loess", lty = "dashed") +
  #geom_text(aes(x = 1, y = 0.20, label = mrate.m.label.1, vjust = 2)) +
  ylab(expression("Mean Mortality Rate (per year)")) + xlab("Water Uptake Depth Index (m)")
ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_mean_Mortality_vs_udi_with_outliers.jpeg")), height = 5, width = 5, units='in')


## Is this because UDI has a relationship with mean growth rate?

sp_size  <- growth_by_si.info$growth.meta$sp_size #.full %>% subset(size == "large" || size == "small") %>% select(-size, -sp, -interval.2) %>% as.matrix()
sp_size.growth.mean.named <- growth_by_si.info$growth %>% lapply(function(x) {mean(x$growth, na.rm = TRUE)}) %>%
  unlist() %>% data.frame()
sp_size.growth.mean <- data.frame(sp_size = row.names(sp_size.growth.mean.named),
                                  growth.mean = sp_size.growth.mean.named$.) %>%
  separate(sp_size, into = c("sp", "size"), remove = FALSE) %>%
  mutate(size = factor(size, levels = c("tiny", "small", "medium", "large"))) %>%
  left_join(ds %>% select(sp_size, udi, udi.se, tlplevel), by = "sp_size")

ggplot(sp_size.growth.mean %>% subset(tlplevel == "comm"), aes(x = udi, y = growth.mean)) +
  geom_point() +
  facet_wrap(. ~ size, scales = "free_y") +
  geom_smooth(method = "lm") +
  ylab(expression("Mean Growth Rate (mm/yr)")) +  xlab("Water Uptake Depth Index (m)")
ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_mean_Growth_vs_udi_size_tlpcomm.jpeg")), height = 6, width = 6, units='in')

ggplot(sp_size.growth.mean %>% subset(tlplevel == "sp"), aes(x = udi, y = growth.mean)) +
  geom_point() +
  facet_wrap(. ~ size, scales = "free_y") +
  geom_smooth(method = "lm") +
  ylab(expression("Mean Growth Rate (mm/yr)")) +  xlab("Water Uptake Depth Index (m)")
ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_mean_Growth_vs_udi_size_tlpsp.jpeg")), height = 6, width = 6, units='in')

sp  <- growth_by_si.info$growth.meta$sp
sp.growth.mean.named <- split(growth_by_si.info$growth, sp) %>% lapply(mean, na.rm = TRUE) %>% unlist() %>% data.frame()
sp.growth.mean <- data.frame(sp = row.names(sp.growth.mean.named), growth.mean = sp.growth.mean.named$.) %>%
  merge(udi, by = "sp") %>%
  transform(sp = reorder(sp, -udi))

growth.m <- lm(growth.mean ~ udi, data = sp.growth.mean %>% subset(tlplevel == "comm"))
summary(growth.m)
summ.growth.m <- summary(growth.m)
growth.m.label = paste0("y = ", round(growth.m$coefficients[1], 2), " + ",
                        round(growth.m$coefficients[2], 2),
                        " * x\nR-squared = ", round(summ.growth.m$r.squared, 4),
                        "\np-val = ", round(summ.growth.m$coefficients[2, 4], 2))

demo.sp_size <- full_join(sp_size.growth.mean, mrate.mean %>% select(sp_size, mrate, avg.abund), by = "sp_size")
head(demo.sp_size)

ggplot(demo.sp_size, aes(y = mrate.mean, x = growth.mean)) +
  geom_point() +
  facet_wrap(. ~ size, scales = "free_y") +
  geom_smooth(method = "loess") +
  #geom_text(aes(x = 0.7, y = 5, label = growth.m.label, vjust = 2)) +
  xlab(expression("Mean Growth Rate (mm/yr)")) +  ylab(expression("Mean Mortality Rate (per year)"))
ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_mean_Growth_vs_mrate.jpeg")), height = 6, width = 6, units='in')

save(demo.sp_size, file = ("results/demo.sp_size.RData"))
