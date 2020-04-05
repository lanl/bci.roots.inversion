rm(list = ls())
gc()
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, data.table)
## Tree canopy classification based on maximum diameter
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

sp.bci.tree <- split(bci.tree, bci.tree$sp)
sp.max.dbh <- lapply(sp.bci.tree, function (x) {
  as.numeric(quantile(x$dbh, 0.99, na.rm = TRUE))
})

size.lims <- data.frame(mode = c("Understorey", "Transient", "Canopy", "Large Canopy"),
                        mean.99 = c(9.8, 14.3, 27.8, 68.4),
                        se.99 = c(2.4, 9.4, 7.0, 18.5)) %>%
  mutate(upr = mean.99 + se.99,
         lwr = mean.99 - se.99,
         mode = factor(mode, levels = c("Understorey", "Transient", "Canopy", "Large Canopy")))
#           mode mean.99 se.99  upr  lwr
# 1  Understorey     9.8   2.4 12.2  7.4
# 2    Transient    14.3   9.4 23.7  4.9
# 3       Canopy    27.8   7.0 34.8 20.8
# 4 Large Canopy    68.4  18.5 86.9 49.9
ggplot(size.lims, aes(x = mode, y = mean.99)) +
  geom_point() +
  geom_errorbar(aes(ymax = upr, ymin = lwr))
cutoff <- c()
sp.mode <- sp.max.dbh %>% map_df(as_tibble) %>%
  rename(max.dbh.99pc = value) %>%
  mutate(sp = names(sp.max.dbh)) %>%
  mutate(mode = cut(max.dbh.99pc, breaks = c(0, 12.2, 23.7, 49.9, max(max.dbh.99pc, na.rm = TRUE)),
                    labels = c("Understorey", "Transient", "Canopy", "Large Canopy")))
head(sp.mode)
save(sp.mode, file = "results/sp.mode.canopy.understorey")
