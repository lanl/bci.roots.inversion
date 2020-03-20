## Estimating parameter b (fates_rootb_par) range for rooting profiles


rm(list=ls())

if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(ncdf4, easyNCDF, tidyverse)

require(scales);
rev_sqrt_trans <- function() {
  scales::trans_new(
    name = "rev_sqrt",
    transform = function(x) -sqrt(abs(x)),
    inverse = function(x) x^2);
}
theme_set(theme_bw())

####--------------------------------------------
#### Rooting profiles for inversion
####--------------------------------------------

load(file = file.path("data-raw/psi.rda"))

soil.depths <- unique(psi$depth)
# max.rds.df <- data.frame(max.rds = soil.depths[-length(soil.depths)]) %>%
#   mutate(lag.max.rds = lag(max.rds, default = 0),
#          diffby2 = (max.rds - lag.max.rds)/2,
#          half.depths = max.rds - diffby2)
# max.rds <- c(max.rds.df$max.rds, max.rds.df$half.depths)
by.param = 0.1
exponents1plus <- data.frame(x = seq(-6, 6, by = by.param)) %>% mutate(y = exp(x))
plot(data = exponents1plus, y~x)
# exponents.full <- rbind(exponents01, exponents1plus)
# plot(data = exponents.full, y~x)
max.rds <- soil.depths

profiles.list <- vector("list", length(max.rds))
profiles.dfs <- vector("list", length(max.rds))

exponents <- exponents1plus$y #seq(0.01, 30, by = 0.1) #c(seq(0.1, 1, length.out = 10), seq(1.2, 2.6, by = 0.2), seq(3, 15, length.out = 20))
nsam <- length(max.rds)*length(exponents)
## each rf.sam is associated with one maxD & power combination
for (i in 1: length(max.rds)) {
  maxD <- max.rds[i]
  profiles.list[[i]] <- vector("list", length(exponents))
  for (j in 1:length(exponents)) {
    rf.sam <- (i-1)*length(exponents) + j
    power <- exponents[j]
    profiles.list[[i]][[j]] <- data.frame(rf.sam = rf.sam, maxD = maxD, power = power,
                                     base = seq(0, 1, length.out = 100)) %>%
      mutate(abs.depth = base*maxD,
             cum.root.frac = base^power,
             ## associating abs.depth to a depth in soil.depths, if present
             # soil.depths[1] needs diferent treatment (signif) to find matches
             # soil.depth.1 = soil.depths[1][match(signif(round(abs.depth, 3), 3),
             #                                 signif(round(soil.depths[1], 3), 3))],
             # soil.depth.rest = soil.depths[-1][match(signif(round(abs.depth, 2), 2),
             #                                        signif(round(soil.depths[-1], 2), 2))],
             # soil.depth = if_else(is.na(soil.depth.1), soil.depth.rest, soil.depth.1))
             soil.depth = soil.depths[match(signif(round(abs.depth, 2), 2),
                                            signif(round(soil.depths, 2), 2))],
             ## only last abs.depth that matches with soil.depth is retained as matching
             # only last will ensure MAX cum.root.frac as 1
             depth = replace(soil.depth, rev(duplicated(rev(soil.depth))), NA))
  }
  profiles.dfs[[i]] <- do.call(rbind, profiles.list[[i]]) %>% as.data.frame()
}

profiles <- do.call(rbind, profiles.dfs)
write.csv(profiles, "results/root_profiles.csv", row.names = FALSE)
profiles <- read.csv(file = file.path("results/root_profiles.csv"), header = TRUE)

n.profiles <- unique(profiles$rf.sam)
length(n.profiles)
ggplot(profiles %>% subset(maxD == max.rds[length(max.rds)]), aes(y = abs.depth, x = cum.root.frac, color = maxD)) +
  ylab("Depths (m)") + xlab("Cumulative Root Fraction (0-1)") +
  geom_line(aes(group = rf.sam)) +
  scale_y_continuous(trans="rev_sqrt", breaks = signif(round(soil.depths, 2), 2))
ggsave(file.path(paste0("figures/rooting_profiles_for_inversion_one depth.jpeg")), height = 8, width = 8, units ='in')

ggplot(profiles,
       aes(y = abs.depth, x = cum.root.frac, color = as.factor(rf.sam))) +
  ylab("Depths (m)") + xlab("Cumulative Root Fraction (0-1)") +
  geom_line(aes(group = rf.sam), alpha = 0.7, show.legend = FALSE) +
  scale_y_continuous(trans="rev_sqrt", breaks = signif(round(soil.depths, 2), 2))
ggsave(file.path(paste0("figures/rooting_profiles_for_inversion.jpeg")), height = 8, width = 8, units ='in')

ggplot(profiles %>% subset(!is.na(depth)),
       aes(y = depth, x = cum.root.frac, color = as.factor(rf.sam))) +
  ylab("Depths (m)") + xlab("Cumulative Root Fraction (0-1)") +
  geom_line(aes(group = rf.sam), alpha = 0.7, show.legend = FALSE) +geom_point(show.legend = FALSE) +
  scale_y_continuous(trans="rev_sqrt", breaks = signif(round(soil.depths, 2), 2))
ggsave(file.path(paste0("figures/rooting_profiles_for_inversion_ELM_FATES_depths.jpeg")), height = 8, width = 8, units ='in')

### only those rf.sam in which root fraction does not increase with depth

ggplot(profiles %>% subset(power < 0.5),
       aes(y = abs.depth, x = cum.root.frac, color = as.factor(rf.sam))) +
  ylab("Depths (m)") + xlab("Cumulative Root Fraction (0-1)") +
  geom_line(aes(group = rf.sam), alpha = 0.7, show.legend = FALSE) +
  scale_y_continuous(trans="rev_sqrt", breaks = signif(round(soil.depths, 2), 2))
ggsave(file.path(paste0("figures/rooting_profiles_for_inversion_power.threshold_0.5.jpeg")), height = 8, width = 8, units ='in')

select.rf.sam <- unique(profiles$rf.sam[profiles$power < 0.5])
write.csv(select.rf.sam, "results/rf.sam_power.threshold_0.5.csv", row.names = FALSE)

profiles.sub <- profiles %>% select(rf.sam, depth, cum.root.frac) %>%
  subset(!is.na(depth))

from <- c(NA, soil.depths); to <- c(NA, 1:length(soil.depths))
library(plyr)
profiles <- profiles %>%
  mutate(soil.levels = plyr::mapvalues(depth, from, to),
  rf.sam.soil.levels = paste(rf.sam, soil.levels, sep = "."))
  # need to find depths that are equivalent to soil depths
  ## find which depths match which soil.depth and get those soil.depths in front of those depths
soil.depths[-length(soil.depths)]
unique(profiles$soil.depth)[-1]
unique(profiles$soil.levels)
View(profiles)
detach("package:plyr")
## but only first match would be enough
# so creating subgroup maxD.power.soildepths so that within those only first soil depth can be retained

profiles.sub1 <- profiles %>%
  subset(!is.na(soil.levels))
View(profiles.sub1)
## creating a table to work with ELM-FATES output, so all the soil depths present there are needed
## *for each rf.sam* and a value for root fraction against it,
## this is to be made for each nsam above.
## cum.root.frac = 0 for depths before teh depth at which cum.root.frac record begin
# cum.root.frac =  1 for depths after the depth cum.root.frac reaches 1

pro.df <- data.frame(rf.sam = rep(1:nsam, each = length(soil.depths)),
                     depth = rep(soil.depths, times = nsam),
                     soil.levels = rep(1:length(soil.depths), times = nsam)) %>%
  mutate(rf.sam.soil.levels = paste(rf.sam, soil.levels, sep = ".")) %>%
  left_join(profiles.sub1 %>% select(rf.sam.soil.levels, cum.root.frac),
            by = "rf.sam.soil.levels")
pro.df.filled <- pro.df %>%
  mutate(cum.root.frac = replace_na(cum.root.frac, 0))

## remove columns that vary within a rf.sam
cum.root.pro <- pro.df %>% select(-depth, -rf.sam.soil.levels) %>%
  pivot_wider(names_from = soil.levels, names_prefix = "level.", values_from = cum.root.frac)
cum.root.mat <- cum.root.pro %>% select(paste0("level.", 1:length(soil.depths)))

## cum.root.frac = 0 only for depths before the depth at which cum.root.frac record begin
## cum.root.frac =  1 for depths after the depth cum.root.frac reaches 1
ncol.mat = ncol(cum.root.mat)
cum.root.mat.filled <- cum.root.mat
for (i in 1: nrow(cum.root.mat)) {
    if (is.na(cum.root.mat.filled[i, 1])) {
      cum.root.mat.filled[i, 1] <- 0
    }
  }
for (i in 1: nrow(cum.root.mat)) {
  for (j in 2: ncol.mat){
    if (is.na(cum.root.mat.filled[i, j])) {
      cum.root.mat.filled[i, j] <- cum.root.mat.filled[i, j - 1]
    }
  }
}

cum.pro.df <- cbind(cum.root.pro$rf.sam, cum.root.mat.filled)
colnames(cum.pro.df) <- colnames(cum.root.pro)


pro.df <- cum.pro.df %>% data.frame() %>%
  group_by(rf.sam) %>% pivot_longer(cols = starts_with("level."),
                                    names_to = "soil.levels",
                                    names_prefix = "level.",
                                    values_to = "cum.root.frac") %>%
  mutate(root.frac = cum.root.frac - lag(cum.root.frac, default = 0)) %>%
  ## after all soil layers with root.frac, at next soil layer, cum.root.frac should be 1, and 0 for the rest
  ungroup(rf.sam) %>%
  mutate(root.frac = if_else(root.frac == -1.0000000000, 0, root.frac))
View(pro.df)
library(plyr)
pro.df <- pro.df %>% mutate(soil.levels = as.numeric(soil.levels),
                            depth = round(plyr::mapvalues(soil.levels, to[-1], from[-1]), 2))

head(pro.df)
detach("package:plyr")
write.csv(pro.df, "results/root.profiles.long.csv", row.names = FALSE)

head(pro.df)
View(pro.df)
root.pro <- pro.df %>% select(-soil.levels, -cum.root.frac) %>%
  pivot_wider(names_from = depth, names_prefix = "depth.", values_from = root.frac)
head(root.pro)
write.csv(root.pro, "results/root.profiles.wide.csv", row.names = FALSE)
## nsam = 533
####---End--Rooting Profiles for inversion-------------------------------------------

