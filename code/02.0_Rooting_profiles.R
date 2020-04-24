
#-------------------------------
# Title: Generate Rooting profiles
# Author : Rutuja Chitra-Tarak
# Original date: 2019-12-20
#-------------------------------

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

# graphics info
theme_set(theme_bw())
theme_update(text = element_text(size=14),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)

file.path.figures <- file.path("figures", "rooting_profiles")
if(!dir.exists(file.path.figures)) {dir.create(file.path.figures)}
####--------------------------------------------
####      Rooting profiles for inversion
####--------------------------------------------

####--------------------------------------------
## --------------in ELM-FATES-------------------
####--------------------------------------------
# FATES uses equation (2) in Zeng et al.
## y = 1 - 0.5*(exp(-a*d) + exp(-b*d)) ....(eq 1)
# where d = depth and Y = root fraction, a = fates_roota_par, b = fates_rootb_par

# fates_roota_par: Parameter a in Eq. 2 in (Zeng 2001) that regulates the shape of the rooting profile.
# Range 5.9 - 7.4 per m
# Range corresponds to this parameter specified for BATS (or IGBP) land cover classified
# Deciduous Broadleaf Trees (5.9 m-1) and Evergreen Broadleaf Trees  (7.4 m-1) as given in
# Table 1 (or 2) of Zeng 2001.

# fates_rootb_par : Parameter b in Eq. 2 in (Zeng 2001) that regulates the depth of the rooting profile.
# Chosen range of b is derived using this equation so as to fit the observed range of rooting depth (dr)
# of 2 - 18 m for Tropical Deciduous Forest (mean ± Standard Error (S.E.);
# 3.7 ± 0.5, n = 5 trees; min = 2, max = 4.7) and
# Tropical Evergreen Forest (mean ± S.E.; 7.3 ± 0.5, n = 3 trees and 3 communities; min = 2, max = 18 m)
# combined (Canadell et al. 1996). Besides the direct observation of roots at 18 m included by (Nepstad et al. 1994)
# in Paragominas, eastern Amazonia that is included in the above study;
# in Tapajos, eastern Amazonia water extraction by roots was also inferred up to 18 m.
# (Davidson et al. Forest Science 2011)

# Rearranging eqn (1)
# exp(-a*d) + exp(-b*d) = 2*(1 - Y)
# exp(-b*d) = 2*(1 - Y) - exp(-a*d)
# b = -log(2*(1 - Y) - exp(-a*d))/d
# evaluate for Y = 0.99
# For a = 5.9 and d = 18 m
# exp(-b*18) = (1-0.99)*2 - exp(-5.9*18)
a <- 5.9; Y <- 0.99
d <- 18
b <- -log((1 - Y)*2 - exp(-a*d))/d
b
# 0.2173346
d <- 2
b <- -log((1 - Y)*2 - exp(-a*d))/d
# 1.956199
# Therefore chosen range of b: 0.2173346 - 1.956199  per m

####---------------------------------------------------------------------
##++++++++++++++ New New Profiles with Zeng et al. form +++++++++++++++++
####---------------------------------------------------------------------
 ## range of b to be within available soil depths

load(file = file.path("data-raw/psi.rda"))

soil.depths <- unique(psi$depth)
# For a = 7.4 and d = 18 m
# For d = soil.depths[13]
Y <- 0.99
a <- 5.9;
d <- soil.depths[13]
b <- -log((1 - Y)*2 - exp(-a*d))/d
b
# at Y = 0.99, b = 0.3009248; at Y = 0.95, b = 0.1771219; a <- 5.9; d <- soil.depths[13]
# at d = 1 m
d <- soil.depths[5]
a <- 19; # a needs to be exp(-a*d) - (1 - Y)*2 >= 0 for a difference of 0.01, exp(-a*d) = 0.01; a = -log(0.01)/d = 9.579037
b <- -log((1 - Y)*2 - exp(-a*d))/d
b
# at Y = 0.99, b = 15.25082, a <- 11.1, d <- soil.depths[6]
# at Y = 0.99, b = 12.46169, at Y = 0.95 b = , a <- 6.345685 , d <- soil.depths[7]
# at Y = 0.95 b = 14.64649, a <- 6.345685, d <- soil.depths[6]
#

### new profiles -----

# `a`` min, max are that from Zeng et al. for Broadleaf Deciduous and Evergreen.
# `b` min max are such that Y=95 at d =  soil.depths[7] (~0.6 m) and soil.depths[13] (~ 13 m)

all.info <- list(par.names = c("a", "b"),
                 # min.param = c(5.9, 1.303), # as in Zeng et al for broad deci b 1.955 & broad evr bmin 1.303
                 # max.param =  c(7.4, 1.956199))
                 min.param = c(5.9, 0.3009248),
                 max.param =  c(19, 30.96))
#min.param = c(5.9, 0.2173346), ## for Y = 0.99 at d = soil.depths[13], b = 0.3026635 at a = 5.9
#max.param =  c(7.4, 1.956)) # 7.4

nsam <- 200

## these should not be regenerated unless you really want to change a parameter of samples.
## because each run will create new random combinations

n.par <- length(all.info$par.names)
grid <- lhs::randomLHS(n.sam, n.par)
all.params <- matrix(0, n.sam, n.par)

## generating ensembles
for (i in 1: n.sam) {
  for (j in 1: n.par) {
    all.params[i, j] <- qunif(grid[i, j], min = all.info$min.param[j], max = all.info$max.param[j])
  }
}
ab.tag <- paste0("ab_min_", paste0(all.info$min.param, collapse = "_"), "_max_", paste0(all.info$max.param, collapse = "_"))

# rf.par <- read.csv(file = file.path(paste0("results/rf.sam_zeng_form_", ab.tag, "_", n.sam, ".csv")), header = TRUE)

rf.par <- data.frame(rf.sam = rep(1:n.sam, each = length(soil.depths)),
                     a = rep(all.params[, 1], each = length(soil.depths)),
                     b = rep(all.params[, 2], each = length(soil.depths)),
                     depth = rep(soil.depths, times = n.sam))
rf.par <- rf.par %>%
  mutate(cum.root.frac = 1 - 0.5*(exp(-a*depth) + exp(-b*depth))) %>%
  arrange(rf.sam, depth) %>%
  group_by(rf.sam) %>%
  mutate(root.frac = cum.root.frac - lag(cum.root.frac, default = NA),
         root.frac = ifelse(is.na(root.frac), 0, root.frac)) %>% ungroup(rf.sam) %>%
  arrange(rf.sam, depth)
summary(rf.par)
# max(rf.par$cum.root.frac[rf.par$depth == soil.depths[7]])
# 0.9999949 ab_min_5.9_0.3009248_max_19_30.96

# write.csv(rf.par, file.path(paste0("results/rf.sam_zeng_form_", ab.tag, "_", n.sam, ".csv")), row.names = FALSE)

g1 <- ggplot(data = rf.par, aes(y = depth, x = root.frac)) +
  xlab("Root fraction") +
  scale_y_continuous(trans="rev_sqrt", breaks = signif(round(soil.depths, 2),2)) #+
g1.1 <- g1 + geom_point(aes(group = rf.sam, color = as.factor(rf.sam)), show.legend = FALSE)
ggsave(file.path(file.path.figures, paste0("points_rooting_profiles_for_inversion_ELM_FATES_depths_zeng_", ab.tag, ".jpeg")),
       plot = g1.1, height = 6.5, width = 5, units = 'in')

g1.2 <- g1 + geom_path(aes(group = rf.sam, color = as.factor(rf.sam)), show.legend = FALSE)
ggsave(file.path(file.path.figures, paste0("path_rooting_profiles_for_inversion_ELM_FATES_depths_zeng_", ab.tag, ".jpeg")),
         plot = g1.2, height = 6.5, width = 5, units = 'in')

g2 <- ggplot(data = rf.par, aes(y = depth, x = cum.root.frac)) +
  xlab("Cumulative Root Fraction") +
  geom_line(aes(group = rf.sam, color = as.factor(rf.sam)), show.legend = FALSE) +
  # scale_colour_discrete(name = "Depth at \nCumulative Root Fraction = 0.9999999") +
  scale_y_continuous(trans="rev_sqrt", breaks = signif(round(soil.depths, 2),2)) +
  # theme(legend.position = c(0.4, 0.5)) +
  scale_x_continuous(breaks = seq(0.1, 1, by = 0.1)) #+
ggsave(file.path(file.path.figures, paste0("cumulative_rooting_profiles_for_inversion_ELM_FATES_depths_zeng_", ab.tag, ".jpeg")),
       plot = g2, height = 6.5, width = 5, units ='in')


####----------------------------------------------------------------
##++++++++++++++ End of New Profiles with Zeng et al. ++++++++++++++
####----------------------------------------------------------------

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
ggsave(file.path(file.path.figures, paste0("rooting_profiles_for_inversion_one depth.jpeg")), height = 8, width = 8, units ='in')

ggplot(profiles,
       aes(y = abs.depth, x = cum.root.frac, color = as.factor(rf.sam))) +
  ylab("Depths (m)") + xlab("Cumulative Root Fraction (0-1)") +
  geom_line(aes(group = rf.sam), alpha = 0.7, show.legend = FALSE) +
  scale_y_continuous(trans="rev_sqrt", breaks = signif(round(soil.depths, 2), 2))
ggsave(file.path(file.path.figures, paste0("rooting_profiles_for_inversion.jpeg")), height = 8, width = 8, units ='in')

ggplot(profiles %>% subset(!is.na(depth)),
       aes(y = depth, x = cum.root.frac, color = as.factor(rf.sam))) +
  ylab("Depths (m)") + xlab("Cumulative Root Fraction (0-1)") +
  geom_line(aes(group = rf.sam), alpha = 0.7, show.legend = FALSE) +geom_point(show.legend = FALSE) +
  scale_y_continuous(trans="rev_sqrt", breaks = signif(round(soil.depths, 2), 2))
ggsave(file.path(file.path.figures, paste0("rooting_profiles_for_inversion_ELM_FATES_depths.jpeg")), height = 8, width = 8, units ='in')

### only those rf.sam in which root fraction does not increase with depth

ggplot(profiles %>% subset(power < 0.5),
       aes(y = abs.depth, x = cum.root.frac, color = as.factor(rf.sam))) +
  ylab("Depths (m)") + xlab("Cumulative Root Fraction (0-1)") +
  geom_line(aes(group = rf.sam), alpha = 0.7, show.legend = FALSE) +
  scale_y_continuous(trans="rev_sqrt", breaks = signif(round(soil.depths, 2), 2))
ggsave(file.path(file.path.figures, paste0("rooting_profiles_for_inversion_power.threshold_0.5.jpeg")), height = 8, width = 8, units ='in')

ggplot(profiles %>% subset(!is.na(depth) & power < 0.5),
       aes(y = depth, x = cum.root.frac, color = as.factor(rf.sam))) +
  ylab("Depths (m)") + xlab("Cumulative Root Fraction (0-1)") +
  geom_line(aes(group = rf.sam), alpha = 0.7, show.legend = FALSE) +
  geom_point(show.legend = FALSE) +
  scale_y_continuous(trans="rev_sqrt", breaks = signif(round(soil.depths, 2), 2))
ggsave(file.path(file.path.figures, paste0("rooting_profiles_for_inversion_ELM_FATES_depths_power.threshold_0.5.jpeg")), height = 8, width = 8, units ='in')

### creating new profiles that exponentially decrease
no.of.bases <- 1000
select.x <- seq(from = 1, to = 5.6, length.out = 1000) ## so that y max is 1000
y <- round(select.x^4, 0)
plot(y ~ select.x)
head(y)

base.df <- data.frame(x = seq(from = 0, to = 500, length = no.of.bases)[y])
plot(base.df$x)
# base.df2 <- base.df %>% mutate(y = x + c(x - lag(x))/2) %>% rename(x = y)
base.df <- base.df %>% #bind_rows(base.df2) %>%
  arrange(x) %>%
  mutate(base = 1.1^(x)) %>% subset(!is.na(base) & !duplicated(base))
# base.df <- base.df[-c(1:4),]
head(base.df); nrow(base.df); summary(base.df)
  # mutate(base = seq(from = 1, to = 10^7, length.out = 100))

## for tick-marks------
#  From https://stackoverflow.com/questions/34533472/insert-blanks-into-a-vector-for-e-g-minor-tick-labels-in-r

ggplot(base.df, aes(x = x, y = base)) +
  geom_point(aes(color = base), show.legend = FALSE) +
  scale_color_distiller(name = "Base",  trans = "log10", direction = -1, palette = "RdYlBu") +
  # scale_color_viridis_c(name = "Base",  trans = "reverse", option = "inferno") +
  ylab("Base") + xlab("x") +
  scale_y_log10() +
  annotate("text", x = 100, y = 0.1*max(base.df$base, na.rm = TRUE), parse = TRUE, label =
             as.character(expression("Base = e"^x)),
            family = "Helvetica", size = 6)
ggsave(file.path(file.path.figures, paste0("base_exponetially_increasing.jpeg")), height = 6, width = 6, units ='in')

pro.df <- data.frame(depth = rep(soil.depths, times = length(base.df$base)),
                      base = rep(base.df$base, each = length(soil.depths)),
                      rf.sam = rep(1:length(base.df$base), each = length(soil.depths))) %>%
                        mutate(exp.fun = base^-depth) %>% group_by(base) %>%
  mutate(cum.exp.fun = sum(exp.fun),
         root.frac = exp.fun/cum.exp.fun,
         cum.root.frac = cumsum(root.frac))
#  max(pro.df$rf.sam) 599
# View(pro.df)
rf.breaks = signif(round(soil.depths, 2), 2)
rf.labels = rf.breaks # c(breaks[1], rep("", 4), breaks[6], "", breaks[8: length(breaks)])
p.rf <- ggplot(pro.df, aes(y = depth, x = root.frac)) +
  geom_point(aes(color = base)) +
  geom_line(aes(color = base, group = base), show.legend = FALSE) +
  ylab("Depth (m)") + xlab("Root Fraction (0-1)") +
  scale_color_distiller(name = "Base",  trans = "log10", direction = -1, palette = "RdYlBu") +
  # scale_color_distiller(name = "Base",  trans = "reverse", palette = "RdYlBu") +
  # scale_color_viridis_c(name = "Base",  trans = "reverse", option = "inferno") +
  theme(legend.position = c(0.75, 0.4)) +
  annotate("text", x = 0.37, y = 1.7, parse = TRUE, label =
             as.character(expression(RootFraction[z] == frac(Base^-Depth[z], sum(Base^-Depth[z], z==1, 13)))),
           family = "Helvetica", size = 5) +
  scale_y_continuous(trans="rev_sqrt", breaks = rf.breaks, labels = rf.labels)
  # scale_y_reverse(breaks = rf.breaks, labels = rf.labels)  # scale_y_continuous(trans="rev_sqrt", breaks = signif(round(soil.depths, 2), 2))
ggsave(file.path(file.path.figures, paste0("rooting_profiles_for_inversion_exponentially_decreasing_root.frac.jpeg")),
       plot = p.rf, height = 6, width = 4.5, units ='in')

p.cum.rf <- ggplot(pro.df, aes(y = depth, x = cum.root.frac)) +
  geom_point(aes(color = base)) +
  geom_line(aes(color = base, group = base), show.legend = FALSE) +
  ylab("Depth (m)") + xlab("Cumulative Root Fraction (0-1)") +
  scale_color_distiller(name = "Base",  trans = "log10", direction = -1, palette = "RdYlBu") +
  # scale_color_distiller(name = "Base",  trans = "reverse", palette = "RdYlBu") +
  theme(legend.position = c(0.3, 0.4)) +
  annotate("text", x = 0.5, y = 10, parse = TRUE, label =
             as.character(expression(RootFraction[z] == frac(Base^-Depth[z], sum(Base^-Depth[z], z==1, N =13)))),
           family = "Helvetica", size = 5) +
  # scale_y_reverse(breaks = rf.breaks, labels = rf.labels)
  scale_y_continuous(trans="rev_sqrt", breaks = rf.breaks, labels = rf.labels)
ggsave(file.path(file.path.figures, paste0("rooting_profiles_for_inversion_exponentially_decreasing.jpeg")),
       plot = p.cum.rf, height = 6, width = 4.5, units ='in')

write.csv(pro.df, "results/rf.sam_exponentially_decreasing.csv", row.names = FALSE)
pro.df <- read.csv(file = file.path("results/rf.sam_exponentially_decreasing.csv"), header = TRUE)

#####---------
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

