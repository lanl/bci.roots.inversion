###------------------------------
### Find out distribution of growth rates and residuals of water-stress on growth
### Author: Rutuja
### Date: April 16, 2020
###------------------------------
rm(list = ls())

gc()
if (!require("pacman")) install.packages("pacman")
p_load(tidyverse, gridExtra)

### Growth rates of dbh residuals

load(file.path("results/4.1GLUEsetup_part2_BCI.RData"))
growth.rate <- as.vector(unlist(growth_by_si.info$growth))
qex <- function(x) qexp((rank(x)-.375)/(length(x)+.25))

jpeg(paste0("figures/growth_rates_distribution_dbh.residuals_", growth_by_si.info$dbh.residuals, ".jpeg"),
     width = 960, height = 960, units = "px", pointsize = 24,
     quality = 100)
par(mfrow = c(2, 2))
qqnorm(growth.rate, main = "Normal Q-Q plot: Growth Rate")
plot(qex(growth.rate), growth.rate)
qqnorm(log(growth.rate), main = "Normal Q-Q plot: log(Growth Rate)")
plot(qex(growth.rate), log(growth.rate))
dev.off()

load("/Users/rutuja/Work_at_LANL/Projects/bci.roots.inversion/results/splevel/
     GLUE.resid_drop.monthsNone_cor0.3_2019-10-14_5000_0_med_size_class_predefined_cc_scaled_on_5.Rdata")
resid.growth.all <- as.numeric(unlist(GLUE.resid))
resid.growth <- sample(resid.growth.all, 10000)
jpeg("figures/residual_growth_rates_distribution.jpeg", width = 960, height = 960, units = "px", pointsize = 24,
     quality = 100)
par(mfrow = c(2, 2))
qqnorm(resid.growth, main = "Normal Q-Q plot: Residuals")
plot(qex(resid.growth), resid.growth)
qqnorm(log(resid.growth), main = "Normal Q-Q plot: log(Residuals)")
plot(qex(resid.growth), log(resid.growth))
dev.off()

