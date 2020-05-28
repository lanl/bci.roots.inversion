
###-------------------------------
### Find out distribution of growth rates, residuals of inverse model, UDI
### Author: Rutuja
### Original Date: April 16, 2020
###-------------------------------

rm(list = ls())

gc()
if (!require("pacman")) install.packages("pacman")
p_load(tidyverse, gridExtra, bci.elm.fates.hydro, spatstat, MASS)
figures.folder <- paste0("figures/Distributions")
if(!dir.exists(file.path(figures.folder))) {dir.create(file.path(figures.folder))}

#*****************
### Growth rates of dbh residuals ------
#*****************

load(file.path("results/GLUEsetup_part2_BCI.RData"))
growth.rate <- as.vector(unlist(growth_by_si.info$growth))
qex <- function(x) qexp((rank(x)-.375)/(length(x)+.25))

jpeg(file.path(figures.folder, "growth_rates_distribution_dbh.residuals_", growth_by_si.info$dbh.residuals, ".jpeg"),
     width = 960, height = 960, units = "px", pointsize = 24,
     quality = 100)
par(mfrow = c(2, 2))
qqnorm(growth.rate, main = "Normal Q-Q plot: Growth Residuals")
plot(qex(growth.rate), growth.rate, main = "Exponential Q-Q plot: Growth Rate",
     ylab = "Growth Residuals", xlab = "Exponential Quantiles")
qqnorm(log(growth.rate), main = "Normal Q-Q plot: log(Growth Rate)")
plot(qex(growth.rate), log(growth.rate), main = "Exponential Q-Q plot: log(Growth Rate)",
     ylab = "Log(Growth Residuals)", xlab = "Exponential Quantiles")
dev.off()
ks.test(growth.rate, "plnorm", alternative = c("two.sided"), exact = FALSE)
shapiro.test(sample(log(growth.rate), 5000))

load("/Users/rutuja/Work_at_LANL/Projects/bci.roots.inversion/results/splevel/GLUE.resid_drop.monthsNone_cor0.3_2019-10-14_5000_0_med_size_class_predefined_cc_scaled_on_5.Rdata")
resid.growth.all <- as.numeric(unlist(GLUE.resid))
resid.growth <- sample(resid.growth.all, 10000)
jpeg(file.path(figures.folder, "residual_growth_rates_distribution.jpeg", width = 960, height = 960, units = "px", pointsize = 24,
     quality = 100))
par(mfrow = c(2, 2))
qqnorm(resid.growth, main = "Normal Q-Q plot: Residuals")
plot(qex(resid.growth), resid.growth)
qqnorm(log(resid.growth), main = "Normal Q-Q plot: log(Residuals)")
plot(qex(resid.growth), log(resid.growth))
dev.off()

shapiro.test(sample(resid.growth.all, 5000))
ks.test(sample(resid.growth.all, 5000), "pnorm",alternative = c("two.sided"), exact = FALSE)
### load udi as well
# rm(ds.bestfit);
file.extension.base3 <- "drop.monthsFebMar_cor0.3_2019-10-14_5000_0_med_size_class_predefined_cc_scaled_on_5_dryseason_on_iso.subset_off"
load(file = paste0("results/splevel/best10.type.udi_", file.extension.base3, ".Rdata"), envir = parent.frame(), verbose = FALSE)

udi <- best10.type.udi[, sample(ncol(best10.type.udi), 10000)]

pdf("figures/udi_distribution.pdf")
for (i in 1:10) {
par(mfrow = c(2, 2))
ind <- apply(udi, 1, function(x) all(is.na(x)))
udi.sp <- as.numeric(udi[!ind,][i, ])
qqnorm(udi.sp, main = "Normal Q-Q plot: Uptake Depth Index")
plot(qex(udi.sp), udi.sp)
qqnorm(log(udi.sp), main = "Normal Q-Q plot: log(Uptake Depth Index)")
plot(qex(udi.sp), log(udi.sp))
}
dev.off()

for (i in 1:10) {
  ind <- apply(udi, 1, function(x) all(is.na(x)))
  udi.sp <- as.numeric(udi[!ind,][i, ])
  shapiro.test(sample(udi.sp, 5000))
  ks.test(sample(udi.sp, 5000), "plnorm", alternative = c("two.sided"), exact = FALSE)
}

#*****************
### PSI ----------
#*****************
psi.mean <- bci.elm.fates.hydro::gpp
hist(psi.mean$psi)
shapiro.test(sample(psi.mean$psi, 5000))
ks.test(sample(psi.mean$psi, 5000), "plnorm", alternative = c("two.sided"), exact = FALSE)

swp <- -psi.mean$psi

# https://daviddalpiaz.github.io/stat3202-sp19/notes/fitting.html
jpeg(file.path(figures.folder, "psi_distribution.jpeg"), width = 900, height = 600, units = "px", pointsize = 24,
     quality = 80)
par(mfrow = c(2, 3))
## Exp
hist(swp, probability = TRUE,
     main = "SWP", xlab = "SWP (-MPa)")
box()
grid()
curve(dexp(x, rate = 1 / mean(swp)),
      add = TRUE, col = "darkorange")


##gamma# calculating sample moments
len_samp_moment_1 = mean(swp)
len_samp_moment_2 = mean(swp ^ 2)

# method of moments estimators
len_alpha_mom = len_samp_moment_1 ^ 2 / (len_samp_moment_2 - len_samp_moment_1 ^ 2)
len_beta_mom  = len_samp_moment_1 / len_alpha_mom

# estimates for this dataset
c(len_alpha_mom, len_beta_mom)

hist(swp, probability = TRUE,
     main = "SWP", xlab = "SWP (-MPa)")
box()
grid()
curve(dgamma(x, shape = len_alpha_mom, scale = len_beta_mom),
      add = TRUE, col = "darkorange", lwd = 2)
## kernel density
hist(swp, probability = TRUE,
     main = "SWP", xlab = "SWP (-MPa)")
box()
grid()
lines(density(swp), col = "darkorange")

## exp qq plot
qqplot(x = qexp(ppoints(swp), rate = 1 / mean(swp)),
       y = swp,
       # xlim = c(0, 400), ylim = c(0, 400),
       main = "QQ-Plot: SWP, Exp Distribution",
       xlab = "Theoretical Quantiles, Exp Distribution",
       ylab = "Sample Quantiles, SWP")
abline(a = 0, b = 1, col = "dodgerblue", lwd = 2)
# gamma qq
qqplot(x = qgamma(ppoints(swp),
                  shape = len_alpha_mom, scale = len_beta_mom),
       y = swp,
       # xlim = c(0, 400), ylim = c(0, 400),
       main = "QQ-Plot: SWP, Gamma Distribution",
       xlab = "Theoretical Quantiles, Gamma Distribution",
       ylab = "Sample Quantiles, SWP")
abline(a = 0, b = 1, col = "dodgerblue", lwd = 2)
grid()
## kernel qq
qqplot(x = quantile.density(density(swp), ppoints(swp)),
       y = swp,
       # xlim = c(0, 400), ylim = c(0, 400),
       main = "QQ-Plot: SWP, KDE",
       xlab = "Theoretical Quantiles, Kernel Density Estimate",
       ylab = "Sample Quantiles, SWP")
abline(a = 0, b = 1, col = "dodgerblue", lwd = 2)
grid()
dev.off()
## So none of the models fit, using non-parametric Kernel Density Estimator

#*****************
### GPP ----------
#*****************

gpp <- bci.elm.fates.hydro::gpp %>% rename(gpp = value)
hist(gpp$gpp)
shapiro.test(sample(gpp$gpp, 5000))
ks.test(sample(gpp$gpp, 5000), "plnorm", alternative = c("two.sided"), exact = FALSE)

gpp <- gpp$gpp

# https://daviddalpiaz.github.io/stat3202-sp19/notes/fitting.html
jpeg(file.path(figures.folder, "gpp_distribution.jpeg"), width = 600, height = 600, units = "px", pointsize = 24,
     quality = 80)
par(mfrow = c(2, 2))
## Normal
hist(gpp, probability = TRUE,
     main = "gpp", xlab = expression('GPP (gC'*m^-2*day^-1*')'))
box()
grid()
curve(dnorm(x, mean = mean(gpp), sd = sd(gpp)),
      add = TRUE, col = "darkorange")
# ## Normal but after log of data
# hist(gpp, probability = TRUE,
#      main = "gpp", xlab = expression('GPP (gC'*m^-2*day^-1*')'))
# box()
# grid()
# curve(dlnorm(x, mean = mean(gpp), sd = sd(gpp)),
#       add = TRUE, col = "darkorange")

## kernel density
hist(gpp, probability = TRUE,
     main = "gpp", xlab = expression('GPP (gC'*m^-2*day^-1*')'))
box()
grid()
lines(density(gpp), col = "darkorange")

## exp qq plot
qqplot(x = qnorm(ppoints(gpp)),
       y = gpp,
       # xlim = c(0, 400), ylim = c(0, 400),
       main = "QQ-Plot: gpp, Normal Distribution",
       xlab = "Theoretical Quantiles, Normal Distribution",
       ylab = "Sample Quantiles, gpp")
abline(a = 0, b = 1, col = "dodgerblue", lwd = 2)
## kernel qq
qqplot(x = quantile.density(density(gpp), ppoints(gpp)),
       y = gpp,
       # xlim = c(0, 400), ylim = c(0, 400),
       main = "QQ-Plot: gpp, KDE",
       xlab = "Theoretical Quantiles, Kernel Density Estimate",
       ylab = "Sample Quantiles, gpp")
abline(a = 0, b = 1, col = "dodgerblue", lwd = 2)
grid()
dev.off()

#*****************
# Growth -------
#*****************

intervals <- 5
growth.selection <- "size_class_predefined_cc_scaled"
load(file = paste0("results/gro.long.cc_", intervals, "_", growth.selection, ".Rdata"))
# growth.var <- "dbh.residuals"
growth.var <- "growth"
if(growth.var == "growth") {
  growth <- gro.long.cc$growth
} else {
  growth <- gro.long.cc$dbh.residuals
}

# https://daviddalpiaz.github.io/stat3202-sp19/notes/fitting.html
jpeg(file.path(figures.folder, paste0(growth.var, "_distribution.jpeg")), width = 600, height = 600, units = "px", pointsize = 24,
     quality = 80)
par(mfrow = c(2, 3))
## Normal
hist(growth, probability = TRUE,
     main = growth.var, xlab = expression('growth (gC'*m^-2*day^-1*')'))
box()
grid()
curve(dnorm(x, mean = mean(growth), sd = sd(growth)),
      add = TRUE, col = "darkorange")
## Log-Normal
hist(growth, probability = TRUE,
     main = growth.var, xlab = expression('growth (gC'*m^-2*day^-1*')'))
box()
grid()
lnorm_fit_params <- fitdistr(growth[growth > 0], "lognormal")
curve(dnorm(x, lnorm_fit_params$estimate['meanlog'], lnorm_fit_params$estimate['sdlog']),
      add = TRUE, col = "darkorange")
## kernel density
hist(growth, probability = TRUE,
     main = growth.var, xlab = expression('growth (gC'*m^-2*day^-1*')'))
box()
grid()
lines(density(growth), col = "darkorange")

## Normal qq plot
qqplot(x = qnorm(ppoints(growth)),
       y = growth,
       # xlim = c(0, 400), ylim = c(0, 400),
       main = "Normal ",
       xlab = "Theoretical Quantiles, Normal Distribution",
       ylab = "Sample Quantiles, growth")
abline(a = 0, b = 1, col = "dodgerblue", lwd = 2)
# Log Normal qq plot
qqplot(x = qlnorm(ppoints(growth)),
       y = growth,
       # xlim = c(0, 400), ylim = c(0, 400),
       main = "Log Normal",
       xlab = "Theoretical Quantiles, Log Normal Distribution",
       ylab = "Sample Quantiles, growth")
abline(a = 0, b = 1, col = "dodgerblue", lwd = 2)

## kernel qq
qqplot(x = quantile.density(density(growth), ppoints(growth)),
       y = growth,
       # xlim = c(0, 400), ylim = c(0, 400),
       main = "KDE",
       xlab = "Theoretical Quantiles, Kernel Density Estimate",
       ylab = "Sample Quantiles, growth")
abline(a = 0, b = 1, col = "dodgerblue", lwd = 2)
grid()
dev.off()

#*****************
## Leaf Fall -----
#*****************
leaf.fall <- read.csv(file.path("data-raw/traits/Wright_Osvaldo_BCI_weekly_leaf-fall_data/Rutuja.csv"))

leaf.fall <- leaf.fall$leaf_gm
leaf.fall.var <- "Leaf Fall"
jpeg(file.path(figures.folder, paste0("leaf.fall_distribution.jpeg")), width = 600, height = 600, units = "px", pointsize = 24,
     quality = 80)
par(mfrow = c(2, 4))
## Normal
hist(leaf.fall, probability = TRUE,
     main = leaf.fall.var, xlab = expression('leaf.fall (gC'*m^-2*day^-1*')'))
box()
grid()
curve(dnorm(x, mean = mean(leaf.fall), sd = sd(leaf.fall)),
      add = TRUE, col = "darkorange")
# ## Log-Normal
# hist(leaf.fall, probability = TRUE,
#      main = leaf.fall.var, xlab = expression('leaf.fall (gC'*m^-2*day^-1*')'))
# box()
# grid()
# lnorm_fit_params <- fitdistr(leaf.fall[leaf.fall> 0], "lognormal")
# curve(dnorm(x, lnorm_fit_params$estimate['meanlog'], lnorm_fit_params$estimate['sdlog']),
#       add = TRUE, col = "darkorange")
#Exp
hist(leaf.fall, probability = TRUE,
     main = leaf.fall.var, xlab = expression('leaf.fall (gC'*m^-2*day^-1*')'))
box()
grid()
curve(dexp(x, rate = 1 / mean(leaf.fall)),
      add = TRUE, col = "darkorange")
# Gamma
##gamma# calculating sample moments
len_samp_moment_1 = mean(leaf.fall)
len_samp_moment_2 = mean(leaf.fall ^ 2)

# method of moments estimators
len_alpha_mom = len_samp_moment_1 ^ 2 / (len_samp_moment_2 - len_samp_moment_1 ^ 2)
len_beta_mom  = len_samp_moment_1 / len_alpha_mom

hist(leaf.fall, probability = TRUE,
     main = leaf.fall.var, xlab = expression('leaf.fall (gC'*m^-2*day^-1*')'))
box()
grid()
curve(dgamma(x, shape = len_alpha_mom, scale = len_beta_mom),
      add = TRUE, col = "darkorange", lwd = 2)
## kernel density
hist(leaf.fall, probability = TRUE,
     main = leaf.fall.var, xlab = expression('leaf.fall (gC'*m^-2*day^-1*')'))
box()
grid()
lines(density(leaf.fall), col = "darkorange")

## Normal qq plot
qqplot(x = qnorm(ppoints(leaf.fall)),
       y = leaf.fall,
       main = "Normal ",
       xlab = "Theoretical Quantiles, Normal Distribution",
       ylab = "Sample Quantiles, leaf.fall")
abline(a = 0, b = 1, col = "dodgerblue", lwd = 2)
# # Log Normal qq plot
# qqplot(x = qlnorm(ppoints(leaf.fall)),
#        y = leaf.fall,
#        # xlim = c(0, 400), ylim = c(0, 400),
#        main = "Log Normal",
#        xlab = "Theoretical Quantiles, Log Normal Distribution",
#        ylab = "Sample Quantiles, leaf.fall")
# abline(a = 0, b = 1, col = "dodgerblue", lwd = 2)
# Exp
qqplot(x = qexp(ppoints(leaf.fall), rate = 1/ mean(leaf.fall, na.rm = TRUE)),
       y = leaf.fall,
       main = "Exp",
       xlab = "Theoretical Quantiles, Exp Distribution",
       ylab = "Sample Quantiles, leaf.fall")
abline(a = 0, b = 1, col = "dodgerblue", lwd = 2)
# Gamma
# gamma qq
qqplot(x = qgamma(ppoints(leaf.fall),
                  shape = len_alpha_mom, scale = len_beta_mom),
       y = leaf.fall,
       main = "Gamma",
       xlab = "Theoretical Quantiles, Gamma Distribution",
       ylab = "Sample Quantiles, leaf.fall")
abline(a = 0, b = 1, col = "dodgerblue", lwd = 2)
grid()
## kernel qq
qqplot(x = quantile.density(density(leaf.fall), ppoints(leaf.fall)),
       y = leaf.fall,
       main = "KDE",
       xlab = "Theoretical Quantiles, Kernel Density Estimate",
       ylab = "Sample Quantiles, leaf.fall")
abline(a = 0, b = 1, col = "dodgerblue", lwd = 2)
grid()
dev.off()




