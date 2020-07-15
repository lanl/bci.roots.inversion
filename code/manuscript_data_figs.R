#******************************************************
## output
#******************************************************
load(file = file.path(results.folder, "psi.stat.4.select.Rdata"))
load(file = file.path(results.folder, "ml.rsq.combine.best.Rdata"))
load(file = file.path(results.folder, "mrate.depth.Rdata"))
load(file = file.path(results.folder, "mrate.mfac.depth.Rdata"))

#****************************
###   Custom Functions   ####
#****************************
rectangles.4 <- data.frame(
  xmin = 120,
  xmax = 335,
  ymin = 0,
  ymax = -3.0
)
rev_sqrt_trans <- function() {
  scales::trans_new(
    name = "rev_sqrt",
    transform = function(x) -sqrt(abs(x)),
    inverse = function(x) x^2);
}
formula <- y ~ x
#****************************

###
ml.rsq.combine.best <- ml.rsq.combine.best %>% left_join(bci.traits %>% dplyr::select(sp, form1), by = "sp") %>%
  mutate(depth = as.numeric(depth))
erd.iso <- ml.rsq.combine.best %>%
  subset(corr >= 0 & form1 == "T" & sp != "guapst" &
           corr.func == "gr.Psi.VPD" & source == "Meinzer et al.1999 Fig. 4" &
           R2 >= 0.1)
## Minimum Soil water potential reached at depth 1.7 + CI

psi.1.7.min <- subset(psi.stat.4.select, depth == 1.7) %>% subset(median == min(median, na.rm = TRUE))


## Graph SWP for depths used in inverse model
plot.psi.stat.6.interval.base.select <- ggplot(psi.stat.4.select %>% subset(depth %in% c(0.5, 1, 1.7))
                                               %>% droplevels()) +
  scale_x_continuous(breaks = c(seq(0, 360, by = 60))) +
  # coord_cartesian(ylim = c(-2.5, 0)) +
  theme(panel.grid.major.y = element_line()) +
  ylab(expression(Psi[soil]*~"(MPa)")) + xlab("Day of the Year") +
  facet_wrap(. ~ interval.yrs.2) +
  geom_rect(data=rectangles.4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill='gray90', alpha = 0.8) +
  geom_line(aes(x = doy, y = median.clim, group = as.factor(depth), color = as.factor(depth)), size = 0.3, linetype = "solid") +
  geom_ribbon(aes(x = doy, ymin = below.q5, ymax = median.clim, group = as.factor(depth_year),
                  fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE) +
  theme(panel.grid.major.y = element_line(size = 0.1)) +
  scale_color_discrete(name = "Depth (m)", labels = c("0.5", "1", "1.7")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  coord_cartesian(ylim = c(-2, 0), xlim = c(0, 200)) +
  theme(legend.position = "top")

plot.psi.stat.7.interval.q2.5.select <- plot.psi.stat.6.interval.base.select %+%
  subset(psi.stat.4.select, depth %in% c(0.5, 1, 1.7) & interval.yrs != "(Missing)") +
  scale_y_continuous(trans = rev_sqrt_trans(), breaks = c(0.01, 0.1, 0.5, 1, 1.5, 2.0),
                     labels = c("0.01", "-0.1", "-0.5", "-1.0", "-1.5", "-2.0")) +
  facet_wrap(. ~ interval.yrs, nrow = 1) +
  geom_ribbon(aes(x = doy, ymin = below.q2.5, ymax = median.clim, group = as.factor(depth_year),
                  fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE)
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_study_period_q2.5_depths_in_inverse_model.jpeg",
       plot = plot.psi.stat.7.interval.q2.5.select, file.path(figures.folder), device = "jpeg", height = 2.5, width = 6, units='in')

mrate.depth.select <- subset(mrate.depth, !is.na(rdi.gr) & avg.abund >= 20)
mrate.mfac.depth.select <- subset(mrate.mfac.depth, !is.na(rdi.gr) & avg.abund >= 20)

mrate.depth.mean <- mrate.depth.select %>%
  group_by(sp, rdi.gr) %>% summarise(se = sd(mrate, na.rm = TRUE)/sqrt(n()),
                                     mrate = mean(mrate, na.rm = TRUE)) %>% droplevels()
## lowest Psi_crit_
# max(subset(data.model.AB, sp %in% unique(mort.erd.to.plot$sp)$psi_kl80, na.rm = TRUE)
mfac.plot.15 <- ggplot(mrate.depth.mean,
                       aes(y = mrate, x = rdi.gr)) +
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymin = mrate - se, ymax = mrate + se), width = 0.2, size = 0.1) +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 1, size = 2.5) +
  ylab(expression('Mean Mortality Rate (%'*'year'^1*')')) +
  xlab("Effective Rooting Depth (m)") +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.95, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 6) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.95, npcy = 0.82, size = 6) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by rdi.gr.jpeg")),
       plot = mfac.plot.15, height = 3, width = 3, units='in')

y.label.1 <- expression(Mort[anomaly]~('%'*yr^{-1}))
mfac.plot.15.1 <- ggplot(mrate.depth.select, aes(y = mrate, x = rdi.gr)) +
  coord_cartesian(xlim = c(0, max(mrate.depth$rdi.gr, na.rm = TRUE))) +
  geom_smooth(method = "lm") +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 0.8, size = 2.5) +
  xlab("Effective Rooting Depth (m)")  + ylab(y.label.1) +
  facet_grid(. ~ censusint.m) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.95, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.95, npcy = 0.82, size = 4)
ggsave(file.path(paste0(figures.folder, "/diff.mortality_by rdi.gr.jpeg")),
       plot = mfac.plot.15.1, height = 2.5, width = 10, units = 'in')

mrate.mfac.depth.gr.mean.mfac <- mrate.mfac.depth.select %>%
  subset(depth == rdi.gr) %>%
  group_by(sp) %>%
  summarise(depth = mean(depth, na.rm = TRUE),
            mfac = sum(mfac, na.rm = TRUE))

mfac.plot.9.0 <- ggplot(mrate.mfac.depth.gr.mean.mfac, aes(x = mfac, y = depth)) +
  scale_y_continuous(trans = "rev_sqrt", breaks =
                       c(0, sort(unique(mrate.mfac.depth.gr.mean.mfac$depth)))) +
  geom_jitter(height = 0.1, width = 0, size = 2, shape = 21, alpha = 0.6, color = "black", aes(fill = sp), show.legend = FALSE) +
  ylab("Effective Rooting Depth (m)") +
  xlab(expression(atop('Time spent below '*Psi['crit'], '(Days over 1990-2015)')))
ggsave(file.path(paste0(figures.folder, "/mean_mfac vs. rdi.gr.jpeg")),
       plot = mfac.plot.9.0, height = 3.5, width = 3.5, units='in')
