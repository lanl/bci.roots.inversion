#******************************************************
## output
#******************************************************
load(file = file.path(results.folder, "psi.stat.4.Rdata"))
load(file = file.path(results.folder, "psi.stat.4.select.Rdata"))
load(file = file.path(results.folder, "psi.stat.4.select.freq.Rdata"))
load(file = file.path(results.folder, "ml.rsq.combine.Rdata"))
load(file = file.path(results.folder, "ml.rsq.combine.best.Rdata"))
load(file = file.path(results.folder, "chosen.model.Rdata"))
load(file = file.path(results.folder, "erd.model.p.Rdata"))
load(file = file.path(results.folder, "mrate.depth.Rdata"))
load(file = file.path(results.folder, "mrate.mfac.depth.Rdata"))
load(file = file.path(results.folder, "erd.stem.traits.Rdata"))
load(file = file.path(results.folder, "df.erd.to.plot.Rdata"))
load(file = file.path(results.folder, "obs.data.model.AB.Rdata"))
load(file = file.path(results.folder, "data.model.AB.sub.Rdata"))
load(file = file.path(results.folder, "obs.sp.vcurves.1.Rdata"))
load(file = file.path(results.folder, "comm.sp.vcurves.1.Rdata"))
load(file = file.path(results.folder, "bci.lifetime.Rdata"))


#****************************
###   Custom Functions   ####
#****************************
rev_sqrt_trans <- function() {
  scales::trans_new(
    name = "rev_sqrt",
    transform = function(x) -sqrt(abs(x)),
    inverse = function(x) x^2);
}
reverselog_trans <- function(base = exp(1)) {
  scales::trans_new(name = paste0("reverselog-", format(base)),
                    log_breaks(base = base),
                    domain = c(1e-100, Inf),
                    transform = function(x) -log(x, base),
                    inverse = function(x) base^(-x))
}

formula <- y ~ x
#****************************
ml.rsq.combine.best <- ml.rsq.combine.best %>%
  ## ERD only for canopy species, so no need to subset
  mutate(depth = as.numeric(depth))
erd.data <- ml.rsq.combine.best %>%
  subset(corr.func == chosen.model) %>%
  distinct(sp, .keep_all = TRUE) # removes the dplicate jac1co due to Jackson dataset
erd.iso <- erd.data %>%
  subset(sp != "guapst" & source == "Meinzer et al.1999 Fig. 4")

#******************************************************
### ERD Species names----
#******************************************************
erd.sp <- erd.data$sp
save(erd.sp, file = file.path("results", "erd.sp.Rdata"))

erd.sp.names <- bci.traits %>%
  subset(sp %in% erd.sp) %>%
  dplyr::rename(Code = sp, Genus = GENUS., Species = SPECIES., Family = FAMILY.) %>%
  dplyr::select(Code, Genus, Species, Family) %>%
  # correct sp names
  mutate(Genus = ifelse(Genus == "Tabebuia", "Handroanthus", as.character(Genus)),
         Genus = ifelse(Genus == "Beilschmiedi", "Beilschmiedia", as.character(Genus)),
         Genus = ifelse(Genus == "Trattinnicki", "Trattinnickia", as.character(Genus)),
         Species = ifelse(Species == "costaricensi", "costaricensis", as.character(Species)))

rownames(erd.sp.names) <- 1: nrow(erd.sp.names)

#****************************
### ELM-FATES param sensitivity----
#****************************


sm1 <- cowplot::ggdraw() + cowplot::draw_image("figures/ELM_FATES_parameter_sensitivity/Sensitivity of daily soil water content at 1 m to parameters_2016-04.jpeg", scale = 1)
sm2 <- cowplot::ggdraw() + cowplot::draw_image("figures/ELM_FATES_parameter_sensitivity/Sensitivity of daily soil water content at 1 m to parameters_2016-07.jpeg", scale = 1)

et1 <- cowplot::ggdraw() + cowplot::draw_image("figures/ELM_FATES_parameter_sensitivity/Sensitivity of monthly ET to parameters_2015-02.jpeg", scale = 1)
et2 <- cowplot::ggdraw() + cowplot::draw_image("figures/ELM_FATES_parameter_sensitivity/Sensitivity of monthly ET to parameters_2016-04.jpeg", scale = 1)

sens.sm.et <- cowplot::plot_grid(sm1, sm2, et1, et2, labels = c('a', 'b', 'c', 'd'),
                              label_size = 14, nrow = 2, ncol = 2, rel_widths = c(1, 1))
ggsave("sens.sm.et.jpeg", plot = sens.sm.et, path =
         file.path("figures/ELM_FATES_parameter_sensitivity"), device = "jpeg", height = 7, width = 7.5, units ='in')

run1 <- cowplot::ggdraw() + cowplot::draw_image("figures/ELM_FATES_parameter_sensitivity/Sensitivity of monthly runoff to parameters_2015-02.jpeg", scale = 1)
run2 <- cowplot::ggdraw() + cowplot::draw_image("figures/ELM_FATES_parameter_sensitivity/Sensitivity of monthly runoff to parameters_2015-07.jpeg", scale = 1)
run3 <- cowplot::ggdraw() + cowplot::draw_image("figures/ELM_FATES_parameter_sensitivity/Sensitivity of monthly runoff to parameters_2015-10.jpeg", scale = 1)
run4 <- cowplot::ggdraw() + cowplot::draw_image("figures/ELM_FATES_parameter_sensitivity/Sensitivity of monthly runoff to parameters_2016-07.jpeg", scale = 1)

sens.run <- cowplot::plot_grid(run1, run2, run3, run4, labels = c('a', 'b', 'c', 'd'),
                                   label_size = 14, ncol = 2, nrow = 2, rel_widths = c(1, 1))
ggsave("sens.run.jpeg", plot = sens.run, path =
         file.path("figures/ELM_FATES_parameter_sensitivity"), device = "jpeg", height = 7, width = 7.5, units ='in')

#****************************
### PSI significant droughts----
#****************************

psi.rectangle.1 <- data.frame(
  xmin = 120,
  xmax = 310,
  ymin = 1,
  ymax = -2.5,
  type = "Wet Season"
)
### Calculate Correlation of growth rates with psi by depth

plot.psi.stat.7.interval.q5.depth.base <-
  ggplot(subset(psi.stat.4, depth %in% c(0.12, 0.62, 1) &
                  interval.yrs != "(Missing)") %>% droplevels()) +
  geom_rect(data = psi.rectangle.1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill = type), color = NA,
            alpha = 0.8, size = 0.3) + #  #cce6ff
  geom_ribbon(aes(x = doy, ymin = q5.clim, ymax = q100.clim, group = as.factor(depth),
                  fill = "95% CI"), alpha = 0.5) +
  theme(panel.grid.major.y = element_line(size = 0.1)) +
  # geom_line(data = psi.stat.4 %>%
  #             subset(extreme.yr.q2.5 & depth %in% c(0.12, 0.62, 1) & interval.yrs != "(Missing)"),
  #           aes(x = doy, y = median, group = as.factor(depth_year),
  #               color = year), size = 0.5, alpha = 1) +
  geom_line(aes(x = doy, y = median.clim, group = as.factor(depth),
                linetype = "Median"), color = "gray30", size = 0.5) +
  geom_line(data = psi.stat.4 %>%
              subset(extreme.yr.q2.5 & depth %in% c(0.12, 0.62, 1) & interval.yrs != "(Missing)"),
            aes(x = doy, y = below.q5, group = as.factor(depth_year),
                color = year), size = 1.3, alpha = 1) +
  scale_linetype_manual(name = "", values = c("solid")) +
  scale_fill_manual(name = "", values = c("gray80", "#e6f2ff")) +
  guides(linetype = guide_legend(order = 1, title = NULL, label.position = "top"),
         fill = guide_legend(order = 2, title = NULL, label.position = "top"),
         color = guide_legend(order = 3, title = "Year",
                              override.aes = list(size = 2))) +
  # scale_x_continuous(breaks = c(seq(0, 360, by = 60))) +
  coord_cartesian(ylim = c(-2, 0), xlim = c(0, 200)) +
  ylab(expression(Psi[soil]*~"(MPa)")) + xlab("Day of the Year")

plot.psi.stat.7.interval.q5.depth <- plot.psi.stat.7.interval.q5.depth.base +
  facet_grid(plot.depth ~ interval.yrs) +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_study_period_q5_by_depth.jpeg",
       plot = plot.psi.stat.7.interval.q5.depth, file.path(figures.folder), device = "jpeg", height = 5, width = 8, units='in')

plot.psi.stat.7.interval.q5.depth.freq.base <-
  ggplot(subset(psi.stat.4.select.freq, depth %in% c(0.5, 1, 1.7) &
                  interval.yrs != "(Missing)") %>% droplevels()) +
  facet_wrap(depth ~ interval.yrs, nrow = 3) +
  geom_count(aes(x = doy, y = freq.below.q5, color = depth), size = 0.5) +
  theme(panel.grid.major.y = element_line(size = 0.1)) +
  scale_linetype_manual(name = "", values = c("solid")) +
  guides(linetype = guide_legend(order = 1, title = NULL, label.position = "top"),
         fill = guide_legend(order = 2, title = NULL, label.position = "top"),
         color = guide_legend(order = 3, title = "Depth (m)",
                              override.aes = list(size = 3))) +
  scale_x_continuous(breaks = c(seq(0, 360, by = 60))) +
  scale_y_continuous(breaks = c(0:5)) +
  coord_cartesian(xlim = c(0, 200)) +
  ylab(expression('Frequency of extreme '*Psi[soil]*~"(yrs)")) + xlab("Day of the Year")
## Minimum Soil water potential reached at depth 1.7 + CI
psi.1.7.min <- subset(psi.stat.4.select, depth == 1.7) %>%
  subset(median == min(median, na.rm = TRUE))

# Heatmap
psi.freq.to.plot <- psi.stat.4.select.freq %>%
  subset(depth %in% c(0.5, 1, 1.7) &
           interval.yrs != "(Missing)") %>%
  mutate(depth.plot = factor(depth, levels = c(1.7, 1, 0.5))) %>%
  droplevels()
psi.freq.rectangle <- data.frame(xmin= 120, xmax=310, ymin=0.5, ymax=3.5, type = "Wet Season")
droughts.psi.heat <- ggplot(psi.freq.to.plot) +
  # scale_y_reverse() +
  facet_grid(. ~ interval.yrs) +
  geom_tile(aes(x = doy, y = depth.plot,
                fill = as.factor(freq.below.q5))) +
  geom_rect(data = psi.freq.rectangle, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, color = type),
            alpha = 0.8, fill = NA, size = 0.8) +
  ylab("Depth (m)") + xlab("Day of Year") +
  scale_color_manual(values = "dodgerblue") +
  scale_fill_manual(values = c("gray90", "#feb24c", "#f03b20"), breaks = c(0, 1, 2)) +
  # scale_fill_manual(values =   c("gray80", "#fec44f", "#d95f0e"), breaks = c(0, 1, 2)) +

  guides(fill = guide_legend(order = 1, title = expression("Frequency of extreme soil drought (years)"),
                             direction = "horizontal",
                              override.aes = list(size = 1)),
         color = guide_legend(order = 2, title = NULL,
                             direction = "horizontal")) +
  theme(legend.position = "top")
ggsave("Frequency of extreme soil droughts_heatmap.jpeg", plot = droughts.psi.heat, path =
         file.path(figures.folder), device = "jpeg", height = 2.5, width = 9, units='in')


#****************************
## LAI seasonality by species----
#****************************

f4 <- ggplot(sp.leaf_cover.for.model %>% subset(sp %in% erd.sp & doy != 366) %>%
               left_join(erd.sp.names %>%
                           mutate(s.names = paste0(substr(Genus, start = 1, stop = 1), ". ", tolower(Species)),
                                  sp = tolower(Code)), by = "sp"),
             aes(x = doy, y = leaf_cover)) +
  facet_wrap(s.names ~ ., scales = "free_y", ncol = 4) +
  ylim(c(0, 1)) +
  geom_line(aes(group = sp, color = deciduousness), size = 0.8) +
  ylab("Leaf Cover Fraction") + xlab("DOY") +
  guides(color = guide_legend(order = 1, title = NULL, direction = "horizontal",
                              override.aes = list(size = 3))) +
  theme(legend.position = "top", legend.title = element_blank()) +
  scale_color_viridis_d(drop = FALSE) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 1, hjust = 1),
        strip.text.x = element_text(face = "italic", size = 11))
ggsave(("leaf.cover_BCI_multi_panel.jpeg"),
       plot = f4, file.path(figures.folder), device = "jpeg", height = 9, width = 7, units='in')

#****************************
## ERD by species----
#****************************

df.erd.to.plot <- df.erd.to.plot %>%
  left_join(erd.sp.names %>%
              mutate(s.names = paste0(substr(Genus, start = 1, stop = 1), ". ", tolower(Species)),
                     sp = tolower(Code)), by = "sp") %>%
  dplyr::select(sp, s.names, depth, depth.se) %>%
  transform(s.names = reorder(s.names, depth)) %>%
  droplevels()

erd.sp.plot <- ggplot(df.erd.to.plot,
                      aes(x = s.names, y = depth)) +
  # geom_point(aes(color = s.names), show.legend = FALSE, size = 3) +
  # geom_point(shape = 21, color = "white", fill = "black", alpha = 1, size = 3) +
  geom_point(aes(color = depth), alpha = 1, size = 3, show.legend = FALSE) +
  scale_color_continuous(trans = 'reverse') +
  geom_errorbar(aes(ymax = depth + depth.se, ymin = depth - depth.se), width = 0.2, size = 0.2) +
  ylab("Effective Rooting Depth (m)") + xlab("Species") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
        axis.text.y = element_text(face = "plain")) +
  scale_y_continuous(trans = "reverse", breaks = unique(ml.rsq.combine$depth))
ggsave("ERD_by_sp_large_canopy.tiff",
       plot = erd.sp.plot, file.path(figures.folder), device = "tiff", height = 3.5, width = 5, units='in')
ggsave("ERD_by_sp_large_canopy.jpeg",
       plot = erd.sp.plot, file.path(figures.folder), device = "jpeg", height = 4.5, width = 5, units='in')

#****************************
## ERD by isotopes----
#****************************

ml.rsq.combine.sub <- ml.rsq.combine.best %>%
  mutate(depth = as.numeric(depth)) %>%
  subset(!sp %in% c("guapst") & !is.na(Xylem_sap_deltaD_permil.mean)) %>%
  left_join(erd.sp.names %>%
              dplyr::rename(sp = Code) %>%
              mutate(s.names = paste0(substr(Genus, start = 1, stop = 1), ". ", Species)) %>%
              dplyr::select(sp, s.names), by = "sp") %>%
  subset(source == "Meinzer et al.1999 Fig. 4") %>%
  droplevels()

formula = y~x
ml.rsq.combine.chosen <- ml.rsq.combine.sub %>% subset(corr.func == chosen.model)
lm.erd.iso <- lm(Xylem_sap_deltaD_permil ~ depth, data = ml.rsq.combine.chosen)
lm.erd.iso.r2 <- round(summary(lm.erd.iso)$r.squared, 2); lm.erd.iso.p <- round(coef(summary(lm.erd.iso))[2,4], 2)

xylem.label <- expression(delta^2*H[xylem]~"( \u2030)")
p4 <- ggplot(ml.rsq.combine.chosen,
             aes(x = Xylem_sap_deltaD_permil, y = depth)) +
  geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + se,
                     xmin = Xylem_sap_deltaD_permil - se, color = s.names),
                 size = 0.5, height = 0.05, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, size = 0.5, formula = formula, color = "gray10") +
  ylab(expression("Effective Rooting Depth (m)")) + xlab(xylem.label) +
  scale_y_continuous(trans="reverse", breaks = unique(ml.rsq.combine$depth)) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.9, npcy = 0.2, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = ifelse(p.value < 0.001, sprintf('italic(p)~"< 0.001"'),
                                     sprintf('italic(p)~"="~%.2f',stat(p.value)))),
                  parse = TRUE,
                  npcx = 0.9, npcy = 0.1, size = 4) +
  geom_point(shape = 21, color = "white", aes(fill = s.names), alpha = 1, size = 3.5) +
  geom_errorbar(aes(ymax = depth + depth.se, ymin = depth - depth.se), color = "black",
                size = 0.3, width = 0.2) +
  guides(fill = guide_legend(title = "Species"), color = FALSE) +
  theme(legend.text = element_text(face = "italic", size = 8))
ggsave("psi.corr_best.depth_xylem_sap_deltaD_phenology_Meinzer_gr.Psi.VPD.jpeg",
       plot = p4, file.path(figures.folder), device = "jpeg", height = 3, width = 4.3, units = 'in')


## as "gr.Psi.VPD.add" is not present
df.pval <- bind_rows(erd.model.p) %>% pivot_longer(cols = everything(), names_to = "corr.func", values_to = "pval")
ml.rsq.combine.sub <- ml.rsq.combine.sub %>%
  transform(models.plot1 = factor(corr.func, levels = c("gr.Psi", "gr.Psi.VPD.add", "gr.Psi.VPD.multi",
                                                        "gr.Psi.leaf", "gr.Psi.VPD.leaf.add", "gr.Psi.VPD.leaf.multi"),
                                  labels = c("a", "b", "c", "d", "e", "f"))) %>%
  left_join(df.pval, by = "corr.func") %>%
  mutate(significant = ifelse(pval < 0.05 | pval == 0, TRUE, FALSE))
erd.iso_sp_N_by_model <- ml.rsq.combine.sub %>%
  group_by(models.plot1) %>%
  summarise(N = n(), .groups = "drop_last")

p3.2 <- ggplot(ml.rsq.combine.sub,
               aes(x = Xylem_sap_deltaD_permil, y = depth)) +
  # coord_cartesian(ylim = c(13, 0.3)) +
  geom_smooth(data = ml.rsq.combine.sub %>% subset(significant), aes(group = models.plot1),
              method = "lm", se = TRUE, size = 0.5, formula = formula, color = "gray10") +
  geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + se,
                     xmin = Xylem_sap_deltaD_permil - se, color = s.names),
                 size = 0.5, height = 0.05) +
  geom_errorbar(aes(ymax = depth + depth.se, ymin = depth - depth.se, color = s.names), size = 0.5, width = 0.2) +
  facet_wrap( ~ models.plot1, nrow = 2) +
  geom_text(data = erd.iso_sp_N_by_model, aes(x = -60, y = 0.01,
                                              label = sprintf('italic(n)~"="~%.2f', N), group = models.plot1),
            parse = TRUE, vjust = "inward", hjust = "inward", inherit.aes = FALSE) +
  ylab(expression("Effective Rooting Depth (m)")) + xlab(xylem.label) +
  scale_y_continuous(trans="reverse", breaks = unique(ml.rsq.combine.sub$depth)) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.95, npcy = 0.15, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = ifelse(p.value < 0.001, sprintf('italic(p)~"< 0.001"'),
                                     sprintf('italic(p)~"="~%.2f',stat(p.value)))),
                  parse = TRUE, npcx = 0.95, npcy = 0.05, size = 4) +
  geom_point(shape = 21, color = "white", aes(fill = s.names), alpha = 1, size = 3) +
  guides(fill = guide_legend(title = "Species", order = 1),
         color = FALSE) +
  theme(legend.text = element_text(face = "italic"),
        strip.text.x = element_text(face = "bold"))
ggsave("psi.corr_best.depth_xylem_sap_deltaD_sp_color_Meinzer.jpeg",
       plot = p3.2, file.path(figures.folder), device = "jpeg", height = 5, width = 7, units = 'in')

#******************************************************

#****************************
## Mrate mfac data prep----
#****************************

mrate.depth.select <- subset(mrate.depth, !is.na(rdi.gr) & avg.abund >= 20) %>%
  subset(sp %in% erd.sp) %>% droplevels()
mrate.mfac.depth.select <- subset(mrate.mfac.depth, !is.na(rdi.gr) &
                                    avg.abund >= 20 & depth == rdi.gr) %>%
  subset(sp %in% erd.sp) %>% droplevels() %>%
  dplyr::select(interval.num, censusint.m, sp, depth, depth.se, avg.abund, trees, mfac, days, mfac.rate, mrate,
                mean.mrate, diff.mrate, mean.grate, grate.se, size, deciduous)

# summarising across interval
mrate.depth.mean <- mrate.depth.select %>%
  group_by(sp, deciduous) %>% summarise(rdi.gr = mean(rdi.gr, na.rm = TRUE),
                                        depth.se = mean(depth.se, na.rm = TRUE),
                                        avg.abund = mean(avg.abund, na.rm = TRUE),
                                        mrate.se = sd(mrate, na.rm = TRUE)/sqrt(n()),
                                        mrate = mean(mrate, na.rm = TRUE), # should be same as mrate = mean(mean.mrate, na.rm = TRUE),
                                        grate = mean(mean.grate, na.rm = TRUE),
                                        grate.se = mean(grate.se, na.rm = TRUE),
                                        .groups = "drop_last") %>% droplevels()
erd.mrate.sp <- unique(mrate.depth.mean$sp)
mrate.mfac.depth.gr.mean.mfac <- mrate.mfac.depth.select %>%
  group_by(sp, deciduous) %>%
  summarise(depth.se = mean(depth.se, na.rm = TRUE),
            avg.abund = mean(avg.abund, na.rm = TRUE),
            depth = mean(depth, na.rm = TRUE),
            mfac = sum(mfac, na.rm = TRUE),
            days = sum(days, na.rm = TRUE),
            mfac.rate = mean(mfac.rate, na.rm = TRUE), .groups = "drop_last")
# save.image("results/manuWorkSpace.RData")

#****************************
### Hydraulic traits vs. ERD---
#****************************

hyd.table <-  erd.stem.traits %>%
  subset(!trait %in% c("lwp.min_Diurnal", "HSM88S")) %>%
  dplyr::select(sp, trait, value) %>%
  pivot_wider(names_from = trait, values_from = value) %>%
  left_join(erd.sp.names %>% dplyr::rename(sp = Code), by = c("sp")) %>%
  subset(sp %in% erd.sp) %>%
  select(-sp) %>%
  relocate(Genus, Species, Family, lwp.min_Predawn, TLP) %>%
  mutate(lwp.min_Predawn = round(lwp.min_Predawn, 2),
         TLP = round(TLP, 2))

erd.stem.traits.only <- erd.stem.traits %>%
  subset(!trait %in% c("lwp.min_Diurnal", "lwp.min_Predawn")) %>%
  left_join(df.erd.to.plot %>%
              dplyr::select(sp, depth, depth.se), by = "sp") %>%
  subset(!is.na(depth)) %>%
  droplevels()
erd.stem.traits.sp <- unique(erd.stem.traits.only$sp)

traits.labels.select <- data.frame(trait = factor(c("KmaxS", "TLP", "p88S", "HSM88S"),
                                                  levels = c("KmaxS", "TLP", "p88S", "HSM88S"), ordered = TRUE),
                                   panel = factor(c("a", "b", "c", "d"), levels = c("a", "b", "c", "d"), ordered = TRUE),
                                   x = c(0, 0, 0, 0),
                                   y = c(7.8, 0, 0.5, 2.4)) %>%
  transform(trait.plot = factor(trait, labels = c(expression(italic('K')['max, stem']~(kg*~m^-1*~s^-1*~MPa^-1)),
                                                  expression(Psi[tlp]~(MPa)),
                                                  expression(Psi['88, stem']~(MPa)),
                                                  expression(Psi[min]*' - '*Psi['88, stem']~(MPa)))),
            trait.plot.chart = factor(trait, labels = c(expression(italic('K')['max,stem']),
                                                  expression(Psi[tlp]),
                                                  expression(Psi['88,stem']),
                                                  expression(Psi[min]*'-'*Psi['88,stem']))))

erd.stem.traits.only.lab <- erd.stem.traits.only %>%
  left_join(traits.labels.select %>%
              dplyr::select(trait, trait.plot, trait.plot.chart), by = "trait") %>%
  droplevels()

#****************************
## Correlation chart: ERD vs. hydraulic traits----
#****************************
# Check correlations (as scatterplots), distribution and print correlation coefficient

erd.pairs <- erd.stem.traits.only.lab %>%
  dplyr::select(sp, depth, trait.plot.chart, value) %>%
  pivot_wider(names_from = trait.plot.chart, values_from = value) %>%
  rename(ERD = depth)

p_load(corrplot, ggcorrplot)
cor.data <- erd.pairs %>% dplyr::select(-sp) %>%
  relocate(ERD, `italic(\"K\")[\"max,stem\"]`)
M <- cor(cor.data, use = "pairwise.complete.obs", method = "spearman")

res1 <- cor.mtest(cor.data, conf.level = 0.95,
                  use = "pairwise.complete.obs",
                  method = c("spearman"), alternative = c("two.sided"))

colnames(res1$p) <- rownames(res1$p) <- colnames(M) <- rownames(M) <- c("ERD", ":K['max,stem']", ":Psi[tlp]",
                                ":Psi['88,stem']", ":Psi[min]*'-'*Psi['88,stem']")

cor.k <- round(M[which(colnames(M) == ":K['max,stem']"), which(rownames(M) == "ERD")], 2)
cor.tlp <- round(M[which(colnames(M) == ":Psi[tlp]"), which(rownames(M) == "ERD")], 2)
cor.stem.88 <- round(M[which(colnames(M) == ":Psi['88,stem']"), which(rownames(M) == "ERD")], 2)
cor.hsm <- round(M[which(colnames(M) ==  ":Psi[min]*'-'*Psi['88,stem']"), which(rownames(M) == "ERD")], 2)

p.k <- round(res1$p[which(colnames(res1$p) == ":K['max,stem']"), which(rownames(res1$p) == "ERD")], 2)
p.tlp <- round(res1$p[which(colnames(res1$p) == ":Psi[tlp]"), which(rownames(res1$p) == "ERD")], 2)
p.stem.88 <- round(res1$p[which(colnames(res1$p) == ":Psi['88,stem']"), which(rownames(res1$p) == "ERD")], 2)
p.hsm <- round(res1$p[which(colnames(res1$p) ==  ":Psi[min]*'-'*Psi['88,stem']"), which(rownames(res1$p) == "ERD")], 2)

jpeg(file.path(figures.folder, "erd.stem.traits_cor.matrix_insig_blank.jpeg"),
     width = 4.5, height = 4.5, units = "in", pointsize = 10,
     quality = 100, res = 300)
corrplot(M, type = "upper", order = "original",
         tl.col = "black", tl.srt = 0, tl.offset = 0.8, tl.cex = 1,
         # p.mat = res1$p, sig.level = 0.05,
         ## show only insignificant p-values >= 0.05
         insig = "blank",
         ## if cor coefficient were to be shown
         addCoef.col = "gray",
         diag = FALSE, cl.ratio = .2, cl.align = "l")
dev.off()

jpeg(file.path(figures.folder, "erd.stem.traits_cor.matrix_all_p-val.jpeg"),
     width = 4.5, height = 4.5, units = "in", pointsize = 10,
     quality = 100, res = 300)
# par(mfrow=c(1,1),omi=c(0,0,0,0),mai=c(0.1,0.1,0.1,0.1))
corrplot(M, type = "upper", order = "original",
         tl.col = "black", tl.srt = 0, tl.offset = 0.8, tl.cex = 1,
         p.mat = res1$p,
         ## show only insignificant p-values >= 0.05
         insig = "p-val", sig.level = -1, pch.col = "grey90",
         diag = FALSE, cl.ratio = .2, cl.align = "l", win.asp= 1)
graphics.off()
chart.erd.pairs <- ggpairs(
  erd.pairs %>% dplyr::select(-sp),
  upper = list(continuous = wrap('cor', method = "spearman")),
  # upper = list(continuous = wrap(cor_func.2,
  #                                method = 'spearman', symbol = expression('\u03C1 ='))),
  lower = list(
    continuous = function(data, mapping, ...) {
      ggally_smooth_lm(data = data, mapping = mapping) +
        theme(panel.background = element_blank())
    }
  ),
  diag = list(
    continuous = function(data, mapping, ...) {
      ggally_densityDiag(data = data, mapping = mapping) +
        theme(panel.background = element_blank())
    }
  ),
  labeller = "label_parsed"
) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  theme(axis.text.x = element_text(
    face = "plain",
    angle = 90,
    vjust = 1,
    hjust = 1
  ))
ggsave(file.path(figures.folder, paste0("erd.stem.traits_cor.chart.jpeg")),
       plot = chart.erd.pairs, height = 5.5, width = 5.5, units ='in')

#****************************
traits.labels.cor.p <- traits.labels.select %>%
  mutate(r = c(cor.k, cor.tlp, cor.stem.88, cor.hsm),
                                   p = c(p.k, p.tlp, p.stem.88, p.hsm),
                                   x.rp = c(7.8, 7.8, 7.8, 0.4),
                                   y.r = c(1, -2.0, -3.7, -0.7),
                                   y.p = c(0, -2.3, -4.3, -1.2))
depth.traits.select.plot <- ggplot(erd.stem.traits.only.lab,
                                   aes(x = depth, y = value)) +
  geom_smooth(method = "lm", formula = formula) +
  geom_text(data = traits.labels.select, aes(x = x, y = y,
                                             label = panel, group = trait.plot),
            fontface = "bold", vjust = "inward", hjust = "inward") +
  geom_text(data = traits.labels.cor.p, aes(x = x.rp, y = y.r,
                                             label = sprintf('italic("r")~"="~%.2f', r),
            group = trait.plot), parse = TRUE, vjust = "inward", hjust = "inward") +
  geom_text(data = traits.labels.cor.p, aes(x = x.rp, y = y.p,
                                            label = sprintf('italic(p)~"="~%.2f', p),
                                            group = trait.plot), parse = TRUE, vjust = "inward", hjust = "inward") +
  geom_errorbarh(aes(xmax = depth + depth.se, xmin = depth - depth.se), size = 0.2) +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 0.8, size = 2.5) +
  coord_cartesian(xlim = c(0, max(erd.stem.traits.only.lab$depth) + 0.5)) +
  xlab("Effective Rooting Depth (m)") + ylab("") +
  facet_wrap(trait.plot ~ ., scales = "free_y", labeller = label_parsed,
             strip.position = 'left') +
  theme(strip.placement = "outside", panel.spacing.x = unit(0, "lines"),
        strip.text.y.left = element_text(size = 10, angle = 90, vjust = -1),
        plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
ggsave(file.path(figures.folder, paste0("erd.stem.traits.tiff")),
       plot = depth.traits.select.plot, height = 4.5, width = 5.5, units ='in')
ggsave(file.path(figures.folder, paste0("erd.stem.traits.jpeg")),
       plot = depth.traits.select.plot, height = 4.5, width = 5.5, units ='in')
#****************************
### Kmax vs. growth----
#****************************
stem.k.gr <- erd.stem.traits %>% left_join(demo.sp, by = "sp") %>%
  left_join(traits.labels.select %>% dplyr::select(trait, trait.plot), by = "trait") %>%
  droplevels()
stem.traits.sp <- length(erd.stem.traits$trait[erd.stem.traits$trait == "TLP"])
grate.adult.stem.traits.plot <- ggplot(stem.k.gr, aes(y = grate.adult, x = value)) +
  # geom_smooth(method = "lm", formula = formula) +
  facet_wrap(. ~ trait.plot, scales = "free_x", labeller = label_parsed, strip.position = 'bottom') +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 1, size = 2.5) +
  ylab(expression("Growth Rate (cm year"^-1*")")) + xlab("") +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.87, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 3) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = ifelse(p.value < 0.001, sprintf('italic(p)~"< 0.001"'),
                                     sprintf('italic(p)~"="~%.2f',stat(p.value)))),
                  parse = TRUE, npcx = 0.87, npcy = 0.8, size = 3) +
  theme(strip.placement = "outside", panel.spacing.y = unit(-0.5, "lines"),
        strip.text.x = element_text(size = 12, vjust = 2.5),
        plot.margin = margin(0.2, 0.2, -1, 0.2, "cm"))
ggsave(file.path(figures.folder, paste0("grate.adult.stem.traits.tiff")),
       plot = grate.adult.stem.traits.plot, height = 4.5, width = 3.8, units ='in')


#****************************
## Mortality rates vs. ERD----
#****************************
lm.mrate.mean <- lm(mrate ~ rdi.gr, data = mrate.depth.mean)
mrate.mean.cf <- as.numeric(round(lm.mrate.mean$coefficients, 2))
mrate.mean.r2 <- round(summary(lm.mrate.mean)$r.squared, 2)
mrate.mean.n <- length(lm.mrate.mean$residuals)
mrate.mean.p <- ifelse(broom::glance(lm.mrate.mean)$p.value < 0.001,
                   paste0("< 0.001"), paste0("= ", signif(broom::glance(lm.mrate.mean)$p.value, 2)))

mrate.plot.15 <- ggplot(mrate.depth.mean,
                        aes(y = mrate, x = rdi.gr)) +
  geom_smooth(method = "lm", formula = formula) +
  geom_errorbar(aes(ymin = mrate - mrate.se, ymax = mrate + mrate.se), width = 0.15, size = 0.1) +
  geom_errorbarh(aes(xmax = rdi.gr + depth.se, xmin = rdi.gr - depth.se), height = 0.15, size = 0.1) +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 1, size = 2.5) +
  ylab(expression('Mean Mortality Rate (%'*~'year'^-1*')')) +
  xlab("Effective Rooting Depth (m)") +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.95, npcy = 0.97, rr.digits = 2,
               formula = formula, parse = TRUE, size = 6) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = ifelse(p.value < 0.001, sprintf('italic(p)~"< 0.001"'),
                                     sprintf('italic(p)~"="~%.2f',stat(p.value)))),
                  parse = TRUE, npcx = 0.95, npcy = 0.84, size = 6) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by rdi.gr.tiff")),
       plot = mrate.plot.15, height = 3.1, width = 3.1, units='in')
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by rdi.gr.jpeg")),
       plot = mrate.plot.15, height = 3.1, width = 3.1, units='in')

mrate.plot.15.sub <- mrate.plot.15 %+%
  subset(mrate.depth.mean, sp %in% erd.stem.traits.sp)
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by rdi.gr_only_with_stem_traits.tiff")),
       plot = mrate.plot.15.sub, height = 3.1, width = 3.1, units = 'in')

mrate.plot.15.evg <- mrate.plot.15 %+% subset(mrate.depth.mean, deciduous == "E")
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by rdi.gr_evergreen.tiff")),
       plot = mrate.plot.15.evg, height = 3.1, width = 3.1, units = 'in')
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by rdi.gr_evergreen.jpeg")),
       plot = mrate.plot.15.evg, height = 3.1, width = 3.1, units = 'in')

#****************************
## Mortality vs. ERD by interval----
#****************************
# https://rpkgs.datanovia.com/ggpubr/reference/stat_cor.html
mrate.p.vals = sapply(unique(mrate.depth.select$censusint.m), function(i) {
  coef(summary(lm(mrate ~ rdi.gr, data=mrate.depth.select[mrate.depth.select$censusint.m==i, ])))[2,4]
})
mrate.p.vals.dat <- mrate.depth.select[mrate.depth.select$censusint.m %in%
                                         names(mrate.p.vals)[round(mrate.p.vals, 2) < 0.1 | round(mrate.p.vals, 1) == 0.1],]
mrate.r2.vals <- sapply(unique(mrate.depth.select$censusint.m), function(i) {
  round(summary(lm(mrate ~ rdi.gr, data=mrate.depth.select[mrate.depth.select$censusint.m==i, ]))$r.squared, 2)*100
})

m.p.90 <- round(mrate.p.vals["1985-90"], 1)
m.p.05 <- round(mrate.p.vals["2000-05"], 2)
m.p.10 <- round(mrate.p.vals["2005-10"], 1)

m.r2 <- c(mrate.r2.vals["1985-90"], mrate.r2.vals["2000-05"], mrate.r2.vals["2005-10"])

mrate.depth.select <- mrate.depth.select %>%
  transform(censusint.m.plot = factor(censusint.m,
                                      labels = c("1982-85", "*1985-90", "1990-95", "1995-00", "*2000-05", "*2005-10", "2010-15")))

mrate.plot.15.1 <- ggplot(mrate.depth.select, aes(y = mrate, x = rdi.gr)) +
  coord_cartesian(xlim = c(0, max(mrate.depth$rdi.gr, na.rm = TRUE))) +
  # scale_x_continuous(breaks = c(0, sort(unique(mrate.mfac.depth.gr.mean.mfac$depth)))) +
  geom_smooth(data = mrate.p.vals.dat, method = "lm", formula = formula) +
  geom_errorbarh(aes(xmax = rdi.gr + depth.se, xmin = rdi.gr - depth.se), height = 0.15, size = 0.1) +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 0.8, size = 2.5) +
  xlab("Effective Rooting Depth (m)")  +
  ylab(expression('Mortality Rate (%'*'yr'^-1*')')) +
  facet_grid(. ~ censusint.m.plot) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.95, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = ifelse(p.value < 0.001, sprintf('italic(p)~"< 0.001"'),
                                     sprintf('italic(p)~"="~%.2f',stat(p.value)))),
                  parse = TRUE, npcx = 0.95, npcy = 0.82, size = 4)
ggsave(file.path(paste0(figures.folder, "/mortality_by rdi.gr.tiff")),
       plot = mrate.plot.15.1, height = 2.5, width = 10, units = 'in')
ggsave(file.path(paste0(figures.folder, "/mortality_by rdi.gr.jpeg")),
       plot = mrate.plot.15.1, height = 2.5, width = 10, units = 'in')


# mrate.depth.select.sub <- subset(mrate.depth.select, sp %in% erd.stem.traits.sp)
# mrate.p.vals.sub = sapply(unique(mrate.depth.select.sub$censusint.m), function(i) {
#   coef(summary(lm(mrate ~ rdi.gr, data=mrate.depth.select.sub[mrate.depth.select.sub$censusint.m==i, ])))[2,4]
# })
# mrate.p.vals.dat.sub <- mrate.depth.select.sub[mrate.depth.select.sub$censusint.m %in% names(mrate.p.vals.sub )[mrate.p.vals.sub < 0.05],]
#
# mrate.plot.15.1.sub <- ggplot(mrate.depth.select.sub, aes(y = mrate, x = rdi.gr)) +
#   coord_cartesian(xlim = c(0, max(mrate.depth$rdi.gr, na.rm = TRUE))) +
#   geom_smooth(data = mrate.p.vals.dat.sub, method = "lm", formula = formula) +
#   geom_text(aes(label = avg.abund), size = 2, nudge_y = 0.4) +
#   geom_errorbarh(aes(xmax = rdi.gr + depth.se, xmin = rdi.gr - depth.se), height = 0.15, size = 0.1) +
#   geom_point(shape = 21, color = "white", fill = "black", alpha = 0.8, size = 2.5) +
#   xlab("Effective Rooting Depth (m)")  +
#   ylab(expression('Mortality Rate (%'*'year'^-1*')')) +
#   facet_grid(. ~ censusint.m) +
#   stat_poly_eq(aes(label = paste(..rr.label..)),
#                npcx = 0.95, npcy = 0.95, rr.digits = 2,
#                formula = formula, parse = TRUE, size = 4) +
#   stat_fit_glance(method = 'lm',
#                   method.args = list(formula = formula),
#                   geom = 'text_npc',
                    # aes(label = ifelse(p.value < 0.001, sprintf('italic(p)~"< 0.001"'),
                    # sprintf('italic(p)~"="~%.2f',stat(p.value)))),
                    # npcx = 0.95, npcy = 0.82, size = 4)
#
# ggsave(file.path(paste0(figures.folder, "/mortality_by rdi.gr_only_with_stem_traits.tiff")),
#        plot = mrate.plot.15.1.sub, height = 2.5, width = 10, units = 'in')


mrate.depth.select.evg <- subset(mrate.depth.select, deciduous == "E")
mrate.p.vals.evg <- sapply(unique(mrate.depth.select.evg$censusint.m), function(i) {
  coef(summary(lm(mrate ~ rdi.gr, data=mrate.depth.select.evg[mrate.depth.select.evg$censusint.m==i, ])))[2,4]
})
mrate.p.vals.dat.evg <- mrate.depth.select.evg[mrate.depth.select.evg$censusint.m %in%
                                                 names(mrate.p.vals.evg)[round(mrate.p.vals.evg, 2) < 0.1 | round(mrate.p.vals.evg, 1) == 0.1],]
mrate.r2.vals.evg <- sapply(unique(mrate.depth.select.evg$censusint.m), function(i) {
  round(summary(lm(mrate ~ rdi.gr, data=mrate.depth.select.evg[mrate.depth.select.evg$censusint.m==i, ]))$r.squared, 2)*100
})
m.evg.p.90 <- round(mrate.p.vals.evg["1985-90"], 1)
m.evg.p.05 <- round(mrate.p.vals.evg["2000-05"], 2)
m.evg.p.10 <- round(mrate.p.vals.evg["2005-10"], 2)

m.evg.r2.90 <- mrate.r2.vals.evg["1985-90"]
m.evg.r2.05 <- mrate.r2.vals.evg["2000-05"]
m.evg.r2.10 <- mrate.r2.vals.evg["2005-10"]

mrate.plot.15.1.evg <- ggplot(mrate.depth.select.evg, aes(y = mrate, x = rdi.gr)) +
  coord_cartesian(xlim = c(0, max(mrate.depth$rdi.gr, na.rm = TRUE))) +
  geom_smooth(data = mrate.p.vals.dat.evg, method = "lm", formula = formula) +
  geom_errorbarh(aes(xmax = rdi.gr + depth.se, xmin = rdi.gr - depth.se), height = 0.15, size = 0.1) +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 0.8, size = 2.5) +
  xlab("Effective Rooting Depth (m)")  +
  ylab(expression('Mortality Rate (%'*'year'^-1*')')) +
  facet_grid(. ~ censusint.m.plot) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.95, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = ifelse(p.value < 0.001, sprintf('italic(p)~"< 0.001"'),
                                     sprintf('italic(p)~"="~%.2f',stat(p.value)))),
                  parse = TRUE, npcx = 0.95, npcy = 0.82, size = 4)

ggsave(file.path(paste0(figures.folder, "/mortality_by rdi.gr_evergreen.tiff")),
       plot = mrate.plot.15.1.evg, height = 2.5, width = 10, units = 'in')
ggsave(file.path(paste0(figures.folder, "/mortality_by rdi.gr_evergreen.jpeg")),
       plot = mrate.plot.15.1.evg, height = 2.5, width = 10, units = 'in')

mrate.plot.15.1.evg.sub <- mrate.plot.15.1.evg %+%
  subset(mrate.depth.select.evg, !census %in% c(1982, 1985, 1990))
ggsave(file.path(paste0(figures.folder, "/mortality_by rdi.gr_evergreen_wo_1982-85-90censuses.jpeg")),
       plot = mrate.plot.15.1.evg.sub, height = 3, width = 9, units = 'in')

mrate.plot.15.1.evg.flipped <- ggplot(mrate.depth.select.evg, aes(x = mrate, y = rdi.gr)) +
  coord_cartesian(ylim = c(10, 0)) +
  # scale_y_reverse(breaks = seq(from = 0, to = 10, by = 2)) +
  scale_y_reverse(breaks = unique(mrate.depth.select.evg$rdi.gr)) +
  geom_smooth(data = mrate.p.vals.dat.evg, method = "lm", formula = formula) +
  geom_errorbar(aes(ymax = rdi.gr + depth.se, ymin = rdi.gr - depth.se), width = 0.15, size = 0.1) +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 0.8, size = 2.5) +
  ylab("Effective Rooting Depth (m)")  +
  xlab(expression('Mortality Rate (% '*'year'^-1*')')) +
  facet_grid(. ~ censusint.m.plot) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.95, npcy = 0.16, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = ifelse(p.value < 0.001, sprintf('italic(p)~"< 0.001"'),
                                     sprintf('italic(p)~"="~%.2f',stat(p.value)))),
                  npcx = 0.95, npcy = 0.03, size = 4) +
  theme(axis.title.y = element_text(size = 12))
ggsave(file.path(paste0(figures.folder, "/mortality_by rdi.gr_evergreen_flipped.jpeg")),
       plot = mrate.plot.15.1.evg.flipped, height = 2.7, width = 10, units = 'in')

mrate.plot.15.1.evg.flipped.sub <- mrate.plot.15.1.evg.flipped %+%
  subset(mrate.depth.select.evg, !census %in% c(1982, 1985, 1990))
ggsave(file.path(paste0(figures.folder, "/mortality_by rdi.gr_evergreen_wo_1982-85-90censuses_flipped.jpeg")),
       plot = mrate.plot.15.1.evg.flipped.sub, height = 3, width = 9, units = 'in')

#****************************
## mfac.rate vs. ERD ----
#****************************
mfac.plot.9.1 <- ggplot(mrate.mfac.depth.gr.mean.mfac,
                        aes(y = mfac.rate, x = depth)) +
  coord_cartesian(xlim = c(0, max(mrate.mfac.depth.gr.mean.mfac$depth, na.rm = TRUE) + 0.3)) +
  # scale_x_continuous(breaks = c(0, sort(unique(mrate.mfac.depth.gr.mean.mfac$depth)))) +
  geom_jitter(height = 0, width = 0.2, size = 2, shape = 21, alpha = 0.6, color = "black", aes(fill = sp), show.legend = FALSE) +
  xlab("Effective Rooting Depth (m)") +
  ylab(expression('Time below '*Psi['crit']*' (% yr'^-1*')'))
ggsave(file.path(paste0(figures.folder, "/mean_mfac vs. rdi.gr.jpeg")),
       plot = mfac.plot.9.1, height = 3.1, width = 3.5, units = 'in')

mfac.plot.9.1.sub <- mfac.plot.9.1 %+% subset(mrate.mfac.depth.gr.mean.mfac,
                                              sp %in% erd.stem.traits.sp)
ggsave(file.path(paste0(figures.folder, "/mean_mfac vs. rdi.gr_only_with_stem_traits.tiff")),
       plot = mfac.plot.9.1.sub, height = 3.1, width = 3.5, units = 'in')

mfac.plot.9.1.evg <- mfac.plot.9.1 %+% subset(mrate.mfac.depth.gr.mean.mfac,
                                              deciduous == "E")
ggsave(file.path(paste0(figures.folder, "/mean_mfac vs. rdi.gr_evergreen.jpeg")),
       plot = mfac.plot.9.1.evg, height = 3.1, width = 3.5, units = 'in')

mfac.plot.9.1.flipped <- ggplot(mrate.mfac.depth.gr.mean.mfac,
                        aes(x = mfac.rate, y = depth)) +
  scale_y_reverse(breaks = c(0, sort(unique(mrate.mfac.depth.gr.mean.mfac$depth)))) +
  # scale_x_continuous(breaks = c(0, sort(unique(mrate.mfac.depth.gr.mean.mfac$depth)))) +
  geom_jitter(height = 0.2, width = 0, size = 2, shape = 21, alpha = 0.6, color = "black", aes(fill = sp), show.legend = FALSE) +
  ylab("Effective Rooting Depth (m)") +
  xlab(expression('Time below '*Psi['crit']*' (% yr'^-1*')'))
ggsave(file.path(paste0(figures.folder, "/mean_mfac vs. rdi.gr_flipped.jpeg")),
       plot = mfac.plot.9.1.flipped, height = 3.1, width = 3.5, units = 'in')

mfac.plot.9.1.evg.flipped <- mfac.plot.9.1 %+% subset(mrate.mfac.depth.gr.mean.mfac,
                                              deciduous == "E")
ggsave(file.path(paste0(figures.folder, "/mean_mfac vs. rdi.gr_evergreen_flipped.jpeg")),
       plot = mfac.plot.9.1.evg, height = 3.1, width = 3.5, units = 'in')

#****************************
## mfac.rate vs. ERD by interval ----
#****************************
mrate.mfac.depth.select <- mrate.mfac.depth.select %>%
  arrange(sp, interval.num) %>%
  group_by(sp) %>%
  mutate(mfac.rate.prior = lag(mfac.rate, n = 1, default = NA)) %>%
  ungroup(sp)

mrate.mfac.depth.select.sp.mean <- mrate.mfac.depth.select %>%
  group_by(depth, censusint.m) %>%
  summarise(mean.mfac.rate = mean(mfac.rate, na.rm = TRUE),
            se = sd(mfac.rate, na.rm = TRUE)/n(), .groups = "drop")

mfac.plot.9.0.int <- ggplot(mrate.mfac.depth.select.sp.mean,
                                 aes(x = depth, y = mean.mfac.rate,
                                 ymax = mean.mfac.rate + se,
                                 ymin = mean.mfac.rate - se,
                                 fill = censusint.m)) +
  xlab("Effective Rooting Depth (m)") +
  ylab(expression('Time spent beyond '*Psi['crit']*' (%yr'^-1*')')) +
  geom_bar(position = position_dodge2(), stat="identity") +
  geom_errorbar(position = position_dodge2(.9, padding = .6)) +
  theme(legend.position = c(0.75, 0.65),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        legend.background = element_rect(fill = "transparent")) +
  scale_x_continuous(breaks = c(sort(unique(mrate.mfac.depth.gr.mean.mfac$depth)))) +
  # scale_fill_brewer(palette = "Greens") + #palette = "Spectral"
  guides(fill = guide_legend(title = "Census Interval"), override.aes = list(size = 1))
ggsave(file.path(paste0(figures.folder, "/mfac vs. rdi.gr.tiff")),
       plot = mfac.plot.9.0.int,  height = 3.2, width = 3.5, units = 'in')
ggsave(file.path(paste0(figures.folder, "/mfac vs. rdi.gr.jpeg")),
       plot = mfac.plot.9.0.int,  height = 3.2, width = 3.5, units = 'in')

mrate.mfac.depth.select.evg <- subset(mrate.mfac.depth.select,
                                      deciduous == "E") %>%
  group_by(depth, censusint.m) %>%
  summarise(mean.mfac.rate = mean(mfac.rate, na.rm = TRUE),
            se = sd(mfac.rate, na.rm = TRUE)/n(), .groups = "drop")

mfac.plot.9.0.int.evg <- mfac.plot.9.0.int %+% mrate.mfac.depth.select.evg

ggsave(file.path(paste0(figures.folder, "/mfac vs. rdi.gr_evergreen.jpeg")),
       plot = mfac.plot.9.0.int.evg,  height = 3.1, width = 3.5, units = 'in')

# dat1 <- mrate.mfac.depth.select %>% subset(deciduous == "E")
# dat1$new_value <- ifelse(dat1$mfac.rate<=5,dat1$mfac.rate,ifelse(dat1$mfac.rate<15,NA,dat1$mfac.rate-10))
# dat1 <- dat1[!is.na(dat1$new_value) ,]
#
# p = ggplot(dat1, aes(x=depth, y=new_value)) +   facet_grid(. ~ censusint.m)
# p + geom_jitter(aes(color = sp), show.legend = FALSE)+
#   theme(text = element_text(size=20),
#         axis.text.x = element_text(angle=90, vjust=1)) +
#   scale_y_continuous(breaks = 1:9, labels = c(1:5,"break",15:25))

#*********************************************
## Composite: mrate, droughts and mfac.rate ----
#*********************************************
composite.mort <- cowplot::plot_grid(mrate.plot.15.1.evg.flipped.sub, droughts.psi.heat,
                                     labels = c('a', 'b'),
                                    label_size = 14, ncol = 1, nrow = 2, rel_heights = c(1, 1.15), hjust = -2.5)
ggsave("erd_mrate_interval_droughts.jpeg", plot = composite.mort, path =
         file.path(figures.folder), device = "jpeg", height = 5, width = 9, units ='in')
ggsave("erd_mrate_interval_droughts.tiff", plot = composite.mort, path =
         file.path(figures.folder), device = "tiff", height = 5, width = 9, units ='in')

composite.mort.mfac <- cowplot::plot_grid(mrate.plot.15.1.evg.flipped.sub, droughts.psi.heat,
                                          mfac.plot.9.0.int.evg,
                                     labels = c('a', 'b', 'c'),
                                     label_size = 14, ncol = 1, nrow = 3,
                                     rel_heights = c(1.05, 1.2, 1), hjust = -2.5)
ggsave("erd_mrate_mfac_interval_droughts.jpeg", plot = composite.mort.mfac, path =
         file.path(figures.folder), device = "jpeg", height = 7, width = 9, units ='in')
ggsave("erd_mrate_mfac_interval_droughts.tiff", plot = composite.mort.mfac, path =
         file.path(figures.folder), device = "tiff", height = 7, width = 9, units ='in')

#****************************
## ERD vs. growth rates----
#****************************

pg.2 <- ggplot(mrate.depth.mean,
               aes(x = rdi.gr, y = grate)) +
  # geom_smooth(method = "lm", formula = formula) +
  geom_errorbar(aes(ymin = grate - grate.se, ymax = grate + grate.se), width = 0.15, size = 0.1) +
  geom_errorbarh(aes(xmax = rdi.gr + depth.se, xmin = rdi.gr - depth.se), height = 0.15, size = 0.1) +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 1, size = 2.5) +
  ylab(expression('Mean Growth Rate (cm yr'^-1*')')) +
  xlab("Effective Rooting Depth (m)") +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.95, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 6) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = ifelse(p.value < 0.001, sprintf('italic(p)~"< 0.001"'),
                                     sprintf('italic(p)~"="~%.2f',stat(p.value)))),
                  npcx = 0.95, npcy = 0.82, size = 6) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/adult_Growth_vs_rdi.gr.jpeg")), plot = pg.2, height = 3, width = 3, units='in')
ggsave(file.path(paste0(figures.folder, "/adult_Growth_vs_rdi.gr.tiff")), plot = pg.2, height = 3, width = 3, units='in')

#****************************
## K by Psi----
#****************************

pt1 <- cowplot::ggdraw() + cowplot::draw_image("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf/kmax_by_psi_color_by_SG100C_AVG_predicted_AB.tiff", scale = 1)
pt2 <- cowplot::ggdraw() + cowplot::draw_image("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf/kmax_by_psi_color_by_LMALAM_AVD_predicted_AB.tiff", scale = 1)

pb1 <- cowplot::ggdraw() + cowplot::draw_image("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf/std.k.spkmax_by_psi_color_by_SG100C_AVG_predicted_AB.tiff", scale = 1)
pb2 <- cowplot::ggdraw() + cowplot::draw_image("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf/std.k.spkmax_by_psi_color_by_LMALAM_AVD_predicted_AB.tiff", scale = 1)

plot.comm.plc <- cowplot::plot_grid(pt1, pt2, pb1, pb2, labels = c('a', 'b', 'c', 'd'),
                                    label_size = 14, ncol = 2, nrow = 2, rel_widths = c(1, 1))
ggsave("plot.comm.plc.tiff", plot = plot.comm.plc, path =
         file.path("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf"), device = "tiff", height = 4.4, width = 4.4, units ='in')
ggsave("plot.comm.plc.jpeg", plot = plot.comm.plc, path =
         file.path("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf"), device = "jpeg", height = 4.4, width = 4.4, units ='in')


pt1 <- cowplot::ggdraw() + cowplot::draw_image("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf/kmax_by_psi_color_by_SG100C_AVG_predicted_AB_for_data_sp.tiff", scale = 1)
pt2 <- cowplot::ggdraw() + cowplot::draw_image("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf/kmax_by_psi_color_by_LMALAM_AVD_predicted_AB_for_data_sp.tiff", scale = 1)

pb1 <- cowplot::ggdraw() + cowplot::draw_image("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf/std.k.spkmax_by_psi_color_by_SG100C_AVG_predicted_AB_for_data_sp.tiff", scale = 1)
pb2 <- cowplot::ggdraw() + cowplot::draw_image("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf/std.k.spkmax_by_psi_color_by_LMALAM_AVD_predicted_AB_for_data_sp.tiff", scale = 1)
plot.obs.plc <- cowplot::plot_grid(pt1, pt2, pb1, pb2, labels = c('a', 'b', 'c', 'd'),
                                   label_size = 14, ncol = 2, nrow = 2, rel_widths = c(1, 1))
ggsave("plot.obs.plc.tiff", plot = plot.obs.plc, path =
         file.path("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf"), device = "tiff", height = 4.4, width = 4.4, units ='in')
ggsave("plot.obs.plc.jpeg", plot = plot.obs.plc, path =
         file.path("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf"), device = "jpeg", height = 4.4, width = 4.4, units ='in')

#****************************
## ERD vs. Psi_crit----
#****************************

erd.data <- erd.data %>% left_join(data.model.AB.sub %>%
                                     dplyr::select(sp, Kmax, psi_kl50, psi_kl80, psi_kl20), by = "sp")
erd.p50.plot <- ggplot(erd.data, aes(y = depth, x = psi_kl50)) +
  geom_errorbar(aes(ymax = depth + depth.se, ymin = depth - depth.se), width = 0.01, size = 0.2) +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 0.8, size = 2.5) +
  # geom_point(shape = 21, color = "white", aes(fill = sp), alpha = 0.8, size = 2.5) +
  scale_y_reverse() +
  coord_cartesian(ylim = c(10, 0)) +
  ylab("Effective Rooting Depth (m)") + xlab(expression(Psi['crit']*~"(MPa)")) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.05, npcy = 0.15, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  aes(label = sprintf('italic(p)~"="~%.2f',stat(p.value))),
                  parse = TRUE,
                  npcx = 0.05, npcy = 0.05, size = 4)
ggsave(file.path(figures.folder, paste0("erd.p0L.jpeg")),
       plot = erd.p50.plot, height = 3, width = 3, units ='in')

#****************************
## plc curves: data vs. model----
#****************************
## those species for which curves are fit to data in solid lines, rest predicrted from model in dashed lines
obs.mod.plc <- obs.sp.vcurves.1 %>%
  subset(sp %in% erd.sp) %>%
  mutate(Fit = "Data")
obs.mod.plc <- obs.mod.plc %>%
  bind_rows(comm.sp.vcurves.1 %>%
              ## only those ERD species that for which data unavailable
              subset(sp %in% erd.sp[!erd.sp %in% unique(obs.mod.plc$sp)]) %>%
              mutate(Fit = "Model"))
sp.n.data <- length(unique(obs.mod.plc$sp[obs.mod.plc$Fit == "Data"]))
sp.n.model <- length(unique(obs.mod.plc$sp[obs.mod.plc$Fit == "Model"]))
obs.mod.vcurves_fits <- ggplot(obs.mod.plc , aes(x = psi, y = k.predict)) +
  ylim(c(0, 10)) +
  xlab("Leaf Water Potential (-MPa)") +
  ylab(expression(atop(italic('K')['max, leaf'], ""^(mmol*~m^-2*~s^-1*~MPa^-1)))) +
  theme(legend.position = c(0.8, 0.65),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent")) +
  scale_linetype_manual(values=c("solid", "longdash")) +
  guides(linetype = guide_legend(order = 2,
                                 title = NULL, direction = "horizontal", label.position = "left"))

obs.mod.plccurves.pre <- ggplot(obs.mod.plc, aes(x = psi, y = k.predict.percent)) +
  xlab("Leaf Water Potential (-MPa)") +
  ylab("Loss of Conductivity (%)") +
  # scale_color_gradient(low = "blue", high = "Red")
  theme(legend.position = c(0.8, 0.4),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent")) +
  guides(linetype = guide_legend(order = 2, title = ""))

obs.mod.plccurves.black <- obs.mod.plccurves.pre +
  geom_line(aes(group = sp, col = Fit), size = 0.25) +
  scale_color_manual(values=c('#E69F00','#999999')) +
  guides(color = guide_legend(override.aes = list(size = 3)))
  #scale_color_grey(start = 0.1, end = 0.7)
ggsave(plot = obs.mod.plccurves.black, file.path("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf/obs.mod_plccurves_fits_ggplot_no_col.tiff"),
       device = "tiff", height = 2.7, width = 3, units='in')
ggsave(plot = obs.mod.plccurves.black, file.path("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf/obs.mod_plccurves_fits_ggplot_no_col.jpeg"),
       device = "jpeg", height = 2.7, width = 3, units='in')


#****************************
## ab.table----
#****************************

ab.table.obs <- obs.data.model.AB %>%
  subset(sp %in% erd.sp) %>%
  mutate(Source = "Data") %>%
  rename(A = data.A, B = data.B) %>%
  dplyr::select(sp, A, B, Kmax, psi_kl50, Source)
ab.table <-  ab.table.obs %>% bind_rows(data.model.AB %>%
              ## only those ERD species that for which data unavailable
              subset(sp %in% erd.sp[!erd.sp %in% unique(ab.table.obs$sp)]) %>%
              mutate(Source = "Model") %>%
              rename(A = model.A, B = model.B) %>%
              dplyr::select(sp, A, B, Kmax, psi_kl50, Source)) %>%
  left_join(erd.sp.names %>% dplyr::rename(sp = Code), by = c("sp")) %>%
  dplyr::rename(Kmax_leaf = Kmax, Psi_50_leaf = psi_kl50) %>%
  mutate(Species = tolower(Species)) %>%
  dplyr::select(Genus, Species, Family, A, B, Kmax_leaf, Psi_50_leaf, Source) %>%
  mutate(Family = as.character(Family)) %>%
  mutate(A = round(A, 4), B = round(B, 4),
         Kmax_leaf = round(Kmax_leaf, 2), Psi_50_leaf = round(Psi_50_leaf, 2))

# Some familynames do not end in ceae. Correcting that
correct.family <- data.frame(misspelt = unique(ab.table$Family[-grep("aceae", ab.table$Family)]),
                             correct = c("Anacardiaceae", "Euphorbiaceae", "Flacourtiaceae", "Nyctaginaceae")) %>%
                             # correct = c("Menispermaceae", "Euphorbiaceae", "Anacardiaceae", "Hippocrateaceae", "Malpighiaceae",
                             #             "Flacourtiaceae", "Rhizophoraceae", "Melastomataceae", "Erythroxylaceae", "Nyctaginaceae", "Sterculiaceae",
                             #             "Lecythidaceae", "Chrysobalanaceae", "Convolvulaceae", "Simaroubaceae", "Elaeocarpaceae", "Staphyleaceae", "Myristicaceae")) %>%
  mutate(misspelt = as.character(misspelt),
         correct = as.character(correct))
## Which row in Family.sub matches with correct.family$misspelt
rows.to.replace <- which(ab.table$Family %in% correct.family$misspelt)
matched.rows <- match(ab.table$Family, correct.family$misspelt)
ab.table$Family[rows.to.replace] <- correct.family$correct[matched.rows[!is.na(matched.rows)]]

correct.family.2 <- data.frame(misspelt = unique(ab.table$Family[grep(":", ab.table$Family)]),
                               correct = c("Fabaceae:Papilionaceae", "Fabaceae:Mimosaceae", "Fabaceae:Papilionaceae")) %>%
  mutate(misspelt = as.character(misspelt),
         correct = as.character(correct))
rows.to.replace.2 <- which(ab.table$Family %in% correct.family.2$misspelt)
matched.rows.2 <- match(ab.table$Family, correct.family.2$misspelt)
ab.table$Family[rows.to.replace.2] <-
  correct.family.2$correct[matched.rows.2[!is.na(matched.rows.2)]]

#****************************
## species with lifespan----
#****************************
erd.sp.with.ll <- length(erd.sp[erd.sp %in% unique(bci.lifetime$sp[!is.na(bci.lifetime$lifetime)])])
erd.sp.wo.ll <- length(erd.sp[!erd.sp %in% unique(bci.lifetime$sp[!is.na(bci.lifetime$lifetime)])])

