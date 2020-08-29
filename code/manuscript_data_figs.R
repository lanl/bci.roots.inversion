#******************************************************
## output
#******************************************************
load(file = file.path(results.folder, "psi.stat.4.select.Rdata"))
load(file = file.path(results.folder, "ml.rsq.combine.Rdata"))
load(file = file.path(results.folder, "ml.rsq.combine.best.Rdata"))
load(file = file.path(results.folder, "mrate.depth.Rdata"))
load(file = file.path(results.folder, "mrate.mfac.depth.Rdata"))
load(file = file.path(results.folder, "erd.stem.traits.Rdata"))
load(file = file.path(results.folder, "depth.traits.kunert.Rdata"))
load(file = file.path(results.folder, "df.erd.to.plot.Rdata"))
load(file = file.path(results.folder, "data.model.AB.sub.Rdata"))

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
reverselog_trans <- function(base = exp(1)) {
  scales::trans_new(name = paste0("reverselog-", format(base)),
                    log_breaks(base = base),
                    domain = c(1e-100, Inf),
                    transform = function(x) -log(x, base),
                    inverse = function(x) base^(-x))
}

formula <- y ~ x
#****************************

###
ml.rsq.combine.best <- ml.rsq.combine.best %>%
  left_join(bci.traits %>% dplyr::select(sp, form1), by = "sp") %>%
  mutate(depth = as.numeric(depth))
erd.data <- ml.rsq.combine.best %>%
  subset(corr >= 0 & form1 == "T" &
           corr.func == "gr.Psi.VPD.multi" &
           R2 >= 0.1)
erd.iso <- erd.data %>% subset(sp != "guapst" &
                                 source == "Meinzer et al.1999 Fig. 4")

#******************************************************
### ERD Species names----
#******************************************************
erd.sp <- erd.data$sp
save(erd.sp, file = file.path("results", "erd.sp.Rdata"))
erd.sp.names <- bci.traits %>%
  subset(sp %in% erd.sp) %>%
  dplyr::rename(Code = sp, Genus = GENUS., Species = SPECIES., Family = FAMILY.) %>%
  dplyr::select(Code, Genus, Species, Family)
rownames(erd.sp.names) <- 1: nrow(erd.sp.names)
df.erd.to.plot <- df.erd.to.plot %>%
  left_join(erd.sp.names %>%
              mutate(s.names = paste0(substr(Genus, start = 1, stop = 1), ". ", tolower(Species)),
                     sp = tolower(Code)), by = "sp") %>%
  dplyr::select(sp, s.names, depth, depth.se) %>%
  transform(s.names = reorder(s.names, depth)) %>%
  droplevels()

erd.sp.plot <- ggplot(df.erd.to.plot,
                      aes(x = s.names, y = depth)) +
  geom_point(aes(color = s.names), show.legend = FALSE, size = 3) +
  geom_errorbar(aes(ymax = depth + depth.se, ymin = depth - depth.se), width = 0.2, size = 0.2) +
  ylab("Effective Rooting Depth (m)") + xlab("Species") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
        axis.text.y = element_text(face = "plain")) +
  scale_y_continuous(trans = reverselog_trans(10), breaks = unique(ml.rsq.combine$depth))
ggsave("ERD_by_sp_large_canopy.tiff",
       plot = erd.sp.plot, file.path(figures.folder), device = "tiff", height = 3.5, width = 5, units='in')
ggsave("ERD_by_sp_large_canopy.jpeg",
       plot = erd.sp.plot, file.path(figures.folder), device = "jpeg", height = 4.5, width = 5, units='in')

xylem.label <- expression('Xylem Sap '*delta~""^2*"H (\u2030)"*'')
ml.rsq.combine.sub <- ml.rsq.combine.best %>%
  mutate(depth = as.numeric(depth)) %>%
  subset(!sp %in% c("guapst") & !is.na(Xylem_sap_deltaD_permil.mean)) %>%
  left_join(bci.traits %>%
              dplyr::rename(Code = sp, Genus = GENUS., Species = SPECIES., Family = FAMILY.) %>%
              mutate(s.names = paste0(substr(Genus, start = 1, stop = 1), ". ", tolower(Species)),
                     sp = tolower(Code)) %>%
              dplyr::select(sp, s.names), by = "sp") %>%
  subset(source == "Meinzer et al.1999 Fig. 4") %>%
  droplevels()


formula = y~x
p4 <- ggplot(ml.rsq.combine.sub %>% subset(corr.func == "gr.Psi.VPD.multi"),
             aes(x = Xylem_sap_deltaD_permil, y = depth)) +
  geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + se,
                     xmin = Xylem_sap_deltaD_permil - se, color = s.names),
                 size = 0.5, height = 0.05, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, color = "black", size = 0.5, formula = formula) +
  ylab(expression("Water Uptake Depth (m)")) + xlab(xylem.label) +
  scale_y_continuous(trans=reverselog_trans(10), breaks = unique(ml.rsq.combine$depth)) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.9, npcy = 0.2, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.9, npcy = 0.1, size = 4) +
  geom_point(shape = 21, color = "white", aes(fill = s.names), alpha = 1, size = 3.5) +
  geom_errorbar(aes(ymax = depth + depth.se, ymin = depth - depth.se), color = "black",
                size = 0.3, width = 0.2) +
  guides(fill = guide_legend(title = "Species"), color = FALSE) +
  theme(legend.text = element_text(face = "italic", size = 8))
ggsave("psi.corr_best.depth_xylem_sap_deltaD_phenology_Meinzer_gr.Psi.VPD.jpeg",
       plot = p4, file.path(figures.folder), device = "jpeg", height = 3, width = 4.3, units = 'in')

ml.rsq.combine.sub <- ml.rsq.combine.sub %>%
  transform(models.plot1 = factor(corr.func, labels = c("A", "B", "C", "D", "E")))
erd.iso_sp_N_by_model <- ml.rsq.combine.sub %>%
  group_by(models.plot1) %>%
  summarise(N = n(), .groups = "drop_last")

p3.2 <- ggplot(ml.rsq.combine.sub,
               aes(x = Xylem_sap_deltaD_permil, y = depth)) +
  coord_cartesian(ylim = c(13, 0.3)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", size = 0.5, formula = formula) +
  geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + se,
                     xmin = Xylem_sap_deltaD_permil - se, color = s.names),
                 size = 0.5, height = 0.05) +
  geom_errorbar(aes(ymax = depth + depth.se, ymin = depth - depth.se, color = s.names), size = 0.5) +
  facet_wrap( ~ models.plot1, nrow = 2) +
  geom_text(data = erd.iso_sp_N_by_model, aes(x = -60, y = 0.3, label = paste0("n = ", N), group = models.plot1),
            vjust = "inward", hjust = "inward", inherit.aes = FALSE) +
  ylab(expression("Water Uptake Depth (m)")) + xlab(xylem.label) +
  scale_y_continuous(trans=reverselog_trans(10), breaks = unique(ml.rsq.combine$depth)) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.98, npcy = 0.15, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.98, npcy = 0.05, size = 4) +
  geom_point(shape = 21, color = "white", aes(fill = s.names), alpha = 1, size = 3) +
  guides(fill = guide_legend(title = "Species"), color = FALSE)
ggsave("psi.corr_best.depth_xylem_sap_deltaD_sp_color_Meinzer.jpeg",
       plot = p3.2, file.path(figures.folder), device = "jpeg", height = 5, width = 7, units = 'in')

#******************************************************

mrate.depth.select <- subset(mrate.depth, !is.na(rdi.gr) & avg.abund >= 20) %>%
  subset(sp %in% erd.sp) %>% droplevels()
mrate.mfac.depth.select <- subset(mrate.mfac.depth, !is.na(rdi.gr) &
                                    avg.abund >= 20 & depth == rdi.gr) %>%
  subset(sp %in% erd.sp) %>% droplevels() %>%
  dplyr::select(censusint.m, sp, depth, depth.se, avg.abund, trees, mfac, mrate,
         mean.mrate, diff.mrate, mean.grate, grate.se, size, deciduous)

# summarising acros interval
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
            mfac = sum(mfac, na.rm = TRUE), .groups = "drop_last")
# save.image("results/manuWorkSpace.RData")

## Minimum Soil water potential reached at depth 1.7 + CI

psi.1.7.min <- subset(psi.stat.4.select, depth == 1.7) %>%
  subset(median == min(median, na.rm = TRUE))

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
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_study_period_q2.5_depths_in_inverse_model.tiff",
       plot = plot.psi.stat.7.interval.q2.5.select, file.path(figures.folder), device = "tiff", height = 2.5, width = 6, units='in')
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_study_period_q2.5_depths_in_inverse_model.jpeg",
       plot = plot.psi.stat.7.interval.q2.5.select, file.path(figures.folder), device = "jpeg", height = 2.5, width = 6, units='in')

erd.stem.traits.only <- erd.stem.traits %>%
  left_join(df.erd.to.plot %>%
              dplyr::select(sp, depth, depth.se), by = "sp") %>%
  subset(!is.na(depth)) %>%
  droplevels()
erd.stem.traits.sp <- unique(erd.stem.traits.only$sp)

traits.labels.select <- data.frame(trait = factor(c("KmaxS", "TLP", "p88S", "HSM88S"),
                                                  levels = c("KmaxS", "TLP", "p88S", "HSM88S"), ordered = TRUE),
                                   panel = factor(c("A", "B", "C", "D"), levels = c("A", "B", "C", "D"), ordered = TRUE),
                                   x = c(0, 0, 0, 0),
                                   y = c(7, 0, 0.5, 2.4)) %>%
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


depth.traits.select.plot <- ggplot(erd.stem.traits.only.lab,
                                   aes(x = depth, y = value)) +
  geom_smooth(method = "lm", formula = formula) +
  geom_text(data = traits.labels.select, aes(x = x, y = y,
                                             label = panel, group = trait.plot), vjust = "inward", hjust = "inward") +
  geom_errorbarh(aes(xmax = depth + depth.se, xmin = depth - depth.se), size = 0.2) +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 0.8, size = 2.5) +
  # geom_point(shape = 21, color = "white", aes(fill = sp), alpha = 0.8, size = 2.5) +
  coord_cartesian(xlim = c(0, max(erd.stem.traits.only.lab$depth) + 0.5)) +
  xlab("Effective Rooting Depth (m)") + ylab("") +
  facet_wrap(trait.plot ~ ., scales = "free_y", labeller = label_parsed, strip.position = 'left') +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.95, npcy = 0.15, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.95, npcy = 0.05, size = 4) +
  theme(strip.placement = "outside", panel.spacing.x = unit(0, "lines"),
        strip.text.y.left = element_text(size = 10, angle = 90, vjust = -1),
        plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
ggsave(file.path(figures.folder, paste0("erd.stem.traits.tiff")),
       plot = depth.traits.select.plot, height = 4.5, width = 5.5, units ='in')

## Correlation chart
# Check correlations (as scatterplots), distribution and print correlation coefficient

erd.pairs <- erd.stem.traits.only.lab %>%
  select(sp, depth, trait.plot.chart, value) %>%
  pivot_wider(names_from = trait.plot.chart, values_from = value) %>%
  rename(ERD = depth)
chart.erd.pairs <- ggpairs(erd.pairs %>% select(-sp),
                       upper = list(continuous = wrap(cor_func,
                                                      method = 'spearman', symbol = expression('\u03C1 ='))),
                       lower = list(continuous = function(data, mapping, ...) {
                         ggally_smooth_lm(data = data, mapping = mapping) +
                           theme(panel.background = element_blank())}),
                       diag = list(continuous = function(data, mapping, ...) {
                         ggally_densityDiag(data = data, mapping = mapping) +
                           theme(panel.background = element_blank())}
                       ), labeller = "label_parsed")
ggsave(file.path(figures.folder, paste0("erd.stem.traits_cor.chart.jpeg")),
       plot = chart.erd.pairs + ggpairs.theme, height = 5.5, width = 5.5, units ='in')

# Kmax vs. growth
stem.k.gr <- erd.stem.traits %>% left_join(demo.sp, by = "sp") %>%
  left_join(traits.labels.select %>% select(trait, trait.plot), by = "trait") %>%
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
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.87, npcy = 0.8, size = 3) +
  theme(strip.placement = "outside", panel.spacing.y = unit(-0.5, "lines"),
        strip.text.x = element_text(size = 12, vjust = 2.5),
        plot.margin = margin(0.2, 0.2, -1, 0.2, "cm"))
ggsave(file.path(figures.folder, paste0("grate.adult.stem.traits.tiff")),
       plot = grate.adult.stem.traits.plot, height = 4.5, width = 3.8, units ='in')

leaf.k.gr <- depth.traits.kunert %>% left_join(demo.sp, by = "sp") %>%
  droplevels()
leaf.traits.sp <- nrow(leaf.k.gr)
grate.adult.leaf.traits.plot <- ggplot(leaf.k.gr %>%
                                         subset(trait %in% c("KmaxL") & !is.na(value)),
         aes(y = grate.adult, x = value)) +
  # geom_smooth(method = "lm", formula = formula) +
  facet_wrap(. ~ trait.plot.chart, scales = "free_x", labeller = label_parsed, strip.position = 'bottom') +
  geom_point(size = 3, alpha = 0.7) +
  ylab(expression("Growth Rate (cm year"^-1*")")) + xlab("") +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.87, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.87, npcy = 0.8, size = 4) +
  theme(strip.placement = "outside", panel.spacing.y = unit(-0.5, "lines"),
        strip.text.x = element_text(size = 12, vjust = 2.5),
        plot.margin = margin(0.2, 0.2, -0.25, 0.2, "cm"))
ggsave(file.path(figures.folder, paste0("grate.adult.leaf.traits.tiff")),
       plot = grate.adult.leaf.traits.plot, height = 3, width = 3, units ='in')

mfac.plot.15 <- ggplot(mrate.depth.mean,
                       aes(y = mrate, x = rdi.gr)) +
  geom_smooth(method = "lm", formula = formula) +
  geom_errorbar(aes(ymin = mrate - mrate.se, ymax = mrate + mrate.se), width = 0.15, size = 0.1) +
  geom_errorbarh(aes(xmax = rdi.gr + depth.se, xmin = rdi.gr - depth.se), height = 0.15, size = 0.1) +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 1, size = 2.5) +
  ylab(expression('Mean Mortality Rate (%'*~'year'^-1*')')) +
  xlab("Effective Rooting Depth (m)") +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.95, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 6) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.95, npcy = 0.82, size = 6) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by rdi.gr.tiff")),
       plot = mfac.plot.15, height = 3, width = 3, units='in')
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by rdi.gr.jpeg")),
       plot = mfac.plot.15, height = 3, width = 3, units='in')

mfac.plot.15.sub <- mfac.plot.15 %+% subset(mrate.depth.mean, sp %in% erd.stem.traits.sp)
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by rdi.gr_only_with_stem_traits.tiff")),
       plot = mfac.plot.15.sub, height = 3, width = 3, units = 'in')

pg.2 <- ggplot(mrate.depth.mean,
               aes(x = rdi.gr, y = grate)) +
  # geom_smooth(method = "lm", formula = formula) +
  geom_errorbar(aes(ymin = grate - grate.se, ymax = grate + grate.se), width = 0.15, size = 0.1) +
  geom_errorbarh(aes(xmax = rdi.gr + depth.se, xmin = rdi.gr - depth.se), height = 0.15, size = 0.1) +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 1, size = 2.5) +
  ylab(expression('Mean Growth Rate (cm year'^-1*')')) +
  xlab("Effective Rooting Depth (m)") +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.95, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 6) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.95, npcy = 0.82, size = 6) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/adult_Growth_vs_rdi.gr.jpeg")), plot = pg.2, height = 3, width = 3, units='in')
ggsave(file.path(paste0(figures.folder, "/adult_Growth_vs_rdi.gr.tiff")), plot = pg.2, height = 3, width = 3, units='in')


mfac.plot.15.1 <- ggplot(mrate.depth.select, aes(y = mrate, x = rdi.gr)) +
  coord_cartesian(xlim = c(0, max(mrate.depth$rdi.gr, na.rm = TRUE))) +
  geom_errorbarh(aes(xmax = rdi.gr + depth.se, xmin = rdi.gr - depth.se), height = 0.15, size = 0.1) +
  geom_smooth(method = "lm", formula = formula) +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 0.8, size = 2.5) +
  xlab("Effective Rooting Depth (m)")  +
  ylab(expression('Mortality Rate (%'*'year'^1*')')) +
  facet_grid(. ~ censusint.m) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.95, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.95, npcy = 0.82, size = 4)
ggsave(file.path(paste0(figures.folder, "/mortality_by rdi.gr.tiff")),
       plot = mfac.plot.15.1, height = 2.5, width = 10, units = 'in')
ggsave(file.path(paste0(figures.folder, "/mortality_by rdi.gr.jpeg")),
       plot = mfac.plot.15.1, height = 2.5, width = 10, units = 'in')

mfac.plot.15.1.sub <- mfac.plot.15.1 %+% subset(mrate.depth.select, sp %in% erd.stem.traits.sp & avg.abund >= 50) +
  geom_text(aes(label = avg.abund), size = 2, nudge_y = 0.4)
ggsave(file.path(paste0(figures.folder, "/mortality_by rdi.gr_only_with_stem_traits.tiff")),
       plot = mfac.plot.15.1.sub, height = 2.5, width = 10, units = 'in')

mfac.plot.9.0 <- ggplot(mrate.mfac.depth.gr.mean.mfac,
                        aes(x = mfac, y = depth)) +
  scale_y_continuous(trans = "rev_sqrt", breaks =
                       c(0, sort(unique(mrate.mfac.depth.gr.mean.mfac$depth)))) +
  geom_jitter(height = 0.1, width = 0, size = 2, shape = 21, alpha = 0.6, color = "black", aes(fill = sp), show.legend = FALSE) +
  ylab("Effective Rooting Depth (m)") +
  xlab(expression(atop('Time spent below '*Psi['crit'], '(Days over 1990-2015)')))
ggsave(file.path(paste0(figures.folder, "/mean_mfac vs. rdi.gr.tiff")),
       plot = mfac.plot.9.0, height = 3.5, width = 3.5, units = 'in')
ggsave(file.path(paste0(figures.folder, "/mean_mfac vs. rdi.gr.jpeg")),
       plot = mfac.plot.9.0, height = 3.5, width = 3.5, units = 'in')

mfac.plot.9.0.sub <- mfac.plot.9.0 %+% subset(mrate.mfac.depth.gr.mean.mfac,
                                              sp %in% erd.stem.traits.sp)
ggsave(file.path(paste0(figures.folder, "/mean_mfac vs. rdi.gr_only_with_stem_traits.tiff")),
       plot = mfac.plot.9.0.sub, height = 3.5, width = 3.5, units = 'in')


pt1 <- cowplot::ggdraw() + cowplot::draw_image("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf/kmax_by_psi_color_by_SG100C_AVG_predicted_AB.tiff", scale = 1)
pt2 <- cowplot::ggdraw() + cowplot::draw_image("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf/kmax_by_psi_color_by_LMALAM_AVD_predicted_AB.tiff", scale = 1)

pb1 <- cowplot::ggdraw() + cowplot::draw_image("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf/std.k.spkmax_by_psi_color_by_SG100C_AVG_predicted_AB.tiff", scale = 1)
pb2 <- cowplot::ggdraw() + cowplot::draw_image("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf/std.k.spkmax_by_psi_color_by_LMALAM_AVD_predicted_AB.tiff", scale = 1)

plot.comm.plc <- cowplot::plot_grid(pt1, pt2, pb1, pb2, labels = c('A', 'B', 'C', 'D'),
                                    label_size = 14, ncol = 2, nrow = 2, rel_widths = c(1, 1))
ggsave("plot.comm.plc.tiff", plot = plot.comm.plc, path =
         file.path("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf"), device = "tiff", height = 4.4, width = 4.4, units ='in')
ggsave("plot.comm.plc.jpeg", plot = plot.comm.plc, path =
         file.path("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf"), device = "jpeg", height = 4.4, width = 4.4, units ='in')


pt1 <- cowplot::ggdraw() + cowplot::draw_image("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf/kmax_by_psi_color_by_SG100C_AVG_predicted_AB_for_data_sp.tiff", scale = 1)
pt2 <- cowplot::ggdraw() + cowplot::draw_image("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf/kmax_by_psi_color_by_LMALAM_AVD_predicted_AB_for_data_sp.tiff", scale = 1)

pb1 <- cowplot::ggdraw() + cowplot::draw_image("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf/std.k.spkmax_by_psi_color_by_SG100C_AVG_predicted_AB_for_data_sp.tiff", scale = 1)
pb2 <- cowplot::ggdraw() + cowplot::draw_image("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf/std.k.spkmax_by_psi_color_by_LMALAM_AVD_predicted_AB_for_data_sp.tiff", scale = 1)
plot.obs.plc <- cowplot::plot_grid(pt1, pt2, pb1, pb2, labels = c('A', 'B', 'C', 'D'),
                                   label_size = 14, ncol = 2, nrow = 2, rel_widths = c(1, 1))
ggsave("plot.obs.plc.tiff", plot = plot.obs.plc, path =
         file.path("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf"), device = "tiff", height = 4.4, width = 4.4, units ='in')
ggsave("plot.obs.plc.jpeg", plot = plot.obs.plc, path =
         file.path("figures/PhenoDemoTraitsPsi/kmax_by_psi/Leaf"), device = "jpeg", height = 4.4, width = 4.4, units ='in')

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
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.05, npcy = 0.05, size = 4)
ggsave(file.path(figures.folder, paste0("erd.p0L.jpeg")),
       plot = erd.p50.plot, height = 3, width = 3, units ='in')


