###---------------------------------------------------------------------------------------------------------
### Parts of PhenoDemoTraits.R after running RED model, but not required for the manucript
### Could be run after running PhenoDemoTraits.R
###---------------------------------------------------------------------------------------------------------

load(file = file.path(results.folder, "chosen.model.Rdata"))
load(file = file.path(results.folder, "psi.stat.4.Rdata"))
load(file = file.path(results.folder, "psi.stat.4.select.Rdata"))

## tlp was used in an earlier version
# tlp.sp <-  traits.long %>% subset(trait == "TLP") %>% select(sp, value) %>% rename(tlp = value)
# n.tlp.sp <- length(tlp.sp$sp)
# tlp.sp <- tlp.sp %>%
#   mutate(psi_threshold = tlp*0.8)
# tlp.sp.ls <- split(tlp.sp, f = list(tlp.sp$sp), drop = TRUE)

grate.gfac <- gfac.interval.long <- vector(mode = "list", length = length(names.gfac))
names(gfac.interval.long) <- names(grate.gfac) <- names.gfac
for (i in 1:length(names.gfac)) {
  gfac.interval.long[[i]] <- gfac.interval[[i]] %>%
    pivot_longer(cols = c(-sp, -interval),
                 names_to = "depth", values_to = "gfac") %>%
    rename(interval.num =  interval) %>%
    mutate(depth = as.numeric(depth),
           interval.num =  as.numeric(as.character(interval.num)),
           size = "large")
  grate.gfac[[i]] <- growth.sub %>% rename(interval.num = interval) %>%
    left_join(gfac.interval.long[[i]], by = c("interval.num", "sp", "size")) %>%
    subset(!is.na(gfac))
}

grate.gfac.best <- dplyr::bind_rows(grate.gfac, .id = "corr.func") %>%
  mutate(censusint.m = recode(interval.num, `1` = "1982-85", `2` = "1985-90",
                              `3` = "1990-95", `4` = "1995-00", `5` = "2000-05", `6` = "2005-10", `7` = "2010-15")) %>%
  transform(corr.func = factor(corr.func, levels = names.gfac)) %>%
  unite(corr.func_sp_depth, corr.func, sp, depth, remove = FALSE) %>%
  group_by(sp, size, corr.func, depth) %>%
  mutate(std.gfac = scale(gfac),
         std.growth = scale(demo.rate)) %>%
  left_join(ml.rsq.combine.best %>%
              dplyr::select(sp, size, corr.func, R2),
            by = c("sp", "size", "corr.func"))
grate.gfac.best.sub <- grate.gfac.best %>%
  subset(#corr.func %in% names.gfac &
    sp %in% iso.1.3.join$sp[iso.1.3.join$source == "Meinzer et al.1999 Fig. 4"] &
      ## those that are likely leafless in Mar-Apr
      !sp %in% as.character(leafless_mar.apr$sp[leafless_mar.apr$leafless_in_mar_apr_from_notes == "Yes"]) &
      corr.func %in% unique(grate.gfac.best$corr.func)[grep("gr.", unique(grate.gfac.best$corr.func))])
grate.gfac.best.plot <- ggplot(grate.gfac.best.sub %>%
                                 subset(corr.func_sp_depth %in% ml.rsq.combine.best$corr.func_sp_depth) %>% droplevels(),
                               aes(x = censusint.m)) +
  geom_line(aes(y = std.growth, group = corr.func_sp_depth, color = sp, linetype = "Std.Growth"), size = 1) +
  geom_line(aes(y = std.gfac, group = corr.func_sp_depth, color = sp, linetype = "Std.Growth Factor"), size = 1) +
  facet_grid(corr.func ~ sp , scales = "free_y") +
  scale_color_brewer(palette = "Dark2") +
  xlab("Census Interval") + ylab("Std.Growth/Growth Factor") +
  guides(color = "none", linetype = guide_legend(order = 1, title = NULL, direction = "horizontal", label.position = "top", override.aes =
                                                   list(linetype = c("Std.Growth" = "solid", "Std.Growth Factor" = "dotted")))) +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_text(aes(label = round(R2, 1), x = "2005-10", y = 1.2))
ggsave("Std.Growth Vs. Std.Growth Factor.jpeg",
       plot = grate.gfac.best.plot, file.path(figures.folder), device = "jpeg", height = 7, width = 14, units='in')

grate.gfac.best.plot.depths <- ggplot(grate.gfac.best.sub, #%>%
                                      # subset(depth %in% c(0.1, 1.0, 2.0)) %>% droplevels(),
                                      aes(x = censusint.m)) +
  geom_line(aes(y = std.gfac, group = corr.func_sp_depth, color = depth), size = 1) +
  geom_line(aes(y = std.growth, group = corr.func_sp_depth,
                linetype = "Std.Growth"), size = 1, color = "red") +
  facet_grid(corr.func ~ sp , scales = "free_y") +
  # scale_color_discrete(name = "Depth (m)", drop = FALSE) +
  scale_color_continuous(name = "Depth (m)", trans = "reverse") +
  xlab("Census Interval") + ylab("Std.Growth/Growth Factor") +
  scale_linetype_manual(name = "", values = c("solid")) +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_text(aes(label = round(R2, 1), x = "2000-05", y = 1.2))
ggsave("Std.Growth Vs. Std.Growth Factor_depths.jpeg",
       plot = grate.gfac.best.plot.depths, file.path(figures.folder), device = "jpeg", height = 7, width = 14, units='in')

grate.gfac.best.plot.depths.best <- ggplot(grate.gfac.best.sub %>%
                                             subset(corr.func_sp_depth %in% ml.rsq.combine.best$corr.func_sp_depth),
                                           aes(x = censusint.m)) +
  geom_line(aes(y = std.gfac, group = corr.func_sp_depth, color = as.factor(depth)), size = 1) +
  geom_line(aes(y = std.growth, group = corr.func_sp_depth,
                linetype = "Std.Growth"), size = 1) +
  facet_grid(corr.func ~ sp , scales = "free_y") +
  scale_color_discrete(name = "Depth (m)", drop = FALSE) +
  # scale_color_continuous(name = "Depth (m)", trans = "reverse") +
  xlab("Census Interval") + ylab("Std.Growth/Growth Factor") +
  scale_linetype_manual(name = "", values = c("solid")) +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_text(aes(label = round(R2, 1), x = "2000-05", y = 1.2))
ggsave("Std.Growth Vs. Std.Growth Factor_depths.best.jpeg",
       plot = grate.gfac.best.plot.depths.best, file.path(figures.folder), device = "jpeg", height = 7, width = 14, units='in')

depth.rsq.isotopes <- ml.rsq.combine.best %>%
  group_by(corr.func, sp, size) %>%
  summarise_at(c("depth", "depth.se","Xylem_sap_deltaD_permil", "se"), mean, na.rm = TRUE, .groups = "drop_last") %>%
  ungroup(corr.func, sp, size)
save(depth.rsq.isotopes, file = file.path(results.folder, "depth.rsq.isotopes.Rdata"))

### Plot best correlated depth against isotopic data and traits-----

# ml.rsq.combine.best <- ml.rsq.combine.best %>%
#   left_join(traits.wide.hyd %>% select(sp, KmaxL, Panama.moist.pref:HSMLWP.TLP), by = "sp")

heat.rsq <- ggplot(ml.rsq.combine %>% droplevels(),
                   aes(y = deci_sp, x = as.factor(depth))) +
  geom_tile(aes(fill = corr)) +
  ylab("Species") + xlab("Depth (m)") +
  facet_wrap(. ~ corr.func, nrow = 1) +
  scale_fill_viridis_c(expression("Pearson's "*italic(rho)), direction = -1, option = "plasma") #+
ggsave("psi.corr_all.depths_phenology_heat_by_corr.func.jpeg",
       plot = heat.rsq, file.path(figures.folder), device = "jpeg", height = 5, width = 12, units='in')

# theme(axis.text.y = element_text(angle = 90, vjust = 0.5)) +
# scale_x_continuous(breaks = soil.depths[c(1,8:13)])
heat.best.rsq <- heat.rsq %+% subset(ml.rsq.combine.best, R2 >= 0.2)
ggsave("psi.corr_best.depth_phenology_heat_by_corr.func.jpeg",
       plot = heat.best.rsq, file.path(figures.folder), device = "jpeg", height = 5, width = 12, units='in')

xylem.label <- expression('Xylem Sap '*delta~""^2*"H (\u2030)"*'')
ml.rsq.combine.best <- ml.rsq.combine.best %>% left_join(bci.traits %>% dplyr::select(sp, form1), by = "sp") %>%
  mutate(depth = as.numeric(depth))
ml.rsq.combine.sub <- ml.rsq.combine.best %>%
  subset(form1 == "T" &
           #corr.func %in% names.gfac[1:5] &
           !sp %in% c("guapst") &
           # !sp %in% c("guapst", "alsebl") &
           #corr.func == "gr.Psi.Rad.VPD" &
           !is.na(Xylem_sap_deltaD_permil.mean)) %>%
  droplevels()

formula = y~x
p0 <- ggplot(ml.rsq.combine.sub,
             aes(x = Xylem_sap_deltaD_permil.mean, y = depth)) + #HSMTLP.80L)) +
  geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil.mean + se.mean,
                     xmin = Xylem_sap_deltaD_permil.mean - se.mean, color = deciduousness),
                 size = 0.5, height = 0.05) +
  geom_errorbar(aes(ymax = depth + depth.se, ymin = depth - depth.se, color = deciduousness), size = 0.5, height = 0.05) +
  facet_wrap( ~ corr.func, nrow = 1) +
  geom_text(aes(x =  Xylem_sap_deltaD_permil.mean, y = depth, label = sp, color = deciduousness), nudge_y = 0.1, nudge_x = 0.2,
            size = 4, show.legend = FALSE) +
  ylab(expression("Best Correlated Depth (m)")) + xlab(xylem.label) +
  scale_y_continuous(trans = reverselog_trans(10), breaks = unique(ml.rsq.combine$depth)) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.05, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.05, npcy = 0.85, size = 4) +
  geom_point(size = 1, show.legend = TRUE, aes(color = deciduousness)) +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.direction = "horizontal") + #, strip.text = element_blank()) +
  scale_color_brewer(palette = "Dark2")
ggsave("psi.corr_best.depth_xylem_sap_deltaD_phenology_mean_isotope_source.jpeg",
       plot = p0, file.path(figures.folder), device = "jpeg", height = 3.5, width = 8, units = 'in')

p1 <- ggplot(ml.rsq.combine.sub,
             aes(x = Xylem_sap_deltaD_permil, y = depth)) + #HSMTLP.80L)) +
  geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + se,
                     xmin = Xylem_sap_deltaD_permil - se, color = deciduousness),
                 size = 0.5, height = 0.05) +
  geom_errorbar(aes(ymax = depth + depth.se, ymin = depth - depth.se, color = deciduousness),
                size = 0.5, height = 0.05) +
  facet_wrap( ~ corr.func, nrow = 1) +
  geom_text(aes(x =  Xylem_sap_deltaD_permil, y = depth, label = sp, color = deciduousness), nudge_y = 0.1, nudge_x = 0.2,
            size = 4, show.legend = FALSE) +
  ylab(expression("Best Correlated Depth (m)")) + xlab(xylem.label) +
  scale_y_continuous(trans=reverselog_trans(10), breaks = unique(ml.rsq.combine$depth)) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.9, npcy = 0.2, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.9, npcy = 0.1, size = 4) +
  geom_point(size = 2, show.legend = TRUE, aes(color = deciduousness)) +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.direction = "horizontal") + #, strip.text = element_blank()) +
  scale_color_brewer(palette = "Dark2")
ggsave("psi.corr_best.depth_xylem_sap_deltaD_phenology_two_isotope_sources.jpeg",
       plot = p1, file.path(figures.folder), device = "jpeg", height = 3.5, width = 8, units = 'in')
p2 <- p1 %+% subset(ml.rsq.combine.sub, source == "Meinzer et al.1999 Fig. 4")
ggsave("psi.corr_best.depth_xylem_sap_deltaD_phenology_Meinzer.jpeg",
       plot = p2, file.path(figures.folder), device = "jpeg", height = 3.5, width = 8, units = 'in')

p3 <- ggplot(ml.rsq.combine.sub %>% subset(source == "Meinzer et al.1999 Fig. 4"),
             aes(x = Xylem_sap_deltaD_permil, y = depth)) + #HSMTLP.80L)) +
  coord_cartesian(ylim = c(13, 0.3)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", size = 0.5, formula = formula) +
  geom_errorbarh(aes(xmax = Xylem_sap_deltaD_permil + se,
                     xmin = Xylem_sap_deltaD_permil - se, color = deciduousness),
                 size = 0.5, height = 0.1) +
  geom_errorbar(aes(ymax = depth + depth.se, ymin = depth - depth.se, color = deciduousness), size = 0.5, height = 0.05) +
  facet_wrap( ~ corr.func, nrow = 1) +
  geom_text(aes(x =  Xylem_sap_deltaD_permil, y = depth, label = sp, color = deciduousness), nudge_y = 0.1, nudge_x = 0.2,
            size = 4, show.legend = FALSE) +
  # position=position_jitter(width=ifelse(ml.rsq.combine.sub$sp=='cordal',1,0),
  #                        height=ifelse(ml.rsq.combine.sub$sp=='cordal',1,0))
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
  geom_point(size = 2, show.legend = TRUE, aes(color = deciduousness)) +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.direction = "horizontal") +
  scale_color_brewer(palette = "Dark2")
ggsave("psi.corr_best.depth_xylem_sap_deltaD_phenology_Meinzer.jpeg",
       plot = p3, file.path(figures.folder), device = "jpeg", height = 4, width = 8, units = 'in')


ml.rsq.combine.sub <- ml.rsq.combine.sub %>% transform(models.plot1 = factor(corr.func,
                                                                             labels = c("A", "B", "C", "D")),
                                                       models.plot2 = factor(corr.func,
                                                                             labels = c(expression(italic(K[italic('leaf')])),
                                                                                        expression(italic(K[italic('leaf')]*'VPD')),
                                                                                        expression(italic(K[italic('leaf')]*'LeafCover')),
                                                                                        expression(italic(K[italic('leaf')]*'VPD'*'LeafCover')))))


g3 <- ggplot(ml.rsq.combine.best %>% subset(R2 >= 0.1 & !duplicated(sp) & !is.na(depth)),
             aes(x = deci_sp.plot, y = depth)) +
  facet_wrap(. ~ corr.func) +
  geom_col(aes(fill = deciduousness)) +
  guides(fill = guide_legend(title = "")) +
  theme(legend.position = "top") +
  ylab("Depth (m)") + xlab("Species") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(trans = reverselog_trans(10), breaks = unique(ml.rsq.combine$depth))
g4 <- g3 + coord_flip() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))
ggsave("psi.corr_best.depth_phenology.jpeg",
       plot = g4, file.path(figures.folder), device = "jpeg", height = 6, width = 9, units = 'in')


### Plot against hydraulic traits------

## Soil preference vs traits
##
traits <- traits %>% left_join(depth.rsq.isotopes %>% ungroup() %>%
                                 subset(corr.func == chosen.model) %>%
                                 select(sp, depth, depth.se), by = "sp") %>%
  left_join(iso.1.3.join %>% subset(source == "Meinzer et al.1999 Fig. 4") %>%
              dplyr::select(sp, Xylem_sap_deltaD_permil, se), by = "sp") %>%
  droplevels()

traits.labels.table.2 <- data.frame(trait = factor(c("depth", "Xylem_sap_deltaD_permil",
                                                     "KmaxL", "lwp.min_Predawn", "lwp.min_Diurnal", "TLP", "p50L", "p80L",
                                                     "HSMLWP.TLP", "HSMLWP.50L", "HSMTLP.50L",
                                                     "HSMLWP.80L", "HSMTLP.80L",
                                                     "Panama.moist.pref", "Plot.swp.pref", "SG100C_AVG", "Chl"),
                                                   levels = c("depth", "Xylem_sap_deltaD_permil",
                                                              "KmaxL", "lwp.min_Predawn", "lwp.min_Diurnal", "TLP", "p50L", "p80L",
                                                              "HSMLWP.TLP", "HSMLWP.50L", "HSMTLP.50L",
                                                              "HSMLWP.80L", "HSMTLP.80L",
                                                              "Panama.moist.pref", "Plot.swp.pref", "SG100C_AVG", "Chl"), ordered = TRUE)) %>%
  transform(trait.plot = factor(trait, labels = c(expression(Depth[italic('Rsq')]), expression(italic(delta)^2*H[Xylem]),
                                                  expression(italic(K)['max, leaf']), expression(Psi[predawn]), expression(Psi[min]),
                                                  expression(Psi[tlp]), expression(Psi['50, leaf']), expression(Psi['80, leaf']),
                                                  expression(Psi[min]*' - '*Psi[tlp]),
                                                  expression(Psi[min]*' - '*Psi['50, leaf']),
                                                  expression(Psi[tlp]*' - '*Psi['50, leaf']),
                                                  expression(Psi[min]*' - '*Psi['80, leaf']),
                                                  expression(Psi[tlp]*' - '*Psi['80, leaf']),
                                                  expression('Panama'[wet]), expression('Plot'[wet]), expression('SG'[100*~degree*C]), "LMA")),
            trait.plot.chart = factor(trait, labels = c(expression(Depth[italic('Rsq')]), expression(italic(delta)^2*H[Xylem]),
                                                        expression(italic(K)['max, leaf']), expression(Psi[predawn]), expression(Psi[min]),
                                                        expression(Psi[tlp]), expression(Psi['50,leaf']), expression(Psi['80,leaf']),
                                                        expression(Psi[min]*'-'*Psi[tlp]),
                                                        expression(Psi[min]*'-'*Psi['50,leaf']),
                                                        expression(Psi[tlp]*'-'*Psi['50,leaf']),
                                                        expression(Psi[min]*'-'*Psi['80,leaf']),
                                                        expression(Psi[tlp]*'-'*Psi['80,leaf']),
                                                        expression('Panama'[wet]), expression('Plot'[wet]), expression('SG'[100*~degree*C]), "LMA")))


## Kunert traits
traits.long <- traits %>% select(-DeciLvl) %>%
  gather(trait, value, -sp, -deciduousness, -deciduous, -form1) %>%
  subset(deciduousness != "NA") %>%
  droplevels()

kruskal.list <- list()
for(i in unique(traits.long$trait)) {
  xx <- traits.long %>% subset(trait == i)
  kruskal.list[[i]] <- cbind(trait = i, kruskal(xx$value, xx$deciduousness, alpha = 0.1, group=TRUE, p.adj="bonferroni")$groups,
                             deciduousness = rownames(kruskal(xx$value, xx$deciduousness, alpha = 0.1, group=TRUE, p.adj="bonferroni")$groups))
}
traits.kruskal.labels <- do.call(rbind.data.frame, kruskal.list)
head(traits.kruskal.labels)

traits.labels <- traits.kruskal.labels
traits.labels.data <- traits.labels %>%
  left_join(traits.long %>% group_by(trait) %>%
              summarise(value = max(value, na.rm = TRUE)), by = c("trait"), .groups = "drop_last") %>%
  subset(deciduousness != "NA") %>%
  droplevels() %>%
  transform(deciduousness = factor(deciduousness,
                                   levels = c("Evergreen", "Brevideciduous",
                                              "Facultative Deciduous", "Obligate Deciduous"), ordered = TRUE)) %>%
  left_join(traits.labels.table.2 %>% select(trait, trait.plot), by = "trait")

traits.long <- traits.long %>%
  left_join(traits.labels.table.2, by = "trait")

## Kunert traits species wise for sp in hyd.traits----
select.traits <- c("lwp.min_Predawn", "lwp.min_Diurnal", "TLP", "KmaxS", "p50S", "p88S",
                   "HSMTLP", "HSM50S", "HSM88S", "HSMTLP.50S", "HSMTLP.88S")

traits.long <- traits.long %>%
  unite("deci_sp", deciduous, sp, remove = FALSE) %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(deciduousness)]), ordered=TRUE),
         deci_sp.plot = factor(deci_sp, levels=unique(deci_sp[order(deciduousness)]), ordered=TRUE))
# just for sp with hyd.traits, but traits.long does not have all those sp, and hab preference and WSG traits will be missed
## so beginning with those other traits
traits.wide <- traits.long %>% select(-trait.plot, -trait.plot.chart) %>%
  pivot_wider(names_from = trait, values_from = value)
traits.long.hyd <- data.frame(sp = unique(c(hyd$sp, traits$sp))) %>%
  full_join(bci.traits %>% select(sp, form1, SG100C_AVG), by = "sp") %>%
  full_join(deci %>% select(-sp4), by = "sp") %>%
  subset(sp %in% unique(c(depth.rsq.isotopes$sp, hyd$sp, traits$sp))) %>%
  left_join(traits.wide %>% select(-form1, -deciduous, -deciduousness,
                                   -SG100C_AVG, -Panama.moist.pref, -Plot.swp.pref), by = "sp") %>%
  pivot_longer(cols = c(-sp, -form1, -deciduous, -deciduousness, -DeciLvl,
                        -deciduousness.label, -sp.plot, -deci_sp, -deci_sp.plot),
               names_to = "trait", values_to = "value") %>%
  unite("deci_sp", deciduous, sp, remove = FALSE) %>%
  left_join(traits.labels.table.2, by = "trait") %>%
  mutate(sp.plot = factor(sp, levels = unique(sp[order(deciduousness)]), ordered=TRUE),
         deci_sp.plot = factor(deci_sp, levels=unique(deci_sp[order(deciduousness)]), ordered=TRUE)) %>%
  droplevels()

save(traits.long, file = file.path(results.folder, "kunert.traits.key.long_depth_isotopes.RData"))
save(traits.long.hyd, file = file.path(results.folder, "kunert.traits.key.long_in_Wolfe_traits_species_list_depth_isotopes.RData"))

d1 <- ggplot(traits.long.hyd %>% subset(!is.na(deciduousness))) +
  facet_wrap(. ~  trait.plot, scales = "free_x", labeller = label_parsed, nrow = 2) +
  geom_col(aes(x = deci_sp.plot, y = value,
               fill = deciduousness),
           position = position_dodge2(width = 0.9, preserve = "single")) +
  theme(axis.text.y = element_text(face = "plain", size = 8)) +
  coord_flip() +
  guides(fill = guide_legend(title = "Deciduousness")) +
  xlab("Deciduousness") + ylab("Value")
ggsave(file.path(figures.folder, paste0("Kunert_traits_vs_deciduousness_sp_bar_depth_isotopes_all.sp.with.depth.jpeg")),
       plot = d1, height = 7, width = 14, units ='in')
d2 <- d1 %+% subset(traits.long.hyd, !is.na(deciduousness) & sp %in% traits$sp)
ggsave(file.path(figures.folder, paste0("Kunert_traits_vs_deciduousness_sp_bar_depth_isotopes.jpeg")),
       plot = d2, height = 7, width = 14, units ='in')

## Correlation chart

# Check correlations (as scatterplots), distribution and print correlation coefficient
select.traits.1 <- c("TLP", "KmaxS","p50S", "p88S", "HSMTLP", "HSM50S",
                     "lwp.min_Predawn", "lwp.min_Diurnal",
                     "HSM88S", "HSMTLP.50S", "HSMTLP.88S")
select.traits.2 <- c("depth", "Xylem_sap_deltaD_permil", "lwp.min_Predawn", "lwp.min_Diurnal",
                     "TLP",  "p88S", "HSMTLP", "HSM88S", "HSMTLP.88S")

load(file = file.path(results.folder, "hyd.long.prepped.Rdata"))
hyd.pairs.1 <- hyd.long %>%
  subset(trait %in% select.traits.1) %>%
  select(sp, deciduousness, trait.plot.chart, value) %>%
  pivot_wider(names_from = trait.plot.chart, values_from = value)

hyd.pairs.2 <- hyd.long %>%
  subset(trait %in% select.traits.2) %>%
  select(sp, deciduousness, trait.plot.chart, value) %>%
  pivot_wider(names_from = trait.plot.chart, values_from = value)

depth.traits.hyd <- hyd.long %>%
  subset(trait %in% c(select.traits.1, select.traits.2)) %>%
  subset(trait != "depth") %>% left_join(hyd.pairs.2 %>% select(sp, `Depth[italic("Rsq")]`), by = "sp")

select.traits.3 <- c("KmaxL", "lwp.min_Predawn",
                     "lwp.min_Diurnal", "TLP", "p50L",
                     "HSMLWP.TLP", "HSMLWP.50L", "HSMTLP.50L")

select.traits.4 <- c("depth", "Xylem_sap_deltaD_permil",
                     "KmaxL", "lwp.min_Predawn", "TLP", "p50L",
                     "HSMLWP.TLP", "Panama.moist.pref", "Plot.swp.pref", "SG100C_AVG", "Chl")

traits.pairs.1 <- traits.long %>%
  subset(trait %in% select.traits.3) %>%
  select(sp, deciduousness, trait.plot.chart, value) %>%
  pivot_wider(names_from = trait.plot.chart, values_from = value)
traits.pairs.2 <- traits.long %>%
  subset(trait %in% select.traits.4) %>%
  select(sp, deciduousness, trait.plot.chart, value) %>%
  pivot_wider(names_from = trait.plot.chart, values_from = value)

depth.traits.kunert <- traits.long %>%
  subset(trait %in% c(select.traits.3, select.traits.4)) %>%
  subset(trait != "depth") %>% left_join(traits.pairs.2 %>% select(sp, `Depth[italic("Rsq")]`), by = "sp")
save(depth.traits.kunert, file = file.path(results.folder, "depth.traits.kunert.Rdata"))

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


chart.hyd.1 <- ggpairs(hyd.pairs.1 %>% select(-sp, -deciduousness),
                       upper = list(continuous = wrap(cor_func,
                                                      method = 'spearman', symbol = expression('\u03C1 ='))),
                       lower = list(continuous = function(data, mapping, ...) {
                         ggally_smooth_lm(data = data, mapping = mapping) +
                           theme(panel.background = element_blank())}),
                       diag = list(continuous = function(data, mapping, ...) {
                         ggally_densityDiag(data = data, mapping = mapping) +
                           theme(panel.background = element_blank())}
                       ), labeller = "label_parsed")
chart.hyd.2 <- ggpairs(hyd.pairs.2 %>% select(-sp, -deciduousness),
                       upper = list(continuous = wrap(cor_func,
                                                      method = 'spearman', symbol = expression('\u03C1 ='))),
                       lower = list(continuous = function(data, mapping, ...) {
                         ggally_smooth_lm(data = data, mapping = mapping) +
                           theme(panel.background = element_blank())}),
                       diag = list(continuous = function(data, mapping, ...) {
                         ggally_densityDiag(data = data, mapping = mapping) +
                           theme(panel.background = element_blank())}
                       ), labeller = "label_parsed")
# chart.hyd.2 <- chart.hyd.1 %+% hyd.pairs.2
# chart.traits.1 <- chart.hyd.1 %+% traits.pairs.1
# chart.traits.2 <- chart.hyd.1 %+% traits.pairs.2
chart.traits.1 <- ggpairs(traits.pairs.1 %>% select(-sp, -deciduousness),
                          upper = list(continuous = wrap(cor_func,
                                                         method = 'spearman', symbol = expression('\u03C1 ='))),
                          lower = list(continuous = function(data, mapping, ...) {
                            ggally_smooth_lm(data = data, mapping = mapping) +
                              theme(panel.background = element_blank())}),
                          diag = list(continuous = function(data, mapping, ...) {
                            ggally_densityDiag(data = data, mapping = mapping) +
                              theme(panel.background = element_blank())}
                          ), labeller = "label_parsed")
chart.traits.2 <- ggpairs(traits.pairs.2 %>% select(-sp, -deciduousness),
                          upper = list(continuous = wrap(cor_func,
                                                         method = 'spearman', symbol = expression('\u03C1 ='))),
                          lower = list(continuous = function(data, mapping, ...) {
                            ggally_smooth_lm(data = data, mapping = mapping) +
                              theme(panel.background = element_blank())}),
                          diag = list(continuous = function(data, mapping, ...) {
                            ggally_densityDiag(data = data, mapping = mapping) +
                              theme(panel.background = element_blank())}
                          ), labeller = "label_parsed")

ggsave(file.path(figures.folder, paste0("Brett_Wolfe_traits_cor.chart.jpeg")),
       plot = chart.hyd.1 + ggpairs.theme, height = 10, width = 10, units ='in')
ggsave(file.path(figures.folder, paste0("Brett_Wolfe_traits_depth_isotopes_cor.chart.jpeg")),
       plot = chart.hyd.2 + ggpairs.theme, height = 10, width = 10, units ='in')
ggsave(file.path(figures.folder, paste0("Kunert_traits_cor.chart.jpeg")),
       plot = chart.traits.1 + ggpairs.theme, height = 8, width = 8, units ='in')
ggsave(file.path(figures.folder, paste0("Kunert_traits_depth_isotopes_cor.chart.jpeg")),
       plot = chart.traits.2 + ggpairs.theme, height = 9, width = 10, units ='in')


depth.traits.hyd.plot <- ggplot(depth.traits.hyd,
                                aes(y = `Depth[italic("Rsq")]`, x = value)) +
  geom_smooth(method = "lm") +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 0.8, size = 2.5) +
  scale_y_reverse() +
  coord_cartesian(ylim = c(10, 0)) +
  ylab("Effective Rooting Depth (m)") + xlab("") +
  facet_wrap(. ~ trait.plot, scales = "free_x", labeller = label_parsed) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.10, npcy = 0.25, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.10, npcy = 0.1, size = 4) +
  theme(panel.spacing = unit(1, "lines"))
ggsave(file.path(figures.folder, paste0("Wolfe_traits_depth.jpeg")),
       plot = depth.traits.hyd.plot, height = 6, width = 6, units ='in')



depth.traits.kunert.plot <- ggplot(depth.traits.kunert %>%
                                     subset(trait != "Xylem_sap_deltaD_permil"),
                                   aes(y = `Depth[italic("Rsq")]`, x = value)) +
  geom_smooth(method = "lm") +
  geom_point(shape = 21, color = "white", fill = "black", alpha = 0.8, size = 2.5) +
  scale_y_reverse() +
  coord_cartesian(ylim = c(10, 0)) +
  ylab("Effective Rooting Depth (m)") + xlab("") +
  facet_wrap(. ~ trait.plot, scales = "free_x", labeller = label_parsed) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.10, npcy = 0.25, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.10, npcy = 0.1, size = 4) +
  theme(panel.spacing = unit(1.5, "lines"))
ggsave(file.path(figures.folder, paste0("Kunert_traits_depth.jpeg")),
       plot = depth.traits.kunert.plot, height = 6.5, width = 7, units ='in')


tlp.hyd.kunert <- hyd.pairs.1 %>% select(sp, `Psi[tlp]`) %>%
  rename(`Wolfe et al. Psi[tlp]` = `Psi[tlp]`) %>%
  left_join(traits.pairs.2 %>% select(sp, `Psi[tlp]`) %>%
              rename(`Kunert et al. Psi[tlp]` = `Psi[tlp]`), by = "sp")
tlp.min.plot = min(c(tlp.hyd.kunert$`Kunert et al. Psi[tlp]`, tlp.hyd.kunert$`Wolfe et al. Psi[tlp]`), na.rm = TRUE)
tlp.max.plot = max(c(tlp.hyd.kunert$`Kunert et al. Psi[tlp]`, tlp.hyd.kunert$`Wolfe et al. Psi[tlp]`), na.rm = TRUE)
tlp.hyd.kunert.plot <- ggplot(tlp.hyd.kunert, aes(x = `Wolfe et al. Psi[tlp]`, y = `Kunert et al. Psi[tlp]`)) +
  geom_point() +
  geom_smooth(method = "lm", formula = formula) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.05, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 5) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.05, npcy = 0.82, size = 5) +
  ylim(c(tlp.min.plot, tlp.max.plot)) + xlim(c(tlp.min.plot, tlp.max.plot))+
  geom_text(data = data.frame(x = -1.2, y = -1), aes(x = x, y = y, label = "1:1"))
ggsave(file.path(figures.folder, paste0("TLP_Kunert_by_Wolfe.jpeg")),
       plot = tlp.hyd.kunert.plot, height = 3, width = 3, units ='in')



#******************************************************
### Plot LWP -----
#******************************************************

ggplot(lwp.all, aes(x = LWP_coll_time, y = lwp.min)) +
  facet_grid(location ~ .) +
  geom_point(aes(color = as.factor(date))) +
  geom_line(aes(group = c(sp_date), color = as.factor(date)), show.legend = FALSE) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, formula = formula) +
  guides(color = guide_legend(title = "Date"))
ggsave(file.path(figures.folder, paste0("LWP_time_series_all_Data.jpeg")), height = 6, width = 5, units ='in')

# Across all the dryseason measurements which diurnal measurement was minimum by species and location
## also get the predawn for the same day
# https://stackoverflow.com/questions/46971945/how-can-i-have-a-greek-symbol-within-one-of-my-facet-labels
plot.lwp.base <- ggplot(lwp.min %>% subset(!is.na(deciduousness)), aes(x = time, y = lwp.min)) +
  facet_grid(. ~ location) +
  geom_point(aes(color = deciduousness)) +
  ylab(expression(psi[min])) + xlab("Time") +
  guides(colour = guide_legend(order = 1, title = "Deciduousness")) +
  scale_x_discrete(name = "",
                   breaks = c("Predawn", "Diurnal"),
                   labels = c("Predawn", expression('Diurnal'['min'])))
plot.lwp <- plot.lwp.base +
  geom_line(aes(group = sp, color = deciduousness)) +
  geom_errorbar(aes(ymax = lwp.min + lwp.se, ymin = lwp.min - lwp.se, color = deciduousness), width = 0.05) +
  ggsave(file.path(figures.folder, paste0("LWP_min.jpeg")), plot = plot.lwp, height = 4, width = 7.5, units ='in')

plot.lwp.diff <- plot.lwp.base %+%
  subset(lwp.diff, !is.na(deciduousness)) +
  geom_line(aes(group = sp_date, color = deciduousness), size = 0.2)
ggsave(file.path(figures.folder, paste0("LWP_min_all_days.jpeg")), plot = plot.lwp.diff, height = 4, width = 7.5, units ='in')

plot.lwp.diff.ts.base <- ggplot(lwp.diff %>% subset(!is.na(deciduousness) & !is.na(lwp.diff) & location != "PA-BCI"),
                                aes(x = date, y = lwp.diff)) +
  facet_grid(. ~ location) +
  ylab(expression(psi[Predawn] - psi['Diurnal'['min']])) + xlab("Time")
plot.lwp.diff.ts <- plot.lwp.diff.ts.base +
  geom_point(aes(color = deciduousness)) +
  geom_line(aes(group = sp, color = deciduousness), size = 0.2) +
  guides(colour = guide_legend(order = 1, title = "Deciduousness"))
ggsave(file.path(figures.folder, paste0("LWP_min_all_days.jpeg")), plot = plot.lwp.diff.ts, height = 4, width = 7.5, units ='in')
plot.lwp.diff.ts.sp <- plot.lwp.diff.ts.base +
  geom_point(aes(shape = deciduousness, color = deci_sp), size = 3) +
  geom_line(aes(group = deci_sp, color = deci_sp), size = 0.2) +
  guides(colour = guide_legend(order = 1, title = "Sp"),
         shape = guide_legend(order = 2, title = "Deciduousness")) +
  theme(legend.direction = "vertical", legend.box = "horizontal")
ggsave(file.path(figures.folder, paste0("LWP_min_all_days_sp_color.jpeg")), plot = plot.lwp.diff.ts.sp, height = 4, width = 8, units ='in')

plot.lwp.ts.sp.pnm <- ggplot(lwp.diff %>% subset(!is.na(deciduousness) & location == "PA-PNM") %>% droplevels(),
                             aes(x = time, y = lwp.min)) +
  facet_grid(location ~ as.factor(date)) +
  ylab(expression(psi[min])) + xlab("Time") +
  ylim(c(-3.5, 0)) +
  scale_x_discrete(name="",
                   breaks = c("Predawn", "Diurnal"),
                   labels = c("Predawn", expression('Diurnal'['min']))) +
  geom_point(aes(color = deci_sp, shape = deciduousness), size = 2) +
  geom_line(aes(group = deci_sp, color = deci_sp), size = 0.2) +
  guides(colour = guide_legend(order = 1, title = "Sp"),
         shape = guide_legend(order = 2, title = "Deciduousness")) +
  theme(legend.direction = "vertical", legend.box = "horizontal")
plot.lwp.ts.sp.snl <- plot.lwp.ts.sp.pnm %+%
  droplevels(subset(lwp.diff, !is.na(deciduousness) & location == "PA-SLZ")) +
  scale_color_brewer(palette = "Set1")
ggsave(file.path(figures.folder, paste0("LWP_min_all_days_diurnal_sp_color.jpeg")),
       plot = arrangeGrob(plot.lwp.ts.sp.pnm, plot.lwp.ts.sp.snl), height = 5.5, width = 12, units ='in')


#******************************************************
####----Plot Phenology by Wolfe hydraulic traits-----
#******************************************************

h1 <- ggplot(hyd.labels.data, aes(x = deciduousness, y = value)) +
  geom_boxplot(data = hyd.long, aes(fill = deciduousness), stat = "boxplot", notch = TRUE) +
  geom_jitter(data = hyd.long, width = 0.05, shape = 21, fill = "darkgray", color = "black", show.legend = FALSE, alpha = 0.7) +
  facet_wrap(. ~  trait.plot, scales = "free_y", labeller = label_parsed) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 1, hjust = 1)) +
  scale_color_brewer(palette = "Greens", direction = -1) +
  scale_fill_brewer(name = "Deciduousness", palette = "Greens", direction = -1) +
  xlab("Deciduousness") + ylab("Value")
ggsave(file.path(figures.folder, paste0("BrettWolfe_traits_vs_deciduousness_isotopes_depth.jpeg")),
       plot = h1, height = 7, width = 10, units ='in')
h1.1 <- h1 + geom_text(aes(label = groups), vjust = 1, hjust = 0, show.legend = FALSE)
ggsave(file.path(figures.folder, paste0("BrettWolfe_traits_vs_deciduousness_kruskal.labels.jpeg")),
       plot = h1.1, height = 7, width = 10, units ='in')

h1.2 <- ggplot(hyd.labels.data %>% subset(!trait %in% c("depth",  "Xylem_sap_deltaD_permil")),
               aes(x = deciduousness, y = value)) +
  geom_boxplot(data = hyd.long %>% subset(!trait %in% c("depth",  "Xylem_sap_deltaD_permil")),
               aes(fill = deciduousness), stat = "boxplot", notch = TRUE) +
  geom_jitter(data = hyd.long %>% subset(!trait %in% c("depth",  "Xylem_sap_deltaD_permil")),
              width = 0.05, shape = 21, fill = "darkgray", color = "black", show.legend = FALSE, alpha = 0.7) +
  facet_wrap(. ~  trait.plot, scales = "free_y", labeller = label_parsed) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 1, hjust = 1)) +
  scale_color_brewer(palette = "Greens", direction = -1) +
  scale_fill_brewer(name = "Deciduousness", palette = "Greens", direction = -1) +
  xlab("Deciduousness") + ylab("Value")
ggsave(file.path(figures.folder, paste0("BrettWolfe_traits_vs_deciduousness.jpeg")),
       plot = h1.2, height = 7, width = 10, units ='in')
h1.3 <- h1.2 + geom_text(aes(label = groups), vjust = 1, hjust = 0, show.legend = FALSE)
ggsave(file.path(figures.folder, paste0("BrettWolfe_traits_vs_deciduousness_kruskal.labels.jpeg")),
       plot = h1.3, height = 7, width = 10, units ='in')

select.traits <- c("depth", "Xylem_sap_deltaD_permil", "lwp.min_Predawn", "lwp.min_Diurnal", "TLP", "p50S", "p88S",
                   "HSMTLP", "HSM50S","HSM88S", "HSMTLP.50S", "HSMTLP.88S")
hyd.long <- hyd.long %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(deciduousness)]), ordered=TRUE),
         deci_sp.plot = factor(deci_sp, levels=unique(deci_sp[order(deciduousness)]), ordered=TRUE))
h2 <- ggplot(hyd.long %>% subset(trait %in% select.traits)) +
  facet_wrap(. ~  trait.plot, scales = "free_x", labeller = label_parsed, nrow = 2) +
  geom_col(aes(x = deci_sp.plot, y = value,
               fill = deciduousness),
           position = position_dodge2(width = 0.9, preserve = "single")) +
  theme(axis.text.y = element_text(face = "plain", size = 8)) +
  coord_flip() +
  guides(fill = guide_legend(title = "Deciduousness")) +
  xlab("Deciduousness") + ylab("Value")
ggsave(file.path(figures.folder, paste0("BrettWolfe_traits_vs_deciduousness_sp_bar_HSM.jpeg")),
       plot = h2, height = 7, width = 10, units ='in')
h3 <- h2 %+% subset(hyd.long, !trait %in% select.traits)
ggsave(file.path(figures.folder, paste0("BrettWolfe_traits_vs_deciduousness_sp_bar_CWR.jpeg")),
       plot = h3, height = 7, width = 12, units ='in')

#******************************************************
####----Plot Phenology by Kunert hydraulic traits-----
#******************************************************

ggplot(traits.labels.data %>% subset(trait != "se"), aes(x = deciduousness, y = value)) +
  facet_wrap(. ~  trait.plot, scales = "free_y", labeller = label_parsed) +
  geom_text(aes(label = groups), vjust = 1, hjust = 0, show.legend = FALSE) +
  geom_boxplot(data = traits.long %>% subset(trait != "se"), aes(fill = deciduousness), stat = "boxplot", notch = TRUE) +
  geom_jitter(data = traits.long %>% subset(trait != "se"), width = 0.05, shape = 21, fill = "darkgray",
              color = "black", show.legend = FALSE, alpha = 0.7) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 1, hjust = 1)) +
  scale_color_brewer(palette = "Greens", direction = -1) +
  scale_fill_brewer(name = "Deciduousness", palette = "Greens", direction = -1)
ggsave(file.path(figures.folder, paste0("kunert_traits_vs_deciduousness.jpeg")), height = 8, width = 11, units ='in')

## species wise for sp in hyd.traits

t2 <- ggplot(traits.long.hyd %>% subset(sp %in% unique(hyd$sp) & trait != "se")) +
  facet_wrap(. ~  trait.plot, scales = "free_x", labeller = label_parsed, nrow = 2) +
  geom_col(aes(x = deci_sp.plot, y = value,
               fill = deciduousness),
           position = position_dodge2(width = 0.9, preserve = "single")) +
  theme(axis.text.y = element_text(face = "plain", size = 8)) +
  coord_flip() +
  guides(fill = guide_legend(title = "Deciduousness")) +
  xlab("Deciduousness") + ylab("Value")
ggsave(file.path(figures.folder, paste0("Kunert_traits_vs_deciduousness_sp_bar.jpeg")),
       plot = t2, height = 7, width = 14, units ='in')

## Thus evergreen species have:
# significantly higher median TLP than F. Deci,
# significantly lower Leaf Kmax than F Deci and Obligate deci
# significantly higher SPAD than F. Deci or O. Deci
# higher lwp.min, greater moist site preference
## Evg, BDeci & FDeci have significantly higher Chl and WD than O Deci

#******************************************************
## Mortality vs growth rate by phenology-----
#******************************************************
n.threshold = 50
growth.type <- "med"
formula = y ~ x
demo.sp_size <- demo.sp_size %>%
  left_join(deci, by = "sp")
p.d0 <- ggplot(demo.sp_size %>% subset(size != "NA"),
               aes(y = mrate, x = grate)) +
  geom_point() +
  facet_wrap(. ~ size, scales = "free_y") +
  geom_smooth(method = "lm", formula = formula) +
  scale_x_log10() + scale_y_log10() +
  xlab(expression("Mean Growth Rate (mm/yr)")) +  ylab(expression("Mean Mortality Rate (% per year)"))
p.d0.1 <- p.d0 +  stat_poly_eq(npcx = 0.1, npcy = 0.2, size = 4, aes(label = paste(..rr.label..)), rr.digits = 2, formula = formula, parse = TRUE) +
  stat_fit_glance(npcx = 0.1, npcy = 0.1, size = 4, method = 'lm',  method.args = list(formula = formula), geom = 'text_npc', aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")))
ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_mean_Growth_vs_mrate.jpeg")),
       plot = p.d0.1, height = 6, width = 6, units='in')

p.d1 <-  p.d0 + facet_wrap(size ~ deciduousness, scales = "free_y") + theme(axis.text.x = element_text(angle = 90)) +
  stat_poly_eq(npcx = 0.05, npcy = 0.2, size = 4, aes(label = paste(..rr.label..)), rr.digits = 2, formula = formula, parse = TRUE) +
  stat_fit_glance( npcx = 0.05, npcy = 0.1, size = 4, method = 'lm',  method.args = list(formula = formula), geom = 'text_npc', aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")))
ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_mean_Growth_vs_mrate_Deci.jpeg")),
       plot = p.d1, height = 9, width = 10, units='in')
p.d2 <-  p.d0 %+% subset(demo.sp_size, avg.abund >= n.threshold  &
                           !deciduousness %in% c("Obligate Deciduous")) +
  facet_wrap(. ~ deciduousness, scales = "free_y") + theme(axis.text.x = element_text(angle = 90)) +
  stat_poly_eq(npcx = 0.85, npcy = 0.2, size = 4, aes(label = paste(..rr.label..)), rr.digits = 2, formula = formula, parse = TRUE) +
  stat_fit_glance( npcx = 0.85, npcy = 0.1, size = 4, method = 'lm',  method.args = list(formula = formula), geom = 'text_npc', aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")))+
  ggtitle("Large Size class (>= 30 cm DBH)") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(paste0("figures/mortality/", growth.type, "/sp_mean_Growth_vs_mrate_Deci_avg.abund_above", n.threshold,"_large.jpeg")),
       plot = p.d2, height = 6, width = 6, units='in')

#******************************************************
### Is leaf phenology linked to vulnerability to different drought intensity and duration?------
#******************************************************
y.label.1 <- expression(atop(Mortality[Interval], '-'~Mean[Mortality]~('%'*yr^{-1})))

m1 <- ggplot(mrate.long %>%
               subset(!is.na(size) & avg.abund >= n.threshold & !is.na(deciduousness)),
             aes(x = deciduous, y = diff.mrate, color = avg.abund)) +
  scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
                       low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
  facet_grid(size.num ~ censusint.m, scales = "free_y") +
  geom_hline(aes(yintercept = 0), color = "blue", size = 0.5) +
  geom_boxplot(aes(fill = deciduousness.label), stat = "boxplot", notch = TRUE) +
  scale_fill_brewer(name = "", palette = "Greens", direction = -1) +
  theme(legend.position = "top") +
  ylab(y.label.1) + xlab("Deciduousness") +
  ggtitle("Mortality Anomaly by Leaf Phenology") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_deci_by_size.jpeg")), plot = m1, height = 8, width = 9, units='in')
m2 <- m1 %+% subset(mrate.long, size == "large" & avg.abund >= n.threshold & !is.na(deciduousness))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_deci_by_size_large.jpeg")), plot = m2, height = 3, width = 9, units='in')
m2 <- m1 %+% subset(mrate.long, size == "large" & avg.abund >= n.threshold & !is.na(deciduousness)) +
  geom_jitter(width = 0.05, shape = 21, fill = "darkgray", color = "black", show.legend = FALSE, alpha = 0.7)
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_deci_by_size_large_points.jpeg")), plot = m2, height = 3, width = 9, units='in')
m2.1 <- ggplot(mrate.long %>%
                 subset(size == "large" & avg.abund >= n.threshold & !is.na(deciduousness)),
               aes(x = censusint.m, y = diff.mrate, color = avg.abund)) +
  scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
                       low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
  facet_grid(size.num ~ ., scales = "free_y") +
  geom_hline(aes(yintercept = 0), color = "blue", size = 0.5) +
  geom_boxplot(stat = "boxplot", notch = TRUE) +
  scale_fill_brewer(name = "", palette = "Greens", direction = -1) +
  theme(legend.position = "top") +
  ylab(y.label.1) + xlab("Census Interval") +
  ggtitle("Mortality Anomaly by Census Interval") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(width = 0.05, shape = 21, fill = "darkgray", color = "black", show.legend = FALSE, alpha = 0.7)
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_by_size_large_points.jpeg")), plot = m2.1, height = 3, width = 6, units='in')

mrate.long.hyd <- subset(mrate.long, sp %in% hyd$sp) %>%
  left_join(hyd.wide %>% select(sp, p88S, HSMTLP.88S, HSM88S, HSM50S), by = "sp") %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(HSMTLP.88S)]), ordered=TRUE))
# show_col(viridis_pal()(4))
m3.base <- ggplot(subset(mrate.long.hyd, size == "large" & sp %in% c("cordal", "luehse", "tab1ro"))) +
  facet_grid(size.num ~ censusint.m, scales = "free_y") +
  theme(legend.position = "top") +
  ylab(y.label.1) + xlab("Deciduousness") +
  guides(fill = guide_legend(title = "Deciduousness")) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y = element_text(face = "plain", size = 8))
m3.1 <- m3.base + geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Species leafless in early wet season, with increasing HSM '*Psi[tlp]*' - '*Psi['88,stem']))
#  does not work: fill = "Facultative Deciduous"
# guides(fill = guide_legend(title = "Deciduousness",
#                            override.aes = list(fill = c("Facultative Deciduous" = "#35B779FF"))))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_deci_HSMTLP.88S_increasing_spp_leafless in early wet season.jpeg")),
       plot = m3.1, height = 3, width = 9, units='in')

m4.1 <- m3.base %+% subset(mrate.long.hyd, sp %in% hyd$sp & size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Mortality for Evergreen Species with increasing HSM '*Psi[tlp]*' - '*Psi['88,stem']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_HSMTLP.88S_increasing_evergreens.jpeg")),
       plot = m4.1, height = 3, width = 9, units='in')

mrate.long.hyd <- mrate.long.hyd %>% mutate(sp.plot = factor(sp, levels=unique(sp[order(-p88S)]), ordered=TRUE))
m3.2 <- m3.base + geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Species leafless in early wet season, with increasingly more negative '*Psi['88,stem']))
ggsave(file.path(paste0(figures.folder, "/sp_Mortality_rate_by_period_deci_p88S_increasing_spp_leafless in early wet season.jpeg")),
       plot = m3.2, height = 3, width = 9, units='in')

m4.2 <- m3.base %+% subset(mrate.long.hyd, sp %in% hyd$sp & size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Mortality for Evergreen Species with increasingly more negative '*Psi['88,stem']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_p88S_increasing_evergreens.jpeg")),
       plot = m4.2, height = 3, width = 9, units='in')
mrate.long.hyd <- mrate.long.hyd %>% mutate(sp.plot = factor(sp, levels=unique(sp[order(-HSM50S)]), ordered=TRUE))
m4.3 <- m3.base %+% subset(mrate.long.hyd, sp %in% hyd$sp & size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Mortality for Evergreen Species with increasing HSM '*Psi[min]*' - '*Psi['50,stem']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_HSMLWP.50S_increasing_evergreens.jpeg")),
       plot = m4.3, height = 3, width = 9, units='in')

mrate.long.traits <- subset(mrate.long, sp %in% traits$sp) %>%
  left_join(traits.wide %>% select(-deciduous, -deciduousness), by = "sp") %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(KmaxL)]), ordered=TRUE))

m4.4 <- m3.base %+% subset(mrate.long.traits, size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Mortality for Evergreen Species with increasing '*italic('K')['max, Leaf']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_KmaxL_increasing_evergreens.jpeg")),
       plot = m4.4, height = 3, width = 9, units='in')

mrate.long.traits <- mrate.long.traits %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(-HSMLWP.50L)]), ordered=TRUE))
m4.5 <- m3.base %+% subset(mrate.long.traits, size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Mortality for Evergreen Species with increasing HSM '*Psi[min]*' - '*Psi['50, Leaf']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_HSMLWP.50L_increasing_evergreens.jpeg")),
       plot = m4.5, height = 3, width = 9, units='in')

## all evergreens
mrate.long <- mrate.long %>% mutate(sp.plot = factor(sp, levels=unique(sp[order(mean.mrate)]), ordered=TRUE))
m5 <- m3.base %+% subset(mrate.long, size == "large" & deciduous == "E") +
  geom_col(aes(x = sp, y = diff.mrate, fill = deciduousness)) +
  ggtitle("Mortality for Evergreen Species") +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 0.5, hjust = 1, size = 4))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_evergreens.jpeg")),
       plot = m5, height = 3, width = 15, units='in')
mrate.long <- mrate.long %>% mutate(sp.plot = factor(sp, levels=unique(sp[order(mean.mrate)]), ordered=TRUE))

m6 <- m3.base %+% subset(mrate.long, size == "large" & deciduous == "DF") +
  geom_col(aes(x = sp, y = diff.mrate, fill = deciduousness)) +
  ggtitle("Mortality for Facultative Deciduous Species") +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 0.5, hjust = 1, size = 4))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_Facultative Deciduous.jpeg")),
       plot = m6, height = 3, width = 15, units='in')

## Ordered along Rooting Depth Index
mrate.long.depth <- mrate.long %>%
  left_join(subset(depth.rsq.isotopes, corr.func == "gr.Psi.Rad.VPD"), by = "sp") %>%
  left_join(bci.traits %>% select(form1, sp), by = "sp")
mutate(sp.plot = factor(sp, levels = unique(sp[order(depth)]), ordered = TRUE))
m4.6 <- m3.base %+% subset(mrate.long.depth, size == "large" &
                             form1 == "T" &!is.na(deciduousness) &!is.na(depth)) +
  geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  # facet_grid(deciduousness ~ censusint.m) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 0.5, hjust = 1, size = 4)) +
  ggtitle(expression('Mortality for Canopy Species with increasing Rooting Depth Index'))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_best_corr_depth_increasing_canopy_sp.jpeg")),
       plot = m4.6, height = 3, width = 10, units='in')

## For each species get % days spent below the kl80 cutoff in the depth of best-correlation (and the depths below)?
## Test whether that explains mortality in that census (in reality this would be cumulative)


#******************************************************
# Is leaf phenology linked to growth vulnerability to different drought intensity and duration?------
# Indeed facultative deciduous species show greater reduction in growth in 2005-2010 period two successive early wet seasons were dry
# a large chunk of their limited growing period
#******************************************************
y.label.2 <- expression(atop(Std.~Growth[interval], '-'~Mean[Std.~Growth]))

g1 <- ggplot(grate.long %>%
               subset(!is.na(size) & !is.na(deciduousness)),
             aes(x = deciduous, y = median)) +
  scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
                       low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
  facet_grid(size.num ~ censusint.m, scales = "free_y") +
  geom_hline(aes(yintercept = 0), color = "blue", size = 0.5) +
  geom_boxplot(aes(fill = deciduousness.label), stat = "boxplot", notch = TRUE) +
  scale_fill_brewer(name = "", palette = "Greens", direction = -1) +
  theme(legend.position = "top") +
  ylab(y.label.2) + xlab("Deciduousness") +
  ggtitle("Growth Rates by Leaf Phenology") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_deci_by_size.jpeg")), plot = g1, height = 8, width = 9, units='in')
g2 <- g1 %+% subset(grate.long, size == "large" & !is.na(deciduousness))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_deci_by_size_large.jpeg")), plot = g2, height = 3, width = 9, units='in')
g2 <- g1 %+% subset(grate.long, size == "large" & !is.na(deciduousness)) +
  geom_jitter(width = 0.05, shape = 21, fill = "darkgray", color = "black", show.legend = FALSE, alpha = 0.7)
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_deci_by_size_large_points.jpeg")), plot = g2, height = 3, width = 9, units='in')

g2.1 <- ggplot(grate.long %>%
                 subset(size == "large" & !is.na(deciduousness)),
               aes(x = censusint.m, y = median)) +
  scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
                       low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
  facet_grid(size.num ~ ., scales = "free_y") +
  geom_hline(aes(yintercept = 0), color = "blue", size = 0.5) +
  geom_boxplot(stat = "boxplot", notch = TRUE) +
  scale_fill_brewer(name = "", palette = "Greens", direction = -1) +
  theme(legend.position = "top") +
  ylab(y.label.2) + xlab("Census Interval") +
  ggtitle("Growth Anomaly by Census Interval") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(width = 0.05, shape = 21, fill = "darkgray", color = "black", show.legend = FALSE, alpha = 0.7)
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_by_size_large_points.jpeg")), plot = g2.1, height = 3, width = 6, units='in')

hyd.wide <- hyd.long %>% pivot_wider(names_from = trait, values_from = value, -trait.plot)
grate.long.hyd <- subset(grate.long, sp %in% hyd$sp) %>%
  left_join(hyd.wide %>% select(sp, p88S, HSMTLP.88S, HSM88S), by = "sp") %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(HSMTLP.88S)]), ordered=TRUE))

# show_col(viridis_pal()(4))
g3.base <- ggplot(subset(grate.long.hyd, size == "large" & sp %in% c("cordal", "luehse", "tab1ro"))) +
  facet_grid(size.num ~ censusint.m, scales = "free_y") +
  theme(legend.position = "top") +
  ylab(expression("Growth Rate - Mean (% per year)")) + xlab("Deciduousness") +
  # ylab(expression('Growth Rate'[interval]*' - '*italic('E')'[Growth Rate'[interval]'] (% yr'^-1')')) +
  # xlab("Deciduousness") +
  guides(fill = guide_legend(title = "Deciduousness")) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 0.5, hjust = 1))
g3.1 <- g3.base + geom_col(aes(x = sp.plot, y = median, fill = deciduousness)) +
  ggtitle(expression('Species leafless in early wet season, with increasing HSM '*Psi[tlp]*' - '*Psi['88,stem']))
#  does not work: fill = "Facultative Deciduous"
# guides(fill = guide_legend(title = "Deciduousness",
#                            override.aes = list(fill = c("Facultative Deciduous" = "#35B779FF"))))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_deci_HSMTLP.88S_increasing_spp_leafless in early wet season.jpeg")),
       plot = g3.1, height = 4, width = 9, units='in')

g4.1 <- g3.base %+% subset(grate.long.hyd, sp %in% hyd$sp & size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = median, fill = deciduousness)) +
  ggtitle(expression('Growth for Evergreen Species with increasing HSM '*Psi[tlp]*' - '*Psi['88,stem']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_HSMTLP.88S_increasing_evergreens.jpeg")),
       plot = g4.1, height = 4, width = 9, units='in')

grate.long.hyd <- grate.long.hyd %>% mutate(sp.plot = factor(sp, levels=unique(sp[order(-p88S)]), ordered=TRUE))
g3.2 <- g3.base + geom_col(aes(x = sp.plot, y = median, fill = deciduousness)) +
  ggtitle(expression('Species leafless in early wet season, with increasingly more negative '*Psi['88,stem']))
ggsave(file.path(paste0(figures.folder, "/sp_Growth_rate_by_period_deci_p88S_increasing_spp_leafless in early wet season.jpeg")),
       plot = g3.2, height = 4, width = 9, units='in')

g4.2 <- g3.base %+% subset(grate.long.hyd, sp %in% hyd$sp & size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = median, fill = deciduousness)) +
  ggtitle(expression('Growth for Evergreen Species with increasingly more negative '*Psi['88,stem']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_p88S_increasing_evergreens.jpeg")),
       plot = g4.2, height = 4, width = 9, units='in')

## all evergreens
grate.long <- grate.long %>%
  left_join(traits.wide %>% select(sp, TLP, KmaxL), by = "sp") %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(KmaxL)]), ordered=TRUE))

g5 <- g3.base %+% subset(grate.long, size == "large" & deciduous == "E") +
  geom_col(aes(x = sp, y = median, fill = deciduousness)) +
  ggtitle("Growth rates (DBH residuals) for Evergreen Species") +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 0.5, hjust = 1, size = 8))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_evergreens.jpeg")),
       plot = g5, height = 4, width = 15, units='in')

g6 <- g3.base %+% subset(grate.long, size == "large" & deciduous == "DF") +
  geom_col(aes(x = sp, y = median, fill = deciduousness)) +
  ggtitle("Growth rates (DBH residuals) for Facultative Deciduous Species") +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 0.5, hjust = 1, size = 12))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Growth_rate_by_period_facultative_deciduous.jpeg")),
       plot = g6, height = 4, width = 15, units='in')

## Plot mortality by time spent below a threshold in the preferred depth-------

# mrate.depth <- mrate.depth %>%
# left_join(subset(ml.rsq.combine.best, corr.func == "mr.Psi.VPD.I") %>%
#             rename(rdi.mr = depth) %>%
#             dplyr::select(sp, size, rdi.mr), by = c("sp", "size")) %>%

mrate.mfac.depth.to.rdi.gr <- mrate.mfac.depth %>%
  group_by(sp, size) %>%
  subset(!depth > rdi.gr) %>%
  ungroup(sp, size) %>%
  group_by(sp, size, censusint.m) %>%
  mutate(mfac.soil.column.gr = sum(mfac, na.rm = TRUE)) %>%
  ungroup(sp, size, censusint.m)
mrate.mfac.depth.to.rdi.gr.total.int <- mrate.mfac.depth.to.rdi.gr %>%
  dplyr::select(sp, size, mfac.soil.column.gr, censusint.m, mrate, diff.mrate, depth) %>%
  group_by(sp, size, censusint.m) %>%
  summarise(mfac.soil.column.gr = mean(mfac.soil.column.gr, na.rm = TRUE),
            mrate = mean(mrate, na.rm = TRUE),
            diff.mrate = mean(diff.mrate, na.rm = TRUE),
            depth = mean(depth, na.rm = TRUE), .groups = "drop_last")
mrate.mfac.depth.to.rdi.gr.total <- mrate.mfac.depth.to.rdi.gr.total.int %>%
  group_by(sp, size) %>%
  summarise(mfac.soil.column.total.gr = sum(mfac.soil.column.gr, na.rm = TRUE), .groups = "drop_last")
mrate.mfac.depth.to.rdi.gr.study <- mrate.mfac.depth.to.rdi.gr %>%
  subset(depth == rdi.gr) %>%
  group_by(sp, size) %>%
  summarise(mfac.total.gr = sum(mfac, na.rm = TRUE),
            mrate.sum = sum(mrate, na.rm = TRUE), .groups = "drop_last")
## mfac vs. mort

# mrate.mfac.depth.gr.mean.mfac.mrate <- mrate.mfac.depth.gr.mean.mfac %>%
#   subset(deciduous == "E") %>%
#   left_join(mrate.depth.mean %>% select(sp, mrate, mrate.se, grate, grate.se), by = "sp")
# formula.2 = y ~ splines::bs(x, 2)
# mfac.plot.mrate <- ggplot(mrate.mfac.depth.gr.mean.mfac.mrate %>% subset(deciduous == "E"),
#                           aes(x = mfac, y = mrate)) +
#   # geom_smooth(method = lm, formula = formula.2, se = FALSE) +
#   geom_errorbar(aes(ymin = mrate - mrate.se, ymax = mrate + mrate.se), width = 0.15, size = 0.1) +
#   geom_point(shape = 21, color = "white", fill = "black", alpha = 1, size = 2.5) +
#   ylab(expression('Mean Mortality Rate (%'*'year'^1*')')) +
#   stat_poly_eq(aes(label = stat(eq.label)),
#                npcx = 0.95, npcy = 0.95, rr.digits = 2,
#                formula = formula.2, parse = TRUE, size = 4) +
#   stat_poly_eq(aes(label = paste(..rr.label..)),
#                npcx = 0.95, npcy = 0.85, rr.digits = 2,
#                formula = formula.2, parse = TRUE, size = 4) +
#   stat_fit_glance(method = 'lm',
#                   method.args = list(formula = formula.2),
#                   geom = 'text_npc',
#                   aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
#                   npcx = 0.95, npcy = 0.75, size = 4) +
#   xlab(expression(atop('Time spent below '*Psi['crit'], '(Days over 1990-2015)')))
# ggsave(file.path(paste0(figures.folder, "/mean_mfac vs. mrate.tiff")),
#        plot = mfac.plot.mrate, height = 3.5, width = 3.5, units = 'in')
# ggsave(file.path(paste0(figures.folder, "/mean_mfac vs. mrate.jpeg")),
#        plot = mfac.plot.mrate, height = 3.5, width = 3.5, units = 'in')
#
# mfac.plot.15.2 <- ggplot(mrate.mfac.depth.select %>% subset(deciduous == "E"),
#                          aes(y = mrate, x = mfac)) +
#   geom_smooth(method = lm, formula = formula.2, se = FALSE) +
#   geom_point(shape = 21, color = "white", fill = "black", alpha = 1, size = 2.5) +
#   ylab(expression('Mean Mortality Rate (%'*'year'^1*')')) +
#   facet_grid(. ~ censusint.m) +
#   stat_poly_eq(aes(label = paste(..rr.label..)),
#                npcx = 0.95, npcy = 0.95, rr.digits = 2,
#                formula = formula.2, parse = TRUE, size = 4) +
#   stat_fit_glance(method = 'lm',
#                   method.args = list(formula = formula.2),
#                   geom = 'text_npc',
#                   aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
#                   npcx = 0.95, npcy = 0.8, size = 4) +
#   xlab(expression(atop('Time spent below '*Psi['crit'], '(Days over 1990-2015)')))
#
# ggsave(file.path(paste0(figures.folder, "/mortality_by mfac.tiff")),
#        plot = mfac.plot.15.2, height = 2.5, width = 10, units = 'in')
# ggsave(file.path(paste0(figures.folder, "/mortality_by mfac.jpeg")),
#        plot = mfac.plot.15.2, height = 2.5, width = 10, units = 'in')

## for rdi.mr
mrate.mfac.depth.to.rdi.mr <- mrate.mfac.depth %>%
  group_by(sp, size) %>%
  subset(!depth > rdi.mr) %>%
  ungroup(sp, size) %>%
  group_by(sp, size, censusint.m) %>%
  mutate(mfac.soil.column.mr = sum(mfac, na.rm = TRUE)) %>%
  ungroup(sp, size, censusint.m)
mrate.mfac.depth.to.rdi.mr.total.int <- mrate.mfac.depth.to.rdi.mr %>%
  dplyr::select(sp, size, mfac.soil.column.mr, censusint.m, mrate, diff.mrate, depth) %>%
  group_by(sp, size, censusint.m) %>%
  summarise(mfac.soil.column.mr = mean(mfac.soil.column.mr, na.rm = TRUE),
            mrate = mean(mrate, na.rm = TRUE),
            diff.mrate = mean(diff.mrate, na.rm = TRUE),
            depth = mean(depth, na.rm = TRUE), .groups = "drop_last") %>%
  ungroup(sp, size, censusint.m)
mrate.mfac.depth.to.rdi.mr.total <- mrate.mfac.depth.to.rdi.mr.total.int %>%
  group_by(sp, size) %>%
  summarise(mfac.soil.column.total.mr = sum(mfac.soil.column.mr, na.rm = TRUE), .groups = "drop_last") %>%
  ungroup(sp, size)
mrate.mfac.depth.to.rdi.mr.study <- mrate.mfac.depth.to.rdi.mr %>%
  subset(depth == rdi.mr) %>%
  group_by(sp, size) %>%
  summarise(mfac.total.mr = sum(mfac, na.rm = TRUE),
            mrate.sum = sum(mrate, na.rm = TRUE), .groups = "drop_last") %>%
  ungroup(sp, size)

mrate.mfac.column <- mrate.mfac.depth %>%
  group_by(sp, size, censusint.m) %>%
  mutate(mfac.soil.column = sum(mfac, na.rm = TRUE)) %>%
  ungroup(sp, size, censusint.m)
mrate.mfac.column.total.int <- mrate.mfac.column %>%
  dplyr::select(sp, size, mfac.soil.column, censusint.m, mrate, diff.mrate, depth) %>%
  group_by(sp, size, censusint.m) %>%
  summarise(mfac.soil.column = mean(mfac.soil.column, na.rm = TRUE),
            mrate = mean(mrate, na.rm = TRUE),
            diff.mrate = mean(diff.mrate, na.rm = TRUE),
            depth = mean(depth, na.rm = TRUE), .groups = "drop_last") %>%
  ungroup(sp, size, censusint.m)

mfac.plot.7 <- ggplot(mrate.mfac.depth.to.rdi.gr.study,
                      aes(x = mfac.total.gr, y = mrate.sum)) +
  # facet_wrap(censusint.m ~ ., nrow = 1) +
  geom_point() +
  geom_smooth(method = "lm", formula = formula) +
  ylab(expression('Total Mortality Rate (% '*'year'^1*')')) +
  xlab(expression('Days '*Psi['Soil,z = ERD']*'<'*Psi['P80,Leaf'])) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4)
ggsave(file.path(paste0(figures.folder,
                        "/sp_pc_during_the_study_period_total days below kl80_in_z through_rdi.mr.jpeg")),
       plot = mfac.plot.7, height = 3, width = 3, units='in')

## Mean mortality rate vs. RDI
# mean.mrate.rdi <- mrate.mfac.depth %>%
#   subset(depth == rdi) %>%
#   group_by(sp, size, rdi) %>%
#   summarise(mrate.mean = mean(mrate, na.rm = TRUE),
#             mrate.se = sd(mrate, na.rm = TRUE)/sqrt(n()))
y.label.1 <- expression(atop(Mortality[Interval], '-'~Mean[Mortality]~('%'*yr^{-1})))
mfac.plot.8 <- ggplot(mrate.mfac.depth %>% subset(depth == rdi.mr),
                      aes(y = diff.mrate, x = mfac)) +
  geom_point() + ylab(y.label.1) +
  xlab(expression('Days '*Psi['Soil, z = RDI.mr']*' < '*Psi['crit'])) +
  geom_smooth(method = "lm", formula = formula) +
  facet_grid(. ~ censusint.m ) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/diff.mortality_rate_mfac in rdi.mr.jpeg")),
       plot = mfac.plot.8, height = 3, width = 9, units='in')

mfac.plot.9 <- mfac.plot.8 %+% subset(mrate.mfac.depth, depth == rdi.gr) +
  xlab(expression('Days '*Psi['Soil, z = ERD']*' < '*Psi['P80, Leaf'])) +
  geom_point(aes(color = depth)) +
  scale_color_continuous(trans = "reverse", guide = guide_colorbar(title = "Depth\n(m)"))
ggsave(file.path(paste0(figures.folder, "/diff.mortality_rate_mfac in rdi.gr.jpeg")),
       plot = mfac.plot.9, height = 3, width = 9, units='in')

mfac.plot.10 <- ggplot(mrate.mfac.depth.to.rdi.mr.total.int,
                       aes(y = diff.mrate, x = mfac.soil.column.mr)) +
  geom_point(aes(color = depth)) +
  scale_color_continuous(trans = "reverse", guide = guide_colorbar(title = "Depth\n(m)")) +
  ylab(y.label.1) +
  xlab(expression('Days '*Psi['Soil, z <= RDI.mr']*' < '*Psi['P80, Leaf'])) +
  geom_smooth(method = "lm", formula = formula) +
  facet_grid(. ~ censusint.m ) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/diff.mortality_rate_mfac upto rdi.mr.jpeg")),
       plot = mfac.plot.10, height = 3, width = 9, units='in')

mfac.plot.11 <- ggplot(mrate.mfac.depth.to.rdi.gr.total.int,
                       aes(y = diff.mrate, x = mfac.soil.column.gr)) +
  geom_point(aes(color = depth)) +
  scale_color_continuous(trans = "reverse", guide = guide_colorbar(title = "Depth\n(m)")) +
  ylab(y.label.1) +
  xlab(expression('Mean Days '*Psi['Soil, z <= ERD']*' < '*Psi['P80, Leaf'])) +
  geom_smooth(method = "lm", formula = formula) +
  facet_grid(. ~ censusint.m ) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/diff.mortality_rate_mfac upto rdi.gr.jpeg")),
       plot = mfac.plot.11, height = 3, width = 9, units='in')

mfac.plot.12 <- ggplot(mrate.mfac.column.total.int,
                       aes(y = mrate, x = mfac.soil.column)) +
  geom_point() + ylab(expression('Mortality Rate (% '*'year'^1*')')) +
  xlab(expression('Mean Days '*Psi['Soil, all z']*' < '*Psi['P80, Leaf'])) +
  geom_smooth(method = "lm", formula = formula) +
  facet_grid(. ~ censusint.m ) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_mfac all z.jpeg")),
       plot = mfac.plot.12, height = 3, width = 9, units='in')

mfac.plot.14 <- ggplot(mrate.mfac.depth %>% subset(depth == rdi.mr),
                       aes(y = mrate, x = rdi.mr)) +
  geom_point() + ylab(expression('Mortality Rate (% '*'year'^1*')')) +
  xlab("Depth best-correlated with Mortality Rate (m)") +
  geom_smooth(method = "lm", formula = formula) +
  facet_grid(. ~ censusint.m ) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by rdi.mr.jpeg")),
       plot = mfac.plot.14, height = 3, width = 9, units='in')


# mfac.plot.16 <- mfac.plot.15 %+% subset(mrate.mfac.depth, depth == rdi.gr) +
#   facet_grid(. ~ censusint.m) +
#   ylab(expression('Mortality Rate (% '*'year'^1*')'))
# ggsave(file.path(paste0(figures.folder, "/mortality_rate_by rdi.gr_interval.jpeg")),
#        plot = mfac.plot.16, height = 3, width = 9, units='in')

adult.mrate.mean <- adult.mrate.long %>%
  group_by(sp, deciduousness) %>%
  summarize_at(vars(mean.mrate, avg.abund, mean.grate), mean, na.rm = TRUE) %>%
  mutate(mean.mrate = ifelse(!is.finite(mean.mrate),
                             rep(NA, length(mean.mrate)), mean.mrate)) %>%
  mutate(size = "large") %>%
  left_join(subset(depth.rsq.isotopes, corr.func == "gr.Psi.VPD.leaf.add") %>%
              dplyr::select(sp, size, depth), by = c("sp", "size")) %>%
  left_join(bci.traits %>% dplyr::select(form1, sp), by = "sp") %>%
  subset(size == "large" & form1 == "T")

pm.2 <- ggplot(adult.mrate.mean %>%
                 subset(avg.abund >= n.threshold & deciduousness %in%
                          c("Facultative Deciduous", "Evergreen")),
               aes(x = depth, y = mean.mrate)) +
  geom_point() +
  scale_x_continuous(trans="sqrt", breaks = soil.depths[-c(2,3, 4, 6, 7, 9)]) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.9, npcy = 0.95, rr.digits = 2,
               formula = formula, parse = TRUE, size = 5) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                  npcx = 0.9, npcy = 0.8, size = 5) +
  ylab(expression('Mean Mortality Rate (% '*'year'^1*')')) +
  xlab("Effective Rooting Depth (m)") +
  facet_grid(. ~ deciduousness)  +
  geom_smooth(method = "lm", se = TRUE, formula = formula) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  ggtitle("Adult tree mortality (>= 10 cm DBH)") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(paste0(figures.folder, "/adult_Mortality_vs_udi_with_outliers_avg.abund_above",
                        n.threshold, "_sp_deci.jpeg")), plot = pm.2, height = 4, width = 6, units='in')

pg.2.deci <- pg.2 %+% subset(adult.mrate.mean,!is.na(deciduousness)) +
  facet_grid(. ~ deciduousness)
ggsave(file.path(paste0(figures.folder, "/adult_Growth_vs_udi_with_outliers_sp_deci.jpeg")), plot = pg.2.deci, height = 4, width = 8, units='in')

#******************************************************
### Yearly psi dynamics versus climatology-------
#******************************************************

## psi does not fit a normal/exp/gamma distribution # so treating non-parametrically
psi.stat.1 <- psi %>%
  group_by(interval.yrs, date, depth) %>%
  summarise(median = -median(psi, na.rm = TRUE),
            upper.CI = -quantile(psi, probs = 0.975),
            lower.CI = -quantile(psi, probs = 0.025), .groups = "drop_last")
psi.stat.1 <- psi.stat.1 %>%
  group_by(interval.yrs, depth) %>%
  mutate(days = 1:n())
rectangles <- data.frame(
  xmin = as.Date(paste0(c(1990:2018), "-04-01")),
  xmax = as.Date(paste0(c(1990:2018), "-11-01")),
  ymin = 0,
  ymax = 2.5
)
plot.psi.stat.1 <- ggplot(psi.stat.1 %>%
                            subset(depth %in% c(0.21,  0.37,  1.00,  1.70)) %>% droplevels()) +
  geom_rect(data=rectangles, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill='gray80', alpha=0.8) +
  # geom_ribbon(aes(x = date, ymin = lower.CI, ymax = upper.CI), alpha = 0.3) +
  theme(panel.grid.major.y = element_line()) +
  geom_line(aes(x = date, y = median, group = depth, color = depth), size = 0.3) +
  scale_color_continuous(trans = "reverse", guide = guide_colorbar(title = "Depth\n(cm)")) +
  geom_vline(xintercept = census.beg) + # as.Date(paste0(c(1990:2015), "-01-01")) +
  scale_y_reverse(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 10, 15), limits = c(2.5, 0)) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%Y")) +
  coord_cartesian(xlim = c(as.Date("1990-01-01"), as.Date("2018-12-31"))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  # facet_grid(interval.yrs ~ .) +
  ylab("-Soil Water Potential [MPa]") +
  theme(text = element_text(size = 12)) + xlab("")
# ggtitle("PSI:Violins of interval medians of best-fit ensembles")
ggsave("psi_model_daily_bestfit_params.top.few_CI_full.jpeg", plot = plot.psi.stat.1,
       file.path(figures.folder), device = "jpeg", height = 3, width = 20, units='in')

psi.2 <- psi %>%
  mutate(interval.yrs.to.plot = forcats::fct_explicit_na(cut(date, include.lowest = TRUE, breaks = cut.breaks.2,
                                                             labels = cut.labels.2, right = TRUE)))

psi.stat.2 <- psi.2 %>%
  group_by(interval.yrs.to.plot, date, depth) %>%
  summarise(median = median(psi, na.rm = TRUE),
            upper.CI = quantile(psi, probs = 0.975),
            lower.CI = quantile(psi, probs = 0.025), .groups = "drop_last")
psi.stat.2 <- psi.stat.2 %>%
  group_by(interval.yrs.to.plot, depth) %>%
  mutate(days = 1:n())
rectangles.2 <- data.frame(
  xmin = c(0:4)*365 + 120,
  xmax = c(0:4)*365 + 335,
  ymin = 0,
  ymax = -2.5
)

plot.psi.stat.2 <- ggplot(psi.stat.2 %>%
                            subset(depth %in% c(0.1, 0.21,  0.37,  1.00,  1.70)) %>% droplevels()) +
  geom_rect(data=rectangles.2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill='gray80', alpha=0.8) +
  geom_ribbon(aes(x = days, ymin = lower.CI, ymax = upper.CI), alpha = 0.3, fill = "red") +
  geom_line(aes(x = days, y = median, group = depth, color = depth), size = 0.3) +
  scale_color_continuous(trans = "reverse", guide = guide_colorbar(title = "Depth\n(cm)")) +
  # scale_y_continuous(breaks = -c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 10, 15)) +
  scale_x_continuous(breaks = c(1:5)*365, labels = paste0("Yr.", c(1:5))) +
  coord_cartesian(ylim = c(-2.5, 0)) +
  theme(panel.grid.major.y = element_line()) +
  facet_grid(interval.yrs.to.plot ~ .) +
  ylab("Soil Water Potential (MPa)") +
  theme(text = element_text(size = 12)) + xlab("Census Year")
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_censuspanels.jpeg", plot = plot.psi.stat.2,
       file.path(figures.folder), device = "jpeg", height = 5, width = 7, units='in')

## by depth panels
plot.psi <- ggplot(psi %>% subset(date < as.Date("1992-12-31") & depth >= 1), aes(x = date, y = -psi)) +
  scale_y_continuous(trans="rev_sqrt") +
  geom_line(aes(group = par.sam, color = as.factor(par.sam)), show.legend = F, size = 0.2) +
  geom_vline(xintercept = census.beg, color = "gray") +
  facet_grid(depth ~ ., scales = "free_y") +
  ylab("-Soil Water Potential [MPa]") + xlab("Date") +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%Y")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Soil Water Potential for best-fit ensembles")
ggsave("psi_model_daily_all_depths_params.top.few_full_until_1992-12-31.jpeg", plot = plot.psi, path =
         file.path(figures.folder), height = 15, width = 8.94, units='in')

psi.stat.3 <- psi %>%
  mutate(doy = format(date, "%j")) %>%
  group_by(doy, depth) %>%
  summarise(median = median(psi, na.rm = TRUE),
            upper.CI = quantile(psi, probs = 0.975),
            lower.CI = quantile(psi, probs = 0.025), .groups = "drop_last") %>%
  ungroup(doy, depth) %>%
  mutate(doy = as.numeric(doy))

source("code/Utilities/plot.ticks.R")
plot.psi.stat.5.base <- ggplot(psi.stat.5 %>% droplevels()) +
  # geom_rect(data=rectangles.3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
  #           fill='gray80', alpha=0.8) +
  scale_x_continuous(breaks = c(seq(0, 360, by = 60))) +
  coord_cartesian(ylim = c(-2.5, 0)) +
  theme(panel.grid.major.y = element_line()) +
  ylab(expression(Psi[soil]*~"(MPa)")) + xlab("Day of the Year")
plot.psi.stat.5 <- plot.psi.stat.5.base +
  geom_ribbon(aes(x = doy, ymin = q2.5.clim, ymax = q97.5.clim), alpha = 0.3, fill = "grey20") +
  geom_line(aes(x = doy, y = median.clim, group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
  guides(color = guide_legend(title = "Depth(m)", order = 2, override.aes = list(size = 3)))
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_censuspanels_climatology.jpeg",
       plot = plot.ticks(plot.psi.stat.5),
       file.path(figures.folder), device = "jpeg", height = 3, width = 5, units='in')

plot.psi.stat.5.over <- plot.psi.stat.5.base %+%
  subset(psi.stat.5, depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)) +
  geom_line(aes(x = doy, y = median.clim, linetype = "climatology", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
  geom_line(data = subset(psi.stat.4, year == "2016" & depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)),
            aes(x = doy, y = median, linetype = "2016", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
  guides(color = guide_legend(title = "Depth(m)", legend.position = "right", order = 1, override.aes = list(size = 3)),
         linetype = guide_legend(order = 2, title = NULL, legend.position = "top", override.aes =
                                   list(linetype = c("climatology" = "dashed", "2016" = "solid")))) +
  coord_cartesian(ylim = c(-3, 0)) + ggtitle("2016")
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_censuspanels_climatology_over.jpeg",
       plot = plot.ticks(plot.psi.stat.5.over),
       file.path(figures.folder), device = "jpeg", height = 3, width = 5, units='in')

pdf(paste0(figures.folder, "/psi_model_daily_bestfit_params.top.few_CI_full_censuspanels_climatology_over_by_year.pdf"), height = 4, width = 7)
for (i in unique(psi.stat.4$year)) {
  plot.psi.stat.5.yr <- plot.psi.stat.5.base %+% subset(psi.stat.5, depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)) +
    geom_line(aes(x = doy, y = median.clim, linetype = "climatology", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
    geom_line(data = subset(psi.stat.4, year == i & depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)),
              aes(x = doy, y = median, linetype = "Year", group = as.factor(depth), color = as.factor(depth)), size = 0.5) +
    guides(color = guide_legend(title = "Depth(m)", order = 1, override.aes = list(size = 3)),
           linetype = guide_legend(order = 2, title = NULL, override.aes =
                                     list(linetype = c("climatology" = "solid", "Year" = "dashed")))) +
    coord_cartesian(ylim = c(-3, 0), xlim = c(0, 200)) + ggtitle(i)
  print(plot.psi.stat.5.yr)
}
dev.off()



# load(file = file.path(results.folder, "psi.stat.4.Rdata"))
# load(file = file.path(results.folder, "psi.stat.4.select.Rdata"))

pct.drought.days <- psi.stat.4 %>%
  mutate(season = ifelse(doy < 120, "Dry Season", "Wet Season")) %>%
  group_by(depth, interval.yrs, interval.yrs.2, season) %>%
  summarise(pct.days.below.q10 = 100*round(sum(!is.na(below.q10))/n(), 3),
            pct.days.below.q5 = 100*round(sum(!is.na(below.q5))/n(), 3),
            pct.days.below.q2.5 = 100*round(sum(!is.na(below.q2.5))/n(), 3),
            pct.days.above.q5 = 100 - pct.days.below.q5,
            pct.days.above.q2.5 = 100 - pct.days.below.q2.5,
            pct.days.below.q5.0.5 = 100*round(sum(!is.na(below.q5) & below.q5 < -0.5)/n(), 3),
            pct.days.above.q5.0.5 = 100 - pct.days.below.q5.0.5,
            pct.days.below.q2.5.0.5 = 100*round(sum(!is.na(below.q2.5) & below.q2.5 < -0.5)/n(), 3),
            pct.days.above.q2.5.0.5 = 100 - pct.days.below.q2.5.0.5, .groups = "drop_last") %>%
  ungroup(depth, interval.yrs) %>% #subset(depth %in% c(0.06, 0.62, 1))
  mutate(depth.fac = factor(depth, levels = sort(unique(psi.stat.4$depth), decreasing = TRUE))) %>%
  mutate(int.ssn = paste(interval.yrs, season, sep = "_"))

heat.fr.drought.days.base <- ggplot(pct.drought.days %>% subset(interval.yrs.2 != "(Missing)"),
                                    aes(x = interval.yrs.2, y = depth.fac)) +
  facet_wrap(~season, nrow = 2) +
  xlab("Census Interval") + ylab("Depth (m)") +
  scale_fill_viridis_c(expression('% Days '*Psi['Soil, DOY']*' < '*italic(Q)[italic(p)*'=0.025, '*Psi['Soil, DOY']]),
                       direction = -1, option = "plasma") +
  theme(legend.position = "top", legend.direction = "horizontal") +
  geom_tile(aes(fill = pct.days.below.q2.5)) +
  theme(axis.text.x = element_text(size = 10))
ggsave("pct.days.below.q2.5.clim_by depth & intervalyrs&ssn_full.jpeg",
       plot = heat.fr.drought.days.base, file.path(figures.folder), device = "jpeg", height = 6, width = 5, units='in')

heat.fr.drought.days.study <- heat.fr.drought.days.base %+%
  subset(pct.drought.days, interval.yrs != "(Missing)")
ggsave("pct.days.below.q2.5.clim_by depth & intervalyrs&ssn_study_period.jpeg",
       plot = heat.fr.drought.days.study, file.path(figures.folder), device = "jpeg", height = 6.5, width = 5, units='in')

heat.fr.drought.days.0.5 <- heat.fr.drought.days.base +
  geom_tile(aes(fill = pct.days.below.q2.5.0.5))
ggsave("pct.days.below.q2.5.clim.2.5_by depth & intervalyrs&ssn_full.jpeg",
       plot = heat.fr.drought.days.2.5, file.path(figures.folder), device = "jpeg", height = 6.5, width = 5, units='in')

## alpha = 0.05
heat.fr.drought.days.base.q5 <-  heat.fr.drought.days.base +
  geom_tile(aes(fill = pct.days.below.q5)) +
  scale_fill_viridis_c(expression('% Days '*Psi['Soil, DOY']*' < '*italic(Q)[italic(p)*'=0.05, '*Psi['Soil, DOY']]),
                       direction = -1, option = "plasma")
ggsave("pct.days.below.q5.clim_by depth & intervalyrs&ssn_full.jpeg",
       plot = heat.fr.drought.days.base.q5, file.path(figures.folder), device = "jpeg", height = 6, width = 5, units='in')
heat.fr.drought.days.study.q5 <- heat.fr.drought.days.base.q5 %+%
  subset(pct.drought.days, interval.yrs != "(Missing)")
ggsave("pct.days.below.q5.clim_by depth & intervalyrs&ssn_study_period.jpeg",
       plot = heat.fr.drought.days.study.q5, file.path(figures.folder), device = "jpeg", height = 6.5, width = 5, units='in')

heat.fr.drought.days.q5.0.5 <- heat.fr.drought.days.base.q5 +
  geom_tile(aes(fill = pct.days.below.q5.0.5))
ggsave("pct.days.below.q5.clim.2.5_by depth & intervalyrs&ssn_full.jpeg",
       plot = heat.fr.drought.days.q5.0.5, file.path(figures.folder), device = "jpeg", height = 6.5, width = 5, units='in')

rectangles.4 <- data.frame(
  xmin = 120,
  xmax = 335,
  ymin = 0,
  ymax = -3.0
)

## add legends for ribbon areas for drought rarity, lines for climatology,
# Could add a line for sd2
plot.psi.stat.6.interval.base <- plot.psi.stat.5.base %+%
  subset(psi.stat.4, depth %in% c(0.12, 0.62, 1)) +
  facet_wrap(. ~ interval.yrs.2) +
  geom_rect(data=rectangles.4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill='gray90', alpha=0.8) +
  geom_line(aes(x = doy, y = median.clim, group = as.factor(depth), color = as.factor(depth)), size = 0.3, linetype = "solid") +
  geom_ribbon(aes(x = doy, ymin = below.q5, ymax = median.clim, group = as.factor(depth_year),
                  fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE) +
  theme(panel.grid.major.y = element_line(size = 0.1)) +
  scale_color_discrete(name = "Depth (cm)", labels = c("10","60","100")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  coord_cartesian(ylim = c(-2, 0), xlim = c(0, 200)) +
  theme(legend.position = "top")
plot.psi.stat.6.interval.q10 <- plot.psi.stat.6.interval.base  +
  geom_ribbon(aes(x = doy, ymin = below.q10, ymax = median.clim, group = as.factor(depth_year),
                  fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE)
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_q10.clim.jpeg",
       plot = plot.psi.stat.6.interval.q10, file.path(figures.folder), device = "jpeg", height = 3.25, width = 5, units='in')

plot.psi.stat.7.interval.q10 <- plot.psi.stat.6.interval.base %+%
  subset(psi.stat.4, depth %in% c(0.12, 0.62, 1) & interval.yrs != "(Missing)") +
  facet_wrap(. ~ interval.yrs, nrow = 1) +
  geom_ribbon(aes(x = doy, ymin = below.q10, ymax = median.clim, group = as.factor(depth_year),
                  fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE)
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_study_period_q10.jpeg",
       plot = plot.psi.stat.7.interval.q10, file.path(figures.folder), device = "jpeg", height = 2.5, width = 7, units='in')
# q5
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_q5.jpeg",
       plot = plot.psi.stat.6.interval.base, file.path(figures.folder), device = "jpeg", height = 3.25, width = 5, units='in')
plot.psi.stat.7.interval.q5 <- plot.psi.stat.6.interval.base %+%
  subset(psi.stat.4, depth %in% c(0.12, 0.62, 1) & interval.yrs != "(Missing)") +
  facet_wrap(. ~ interval.yrs, nrow = 1)
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_study_period_q5.jpeg",
       plot = plot.psi.stat.7.interval.q5, file.path(figures.folder), device = "jpeg", height = 2.5, width = 6, units='in')
# q2.5
plot.psi.stat.6.interval.q2.5 <- plot.psi.stat.6.interval.base  +
  geom_ribbon(aes(x = doy, ymin = below.q2.5, ymax = median.clim, group = as.factor(depth_year),
                  fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE)
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_q2.5.jpeg",
       plot = plot.psi.stat.6.interval.q2.5, file.path(figures.folder), device = "jpeg", height = 3.25, width = 5, units='in')
plot.psi.stat.7.interval.q2.5 <- plot.psi.stat.6.interval.base %+%
  subset(psi.stat.4, depth %in% c(0.12, 0.62, 1) & interval.yrs != "(Missing)") +
  facet_wrap(. ~ interval.yrs, nrow = 1) +
  geom_ribbon(aes(x = doy, ymin = below.q2.5, ymax = median.clim, group = as.factor(depth_year),
                  fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE)
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_study_period_q2.5.jpeg",
       plot = plot.psi.stat.7.interval.q2.5, file.path(figures.folder), device = "jpeg", height = 2.5, width = 6, units='in')

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

## plot over observed
plot.psi.stat.7.interval.median <- plot.psi.stat.6.interval.base %+%
  subset(psi.stat.4, depth %in% c(0.12, 0.62, 1) & interval.yrs.2 != "2015-2020") +
  facet_wrap(. ~ interval.yrs, nrow = 1) +
  geom_ribbon(aes(x = doy, ymin = median, ymax = median.clim, group = as.factor(depth_year),
                  fill = as.factor(depth)), alpha = 0.7, show.legend = FALSE)
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_interval_panels_climatology_over_study_period_all_observed.jpeg",
       plot = plot.psi.stat.7.interval.median, file.path(figures.folder), device = "jpeg", height = 2.5, width = 7, units='in')

plot.list.grid <- list()
for (i in 1:length(unique(psi.stat.4$year))) {
  year.on = unique(psi.stat.4$year)[i]
  plot.psi.stat.5.yr <-  ggplot(psi.stat.5 %>% subset(depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)) %>%
                                  droplevels()) +
    # geom_rect(data=rectangles.3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
    #           fill='gray80', alpha=0.8) +
    scale_x_continuous(breaks = c(seq(0, 360, by = 60))) +
    coord_cartesian(ylim = c(-2.5, 0)) +
    theme(panel.grid.major.y = element_line()) +
    ylab("Soil Water Potential (MPa)") + xlab("Day of the Year") +
    theme(text = element_text(size = 12)) +
    # geom_ribbon(aes(x = doy, ymin = lower.CI, ymax = upper.CI), alpha = 0.3, fill = "grey20") +
    geom_line(aes(x = doy, y = median.clim, linetype = "climatology", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
    geom_line(data = subset(psi.stat.4, year == year.on & depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)),
              aes(x = doy, y = median, linetype = "Year", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
    guides(color = guide_legend(title = "Depth(m)", order = 1, override.aes = list(size = 3)),
           linetype = guide_legend(order = 2, title = NULL, override.aes =
                                     list(linetype = c("climatology" = "solid", "Year" = "dashed")))) +
    coord_cartesian(ylim = c(-3, 0), xlim = c(0, xlim.in.wet.season)) +
    ggtitle(year.on) + xlab("") + ylab("")
  if(year.on != "2018") {
    plot.psi.stat.5.yr <- plot.psi.stat.5.yr + theme(legend.position = "none")
  }
  # if (as.numeric(year.on)%%5 != 0) {
  #   plot.psi.stat.5.yr <- plot.psi.stat.5.yr + ylab("")
  # }
  # if (as.numeric(year.on)[i] %in% c(2015:2018)) {
  #   plot.psi.stat.5.yr <- plot.psi.stat.5.yr + xlab("")
  # }
  plot.list.grid[[i]] <- plot.psi.stat.5.yr
}

plot.list <- lapply(plot.list.grid, plot.ticks)
# pdf("all.pdf")
# invisible(lapply(plot.list, print))
# dev.off()

ggsave(file.path(figures.folder, "arrange.pdf"), arrangeGrob(grobs = plot.list.grid, ncol = 5), width = 15, height = 10, units = "in")

### Does % of the growing season < p50 or p80 explain inter-census growth or mortality?
### Use Meinzer data for rooting depths, and assign the depth whose water potential to track

head(psi.stat.4)

### Plot pct.drought.days vs. mortality -------
###***************************************
pct.drought.days <- pct.drought.days %>%
  mutate(interval.num = as.numeric(as.character(recode(interval.yrs.2, `1982-1985` = "1",
                                                       `1985-1990` = "2", `1990-1995` = "3",
                                                       `1995-2000` = "4", `2000-2005` = "5",
                                                       `2005-2010` = "6", `2010-2015` = "7", `2015-2018` = "8"))))
pct.drought.days.mean <- pct.drought.days %>%
  dplyr::select(-season, -int.ssn) %>%
  group_by(interval.num, interval.yrs, interval.yrs.2, depth, depth.fac) %>%
  summarise_all(list(~sum(., na.rm = TRUE)), .groups = "drop_last")

mrate.drought <- mrate.long %>% left_join(pct.drought.days.mean, by = "interval.num")  %>%
  left_join(bci.traits %>% dplyr::select(sp, form1), by = "sp")

mrate.drought.ssn <- mrate.long %>% left_join(pct.drought.days, by = "interval.num")  %>%
  left_join(bci.traits %>% dplyr::select(sp, form1), by = "sp")

mrate.drought.plot <- ggplot(mrate.drought %>%
                               subset(size == "large" & form1 == "T" &
                                        depth %in% c(0.12, 0.37, 1.00)),
                             aes(y = mrate, x = pct.days.below.q2.5)) +
  geom_point(aes(color = sp), show.legend = FALSE) + ylab(expression('Mortality Rate (% '*'year'^1*')')) +
  # xlab('% Drought Days (below 2.5% Quantile)') +
  xlab(expression('% Drought Days ('*Psi['z, DOY']*' < '*italic(Q)[italic(p)*'=0.025, '*Psi['Soil, DOY']]*')')) +
  geom_smooth(method = "loess", formula = formula) +
  facet_grid(depth  ~ .) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'loess',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by_pct.drought.days_by_depth.jpeg")),
       plot = mrate.drought.plot, height = 3, width = 9, units='in')

mrate.drought.plot.2 <- ggplot(mrate.drought %>%
                                 subset(size == "large" & form1 == "T" &
                                          depth %in% c(0.12, 0.37, 1.00)),
                               aes(y = mrate, x = censusint.m)) +
  geom_point(aes(color = pct.days.below.q2.5)) +
  ylab(expression('Mortality Rate (% '*'year'^1*')')) + xlab('Census Interval') +
  guides(color = guide_legend(title = '% Drought Days (below 2.5% Quantile)')) +
  scale_color_continuous(trans = "reverse") +
  theme(legend.position = "top") +
  geom_smooth(method = "loess", formula = formula) +
  facet_grid(depth  ~ .) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'loess',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by_interval_colored_by_pct.drought.days_by_depth.jpeg")),
       plot = mrate.drought.plot.2, height = 3, width = 9, units='in')

###
mrate.drought.plot.3 <- ggplot(mrate.drought %>%
                                 subset(size == "large" & form1 == "T" &
                                          depth %in% c(0.12, 0.37, 1.00)),
                               aes(y = mrate, x = pct.days.below.q2.5.0.5)) +
  geom_point(aes(color = sp), show.legend = FALSE) + ylab(expression('Mortality Rate (% '*'year'^1*')')) +
  # xlab('% Drought Days (below 2.5% Quantile)') +
  xlab(expression('% Drought Days ('*Psi['z, DOY']*' < '*italic(Q)[italic(p)*'=0.025, '*Psi['Soil, DOY']]*' & < 0.5 MPa)')) +
  geom_smooth(method = "loess", formula = formula) +
  facet_grid(depth  ~ .) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'loess',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by_pct.drought.days_below0.5_by_depth.jpeg")),
       plot = mrate.drought.plot.3, height = 3, width = 9, units='in')

mrate.drought.plot.4 <- ggplot(mrate.drought %>%
                                 subset(size == "large" & form1 == "T" &
                                          depth %in% c(0.12, 0.37, 1.00)),
                               aes(y = diff.mrate, x = censusint.m)) +
  geom_point(aes(color = pct.days.below.q2.5.0.5)) +
  geom_line(aes(group = sp, color = pct.days.below.q2.5.0.5)) +
  geom_hline(aes(yintercept = 0)) +
  ylab(y.label.1) + xlab('Census Interval') +
  guides(color = guide_legend(title = '% Drought Days (<2.5% Quantile & <0.5MPa)')) +
  scale_color_continuous(trans = "reverse") +
  theme(legend.position = "top") +
  geom_smooth(method = "loess", formula = formula) +
  facet_grid(depth  ~ .) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'loess',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.8, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by_interval_colored_by_pct.drought.days_below0.5_by_depth.jpeg")),
       plot = mrate.drought.plot.4, height = 4, width = 6.5, units='in')

###
mrate.drought.plot.5 <- ggplot(mrate.drought.ssn %>%
                                 subset(size == "large" & form1 == "T" &
                                          depth %in% c(0.12, 0.37, 1.00)),
                               aes(y = mrate, x = pct.days.below.q2.5)) +
  geom_point(aes(color = sp), show.legend = FALSE) +
  ylab(expression('Mortality Rate (% '*'year'^1*')')) +
  # xlab('% Drought Days (below 2.5% Quantile)') +
  xlab(expression('% Drought Days ('*Psi['z, DOY']*' < '*italic(Q)[italic(p)*'=0.025, '*Psi['Soil, DOY']]*')')) +
  geom_smooth(method = "lm", formula = formula) +
  facet_grid(depth  ~ season) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.6, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by_pct.drought.days_by_depth_ssn.jpeg")),
       plot = mrate.drought.plot.5, height = 3, width = 9, units='in')

mrate.drought.plot.6 <- ggplot(mrate.drought.ssn %>%
                                 subset(size == "large" & form1 == "T" &
                                          depth %in% c(0.12, 0.37, 1.00)),
                               aes(y = diff.mrate, x = censusint.m)) +
  geom_point(aes(color = pct.days.below.q2.5)) +
  geom_line(aes(group = sp, color = pct.days.below.q2.5)) +
  geom_hline(aes(yintercept = 0)) +
  ylab(y.label.1) + xlab('Census Interval') +
  guides(color = guide_legend(title = '% Drought Days (below 2.5% Quantile)')) +
  scale_color_continuous(trans = "reverse") +
  theme(legend.position = "top") +
  geom_smooth(method = "loess", formula = formula) +
  facet_grid(depth  ~ season) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.8, npcy = 0.9, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", round(..p.value.., digits = 3), sep = "")),
                  npcx = 0.8, npcy = 0.6, size = 4) #+ scale_y_sqrt()
ggsave(file.path(paste0(figures.folder, "/mortality_rate_by_interval_colored_by_pct.drought.days_by_depth_ssn.jpeg")),
       plot = mrate.drought.plot.6, height = 4, width = 9, units='in')

