##-------------------------
## Date : 22 Jul 2018
## Author : Rutuja
## Title : Uptake depth and hydraulic traits
##-------------------------
## Do species with shallow uptake depth have greater drought tolerance? (more -ve TLP)

rm(list = ls())
gc()
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, ggcorrplot, Hmisc, PerformanceAnalytics,
               corrplot, devtools, multcompView, agricolae, readxl, ggpmisc)
# install_github("vqv/ggbiplot")
  # graphics info
theme_set(theme_bw())
theme_update(text = element_text(size=14),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)
require(scales)
rev_sqrt_trans <- function() {
  scales::trans_new(
    name = "rev_sqrt",
    transform = function(x) -sqrt(abs(x)),
    inverse = function(x) x^2);
}
## adding uptake depth index
load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
load("results/4.1GLUEsetup_part2_BCI.RData") # has n.ensembles and growth and si matrix
load("results/demo.sp_size.RData")
load("results/demo.sp.RData") ## Mortality rate only for those cases for which avg.abundance greater than 10
intervals <- info$intervals
n.ensembles <- growth_by_si.info$n.ensembles
growth.type <- growth_by_si.info$growth.type
growth.selection <- growth_by_si.info$growth.selection
si.type <- growth_by_si.info$si.type
goodness.fit <- 0.3
soil.depths <- unique(info$root.param.long$depth)
##
dbh.residuals <- "on"#growth_by_si.info$dbh.residuals
dryseason <- "on"
root.selection <- "on"
iso.subset <- "off"
##
file.extension.base4 <- paste0(goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", growth.selection,
                               "_", dbh.residuals, "_", intervals, "_dryseason_", dryseason, "_iso.subset_",
                               iso.subset, "_root.selection_", root.selection)
tlplevels <- "sp"
level.folder <- "splevel"
load(file = paste0("results/", level.folder, "/ds.bestfit_cor", file.extension.base4, ".Rdata"))
ds <- ds.bestfit
length(unique(ds$sp))

figures.folder <- paste0("figures/traits/", growth.type, "/", level.folder)
if(!dir.exists(figures.folder)) {dir.create(figures.folder)}

# load(file = paste("results/commlevel/ds.bestfit_cor", file.extension.base4, ".Rdata", sep = ""))
# ds <- rbind(ds, ds.bestfit)
level_key <- c("tiny" = 1, "small" = 2, "medium" = 3, "large" = 4) #c(a = "apple", b = "banana", c = "carrot")

ds <- ds %>% mutate(tlplevel = as.factor(tlplevel)) %>%
  subset(!is.na(udi.best)) %>%
  mutate(size = factor(size, levels = c("tiny", "small", "medium", "large"))) %>%
  mutate(size.level = as.numeric(recode_factor(size, !!!level_key))) %>%
  droplevels()

deci <- read.csv("data-raw/traits/HydraulicTraits_Kunert/deciduous_species_Meakem.csv")
deci <- deci %>% mutate(sp = as.character(Species.code), Deciduousness = as.character(Deciduousness)) %>%
  select(sp, Deciduousness)
deci.level_key <- c("Evg" = 1, "DF" = 2, "DB" = 3, "DO" = 4, "D" = "4") #c(a = "apple", b = "banana", c = "carrot")

lwp <- read_excel("data-raw/traits/2016ENSO_Panama_LWP_20170407181307/2016ENSO_Panama_LWP.xlsx",
             sheet = "Panama_LWP")
lwp <- lwp %>% mutate(month = format(Date, "%b")) %>% subset(LWP_bar != -9999) ## same sp is not measured at two different locations
lwp.min <- lwp %>% group_by(Species) %>%
  summarise(lwp.min = -min(LWP_bar/10, na.rm = TRUE)) %>% ## in MPa
  mutate(sp = tolower(Species)) %>% select(-Species)

## Panama rainfall gradient preference
moist <- read.csv("data-raw/Condit_et_al_2013/TreeCommunityDrySeasonSpeciesResponse.csv")
moist <- moist %>% mutate(sp = paste0(tolower(str_sub(species, 1, 4)), str_sub(genus, 1, 2))) %>%
  rename(moist = Moist, moist.2 = Moist.2) %>% select(sp, moist, moist.2)
which(moist$sp == "pipeco") ## appears twice, so removing one.
moist <- moist[-372,]
hab.swp <- read.csv(file.path("data-raw/sp.plot.hab.swp.csv"))
sp.hab <- moist %>% full_join(hab.swp) %>%
  rename(Panama.Moist = moist, Panama.Moist2 = moist.2, Plot.swp = med.swp.reg, Plot.swp.ENSO = med.swp.dry) %>%
  select(-sd.swp.reg, -sd.swp.dry)
panama.traits <- read.csv("data-raw/traits/BCITRAITS_20101220.csv") %>%
  rename(form1 = GRWFRM1., form2 = GRWFRM2., sp = SP.) %>% mutate(sp = tolower(sp))
canopy.bci.traits <- panama.traits %>% subset(form1 %in% c("T")) %>% droplevels()
canopy.sp <- unique(canopy.bci.traits$sp)
udi.panama.traits <- panama.traits %>% select(-GENUS., -SPECIES., -FAMILY., -form2, -form1) %>%
  # full_join(lwp.min, by = "sp") %>%
  # full_join(deci, by = "sp") %>%
  # mutate(Deciduousness = replace_na(Deciduousness, "Evg")) %>%
  # mutate(DeciLvl = as.numeric(recode_factor(Deciduousness, !!!deci.level_key))) %>%
  # full_join(demo.sp %>% select(-mrate, -grate, -avg.abund), by = "sp") %>%
  # full_join(sp.hab, by = "sp") %>%
  full_join(ds %>%
              subset(size == "large" & sp %in% canopy.sp) %>%
              select(sp, udi.best), by = "sp") %>%
  subset(!is.na(udi.best))

plot.dst <- udi.panama.traits %>%
  select(-starts_with("LEAFLET")) %>%
  # gather(trait, value, -sp, -Deciduousness, -size, -udi.best) %>%
  gather(trait, value, -sp, -udi.best) %>%
  group_by(trait) %>%
  mutate(pval = summary(lm(udi.best ~ value))$coefficients[2,4],
         stars = cut(pval, breaks = c(0, 0.001, 0.01, 0.05, 0.1), labels = c("***", "**", "*", "`"),
                     include.lowest = FALSE, right = TRUE)) %>%
  subset(!is.na(udi.best))
title.piece = "Panama traits data"

formula <- y ~ x
p0 <- ggplot(plot.dst, aes(y = udi.best, x = value)) +
  geom_point() + ylab("Water Uptake Depth (m)") +
  geom_smooth(method = "lm", formula = formula) +
  facet_wrap( ~ trait, scales = "free_x") +
  ggtitle(paste0(title.piece, " | No. of Species = ", length(unique(plot.dst$sp)))) +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               npcx = 0.1, npcy = 0.1, rr.digits = 2,
               formula = formula, parse = TRUE, size = 4) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text_npc',
                  aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                  npcx = 0.85, npcy = 0.1, size = 4) +
  geom_text_npc(inherit.aes = FALSE,
                aes(npcx = 0.95, npcy = 0.05, label = stars), size = 6, color = "red", show.legend = FALSE) +
  ylim(c(max(plot.dst$udi.best, na.rm = TRUE) + 1), 0) + scale_y_reverse()
ggsave(file.path(paste0(figures.folder, "/udi_vs_panama.traits",
                      "splevel_large.jpeg")), plot = p0, height = 18, width = 24, units ='in')
p1 <- p0 %+% subset(plot.dst, !is.na(stars) & !trait %in% c("Chl", "root.95"))
ggsave(file.path(paste0(figures.folder ,"/udi_vs_panama.traits",
                        "splevel_large_signif.jpeg")), plot = p1, height = 9, width = 12, units ='in')

hyd <- read.csv("data-raw/traits/HydraulicTraits_BrettWolfe/ht1_20200103.csv") #  # Brett's data
hyd <- hyd %>% select(-genus, -species, -deciduousness, -site) %>%
  rename(AlAsw = laxa_m2m2_m, LMABrett = lma_gm2_m, WDBrett = xylem_den_m, TLPBrett = tlp_m, p80S = p80, p88S = p88,
         Fcap = cwr_xylem_elbow_origin_slope, CWR = cwr_xylem_cwr_at_elbow, Felbow = cwr_xylem_elbow,
         LDMC = ldmc_m, BarkThick = barkthickness10mm) %>%
  ## Given Fcap slope values are -ve but they should be positive
  mutate(sp = tolower(sp), Fcap = -Fcap) %>%
  left_join(lwp.min, by = "sp") %>%
  mutate(HSMTLP.88S =  TLPBrett - p88S, HSM88S = lwp.min - p88S,
         HSMTLPBrett = lwp.min - TLPBrett, HSMFelbow = lwp.min - Felbow) %>%
  left_join(deci, by = "sp") %>%
  mutate(Deciduousness = replace_na(Deciduousness, "Evg")) %>%
  mutate(DeciLvl = as.numeric(recode_factor(Deciduousness, !!!deci.level_key))) %>%
  left_join(demo.sp %>% select(-mrate, -grate, -avg.abund), by = "sp") %>%
  left_join(sp.hab, by = "sp")
length(unique(hyd$sp)) # 27 sp across BCI, PNM, San Lorenzo
## which of these species present in BCI so with splevel estimates
length(match(unique(hyd$sp), unique(ds$sp))[!is.na(match(unique(hyd$sp), unique(ds$sp)))]) # 15

traits.indi <- read.csv("data-raw/traits/HydraulicTraits_Kunert/hydraulic_traits_panama_kunert.csv") # Nobby's data
traits <- traits.indi %>% group_by(sp) %>%
  select(-idividual, -ind_ID,
         -leaf_area_dry_cm2, -leaf_area_fresh_cm2,
         -sum_dry_mass_leaf_blade_g, -sum_fresh_mass_leaf_blade_g, -PLA_dry_percent) %>%
  summarise_all(mean, na.rm = TRUE)
colnames(traits) <- c("sp", "LA", "LMANob", "LT", "LD", "WMA", "SPAD", "Chl", "TLPNob", "WDNob")
traits <- traits %>%
  left_join(lwp.min, "sp") %>%
  left_join(demo.sp %>% select(-mrate, -grate, -avg.abund), by = "sp")
leaf_cond.models <- read.csv("data-raw/traits/HydraulicTraits_Kunert/Panama_fits_leaf_K_p50_Kunert.csv")
leaf.k.p80 <- leaf_cond.models %>% subset(model_type == "Exponential") %>%
  mutate(sp = data.type, p80L = -psi_kl80, KmaxL = Kmax) %>% # these are Kmax that are extrapolated from the exponential curve
  select(sp, KmaxL, p80L) # 21 species
traits <- traits %>% left_join(leaf.k.p80 %>% mutate(sp = as.character(sp)), by = "sp") %>%
  mutate(HSMTLP.80L = TLPNob - p80L, HSM80L = lwp.min - p80L, HSMTLPNob = lwp.min - TLPNob) %>%
  left_join(deci, by = "sp") %>%
  mutate(Deciduousness = replace_na(Deciduousness, "Evg"))  %>%
  mutate(DeciLvl = as.numeric(recode_factor(Deciduousness, !!!deci.level_key))) %>%
  left_join(sp.hab, by = "sp")

#####----------traits charts alone-----------
with_kmax <- "yes"
df <- dplyr::select_if(traits, is.numeric)
if(with_kmax == "no"){
  df <- df %>% select(-KmaxL, -p80L, -HSM80L)
}
cor.mat <- cor(df, use = "complete.obs")
ggcorrplot(cor.mat,
           hc.order = TRUE,
           type = "full",
           lab = TRUE,
           title = paste0("No. of species = ", nrow(df[complete.cases(df),])))
ggsave(file.path(paste0("kunert_traits_corrplot_chart_with_kmax_", with_kmax,".jpeg")))
## Correlation chart
pdf(file.path(paste0("kunert_traits_corrplot_chart_with_kmax_", with_kmax,"_chart.pdf")))
chart.Correlation(df, histogram=TRUE, pch=19) # does not take title
dev.off()
#####----------traits-----------

## Heatmap by deciduousness
# Heatmap
traits.long <- traits %>% select(-DeciLvl) %>%
  left_join(ds %>% subset(size %in% c("large")) %>% select(sp, udi.best), by = "sp") %>%
  gather(trait, value, -sp, -Deciduousness)

# traits.long <- traits %>% mutate_all(as.character) %>% pivot_longer(traits, cols = -c("sp", "Deciduousness"), names_to = "trait",
#                                        values_to = "value")
# ggplot(traits.long, aes(x = value, color = Deciduousness)) +
#   geom_histogram(stat = "identity") +
#   #ylab("Trait") + xlab("Deciduousness") +
#   facet_grid(trait ~ ., scales = "free_y") +
#   ggtitle("Kunert's trait data by deciduousness")
# ggsave(file.path(paste0(figures.folder ,"/sp_udi.best_vs_deciduousness.jpeg")), height = 5, width = 5, units ='in')
#
### With Tukey labels:
# xx <- traits.long %>% subset(trait == "SPAD")
# kruskal(xx$value, xx$Deciduousness, group=TRUE, p.adj="bonferroni")$groups
kruskal.list <- list()
for(i in unique(traits.long$trait)) {
  xx <- traits.long %>% subset(trait == i)
  kruskal.list[[i]] <- cbind(trait = i, kruskal(xx$value, xx$Deciduousness, alpha = 0.1, group=TRUE, p.adj="bonferroni")$groups,
                             Deciduousness = rownames(kruskal(xx$value, xx$Deciduousness, alpha = 0.1, group=TRUE, p.adj="bonferroni")$groups))
}
traits.kruskal.labels <- do.call(rbind.data.frame, kruskal.list)
head(traits.kruskal.labels)
# generate_label_df <- function(TUKEY, variable){
#   # Extract labels and factor levels from Tukey post-hoc
#   Tukey.levels <- TUKEY[[variable]][,4]
#   Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
#   #I need to put the labels in the same order as in the boxplot :
#   Tukey.labels$treatment=rownames(Tukey.labels)
#   Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
#   return(Tukey.labels)
# }
# traits.labels <- traits.long %>% group_by(trait) %>% arrange(trait) %>%
#   do(mod = lm(value ~ Deciduousness, data = .)) %>%
#   do(aov.mod = aov(.$mod)) %>%
#   do(TUKEY = TukeyHSD(.$aov.mod, conf.level = 0.95)) %>%
#   do(labels = generate_label_df(.$TUKEY, "Deciduousness")) %>%
#   ## since we arranged tibble by trait in the first row
#   add_column(trait = sort(unique(traits.long$trait))) %>%
#   select(trait, labels) %>% unnest(cols = c(labels)) %>%
#   rename(Deciduousness = treatment)
traits.labels <- traits.kruskal.labels
traits.labels.data <- traits.labels %>%
  left_join(traits.long %>% group_by(trait) %>%
              summarise(value = max(value, na.rm = TRUE)), by = c("trait"))
ggplot(traits.labels.data, aes(x = Deciduousness, y = value)) +
  facet_wrap(. ~  trait, scales = "free_y") +
  geom_text(aes(label = groups, color = Deciduousness), vjust = 1, hjust = 0) +
  geom_boxplot(data = traits.long, aes(fill = Deciduousness), stat = "boxplot", notch = TRUE)
ggsave(file.path(figures.folder, paste0("traits_vs_deciduousness.jpeg")), height = 9, width = 12, units ='in')

hyd.long <- hyd %>%
  select(sp, Deciduousness, AlAsw, TLPBrett, p80S, p88S, Fcap, CWR, Felbow, WDBrett, LMABrett, LDMC, BarkThick,
         HSM88S,  HSMTLPBrett, HSMFelbow, HSMTLP.88S,
         grate.adult, grate.juve, mrate.adult, mrate.juve, Panama.Moist, Panama.Moist2, Plot.swp, Plot.swp.ENSO, lwp.min) %>%
  left_join(ds %>% subset(size %in% c("large")) %>% select(sp, udi.best), by = "sp") %>%
  gather(trait, value, -sp, -Deciduousness)

# test.aov <- generate_label_df(TukeyHSD(aov(lm(value ~ Deciduousness, data= hyd.long %>% subset(trait == "BarkThick")))), "Deciduousness")
# hyd.labels <- hyd.long %>% group_by(trait) %>% arrange(trait) %>%
#   do(mod = lm(value ~ Deciduousness, data = .)) %>%
#   do(aov.mod = aov(.$mod)) %>%
#   do(TUKEY = TukeyHSD(.$aov.mod, conf.level = 0.95)) %>%
#   do(labels = generate_label_df(.$TUKEY, "Deciduousness")) %>%
#   ## since we arranged tibble by trait in the first row
#   add_column(trait = sort(unique(hyd.long$trait))) %>%
#   select(trait, labels) %>% unnest(cols = c(labels)) %>%
#   rename(Deciduousness = treatment)
kruskal.list <- list()
for(i in unique(hyd.long$trait)) {
  xx <- hyd.long %>% subset(trait == i)
  kruskal.list[[i]] <- cbind(trait = i, kruskal(xx$value, xx$Deciduousness, alpha = 0.1, group=TRUE, p.adj="bonferroni")$groups,
                             Deciduousness = rownames(kruskal(xx$value, xx$Deciduousness, alpha = 0.1, group=TRUE, p.adj="bonferroni")$groups))
}
hyd.kruskal.labels <- do.call(rbind.data.frame, kruskal.list)
head(hyd.kruskal.labels)
hyd.labels <- hyd.kruskal.labels
hyd.labels.data <- hyd.labels %>%
  left_join(hyd.long %>% group_by(trait) %>%
              summarise(value = max(value, na.rm = TRUE)), by = c("trait"))

ggplot(hyd.labels.data, aes(x = Deciduousness, y = value)) +
  facet_wrap(. ~  trait, scales = "free_y") +
  geom_text(aes(label = groups, color = Deciduousness), vjust = 1, hjust= 0) +
  geom_boxplot(data = hyd.long, aes(fill = Deciduousness), stat = "boxplot", notch = TRUE)
ggsave(file.path(figures.folder, paste0("hyd_traits_vs_deciduousness.jpeg")), height = 9, width = 12, units ='in')


#####----------traits plot end-----------
size.class <- levels(ds$size)
demo <- c("off", "on")

for (j in 1: length(tlplevels)){
  for (i in 1: length(size.class)){
    ds.tlp <- ds %>% left_join(traits %>% select(TLP, sp), by = "sp") %>% subset(size.class == size.class[i])
    lm.1 <- lm(udi ~ TLP, data = ds.tlp); summ.lm.1 <- summary(lm.1)
    lm.1.label <- paste0("Water Uptake Depth = ", round(summ.lm.1$coefficients[2, 1], 2), " * TLP + ", round(summ.lm.1$coefficients[1, 1], 2),
                     "\nR-squared = ", round(summ.lm.1$r.squared, 3),
                  ", p-val = ", round(summ.lm.1$coefficients[2, 4], 2))
    ggplot(ds.tlp,
      aes(x = TLP, y = udi.best)) +
      geom_point(size = 2, aes(shape = Deciduousness, color = Deciduousness)) +
      scale_y_continuous(trans="rev_sqrt", breaks = soil.depths) +
      geom_smooth(aes(group = Deciduousness, color = Deciduousness), method = "lm", se= FALSE, size = 0.5) +
      geom_smooth(method = "lm", se = FALSE, color = "black", show.legend = FALSE, size = 0.5) +
      ylab("Water Uptake Depth (m)") + xlab("Turgor Loss Point [-MPa]") +
      geom_text(aes(x = TLP + 0.06, label = sp, color = Deciduousness), size = 2) +
      ggtitle(paste0("Species Uptake Depth Vs TLP\nSize Class = ", size.class[i], "\n", lm.1.label))
    ggsave(file.path(paste0(figures.folder ,"/sp_udi.best_vs_tlp_by_dec_", size.class[i], ".jpeg")), height = 6, width = 6, units ='in')
  }
}

length(unique(traits$sp)) # 51 species
figures.folder.deci <- paste0(figures.folder, "/deciduousness/")
if(!dir.exists(figures.folder.deci)) {dir.create(figures.folder.deci)}
if(!dir.exists(paste0(figures.folder, "/pca/"))) {dir.create(paste0(figures.folder, "/pca/"))}
if(!dir.exists(paste0(figures.folder, "/chart/"))) {dir.create(paste0(figures.folder, "/chart/"))}

with_kmax <- "yes"
for (j in 1: length(tlplevels)){
  for (i in 1: length(size.class)){
    for (k in 1:length(demo)){
      ds.sp <- ds %>% select(sp, size, udi.best, root.95, tlplevel) %>% subset(tlplevel == tlplevels[j] & size == size.class[i])
    ## when udi.best for large size is absent, using medium size data
      if (i == 4) {
        ds.sp.smaller <- ds %>%
          subset(tlplevel == tlplevels[j] & size == size.class[i - 1] &
                   ## only retaining species for which ds.sp did not have a udi.best
                   sp %in% traits$sp[-match(ds.sp$sp, traits$sp)]) %>%
          select(sp, size, udi.best, root.95,  tlplevel)
        ds.sp <- ds.sp %>% rbind(ds.sp.smaller)
      }
      dst <- left_join(traits, ds.sp, by = "sp") %>% rename(udi.best = udi.best)
       # 51 sp for i = 4, but 17 with complete cases
       if(with_kmax == "no"){
        dst <- dst %>% select(-KmaxL, -p80L, -HSM80L)
       }
      ### plot of udi vs other traits-----
      if( with_kmax == "no") { title.piece = "Nobby's Data without Kmax"} else { title.piece = "Nobby's Data with Kmax"}
      plot.dst <- dst %>%
        gather(trait, value, -sp, -Deciduousness, -size, -tlplevel, -udi.best) %>%
        group_by(trait) %>%
        mutate(pval = summary(lm(udi.best ~ value))$coefficients[2,4],
             stars = cut(pval, breaks = c(0, 0.001, 0.01, 0.05, 0.1), labels = c("***", "**", "*", "`"),
                         include.lowest = FALSE, right = TRUE)) %>%
        subset(!is.na(udi.best))
      formula <- y ~ x
      p0 <- ggplot(plot.dst, aes(y = udi.best, x = value)) +
        geom_point() + ylab("Water Uptake Depth (m)") +
        geom_smooth(method = "lm", formula = formula) +
        facet_wrap( ~ trait, scales = "free_x") +
        ggtitle(paste0(title.piece, " | No. of Species = ", length(unique(plot.dst$sp)))) +
        stat_poly_eq(aes(label = paste(..rr.label..)),
                     npcx = 0.1, npcy = 0.1, rr.digits = 2,
                     formula = formula, parse = TRUE, size = 4) +
        stat_fit_glance(method = 'lm',
                        method.args = list(formula = formula),
                        geom = 'text_npc',
                        aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                        npcx = 0.85, npcy = 0.1, size = 4) +
        geom_text_npc(inherit.aes = FALSE,
                      aes(npcx = 0.95, npcy = 0.05, label = stars), size = 6, color = "red", show.legend = FALSE) +
        ylim(c(max(plot.dst$udi.best, na.rm = TRUE) + 1), 0) + scale_y_reverse()
      ggsave(file.path(paste0(figures.folder ,"/udi_vs_traits",
                              tlplevels[j], "_", size.class[i], "_demo_", demo[k], "_chart_with_kmax_",
                              with_kmax, ".jpeg")), plot = p0, height = 9, width = 12, units ='in')
      p1 <- p0 %+% subset(plot.dst, !is.na(stars) & !trait %in% c("Chl", "root.95"))
      ggsave(file.path(paste0(figures.folder ,"/udi_vs_traits",
                              tlplevels[j], "_", size.class[i], "_demo_", demo[k], "_chart_with_kmax_",
                              with_kmax, "_signif.jpeg")), plot = p1, height = 6, width = 7, units ='in')

      ####---traits by deciduousness-----
      traits.long <- dst %>%
        gather(trait, value, -sp, -Deciduousness, -size, -tlplevel, -DeciLvl)
      kruskal.list <- list()
      for(n in unique(traits.long$trait)) {
        xx <- traits.long %>% subset(trait == n)
        kruskal.list[[n]] <- cbind(trait = n, kruskal(xx$value, xx$Deciduousness, alpha = 0.1, group=TRUE, p.adj="bonferroni")$groups,
                                   Deciduousness = rownames(kruskal(xx$value, xx$Deciduousness, alpha = 0.1, group=TRUE, p.adj="bonferroni")$groups))
      }
      traits.kruskal.labels <- as.data.frame(do.call(rbind.data.frame, kruskal.list))
      traits.labels.data <- traits.kruskal.labels %>%
        left_join(traits.long %>% group_by(trait) %>%
                    summarise(value = max(value, na.rm = TRUE)), by = c("trait"))
      ggplot(traits.labels.data, aes(x = Deciduousness, y = value)) +
        facet_wrap(. ~  trait, scales = "free_y") +
        geom_text(aes(label = groups, color = Deciduousness), vjust = 1, hjust = 0) +
        geom_boxplot(data = traits.long, aes(fill = Deciduousness), stat = "boxplot")
      ggsave(file.path(paste0(figures.folder.deci, "traits_vs_deciduousness_notch_",
                              tlplevels[j], "_", size.class[i], "_demo_", demo[k], "_chart_with_kmax_", with_kmax, ".jpeg")),
             height = 9, width = 12, units ='in')
      ###
      ### PCA ---------
      df.pca <- dst %>% remove_rownames %>% column_to_rownames(var = "sp") %>% select_if(is.numeric) ## more species without Kmax data %>% select(-SafetyMargin.p50, -p50_Kmax, -Kmax)
      # if (diff(range(df.pca$root.95, na.rm = TRUE)) == 0) {
      df.pca <- df.pca %>% subset(complete.cases(df.pca)) %>% select(-root.95)
      result.pca <- prcomp(df.pca, center = TRUE, scale = TRUE)
      if( with_kmax == "no") { title.piece = "Nobby's Data without Kmax"} else { title.piece = "Nobby's Data with Kmax"}
      pdf(file.path(paste0(figures.folder ,"/pca/sp_udi.best_traits_demo",
                           tlplevels[j], "_", size.class[i], "_with_kmax_", with_kmax, "_axes_1_2.pdf")), height = 8, width = 8)
      biplot(result.pca, choices = 1:2, pc.biplot = TRUE, main = paste0(title.piece, " | No. of Species = ", nrow(df.pca)))
      dev.off()
      pdf(file.path(paste0(figures.folder ,"/pca/sp_udi.best_traits_demo",
                           tlplevels[j], "_", size.class[i], "_with_kmax_", with_kmax, "_axes_3_4.pdf")), height = 8, width = 8)
      biplot(result.pca, choices = 3:4, pc.biplot = TRUE, main = paste0(title.piece, " | No. of Species = ", nrow(df.pca)))
      dev.off()
      ### Correlation pair plots ---------
      df <- dplyr::select_if(dst, is.numeric) %>% subset(complete.cases(dst))## more species without Kmax data %>% select(-SafetyMargin.p50, -p50_Kmax, -Kmax)
      cor.mat <- cor(df, use = "complete.obs")
      g.deci <- ggcorrplot(cor.mat,
                 hc.order = FALSE,
                 type = "full",
                 lab = TRUE,
                 title = paste0("TLPlevel = ", tlplevels[j], "\n Tree Size Class = ", size.class[i],
                                "\nNo. of species = ", nrow(df[complete.cases(df),])))
      ggsave(plot = g.deci, file.path(paste0(figures.folder ,"/sp_udi.best_vs_traits_corrplot_tlp",
                              tlplevels[j], "_", size.class[i], "_demo_", demo[k], "_with_kmax_", with_kmax,".jpeg")), height = 10, width = 10, units ='in')
      ## Correlation chart
      pdf(file.path(paste0(figures.folder ,"/chart/sp_udi.best_vs_traits_corrplot_tlp",
                           tlplevels[j], "_", size.class[i], "_demo_", demo[k], "_chart_with_kmax_", with_kmax, ".pdf")))
      chart.Correlation(df, histogram=TRUE, pch=19)# does not take title
      dev.off()
    }
  }
}

figures.folder.hyd <- paste0(figures.folder, "/hyd/")
if(!dir.exists(figures.folder.hyd)) {dir.create(figures.folder.hyd)}
if(!dir.exists(paste0(figures.folder.hyd, "/pca/"))) {dir.create(paste0(figures.folder.hyd, "/pca/"))}
if(!dir.exists(paste0(figures.folder.hyd, "/chart/"))) {dir.create(paste0(figures.folder.hyd, "/chart/"))}
if(!dir.exists(paste0(figures.folder.hyd, "/sub/"))) {dir.create(paste0(figures.folder.hyd, "/sub/"))}
if(!dir.exists(paste0(figures.folder.hyd, "/full/"))) {dir.create(paste0(figures.folder.hyd, "/full/"))}
if(!dir.exists(paste0(figures.folder.hyd, "/sub/deciduousness/"))) {dir.create(paste0(figures.folder.hyd, "/sub/deciduousness/"))}
if(!dir.exists(paste0(figures.folder.hyd, "/full/deciduousness/"))) {dir.create(paste0(figures.folder.hyd, "/full/deciduousness/"))}

data.subset <- c("sub", "full")
q <- 1 #for data.subset[q]
all.traits <- c("off", "on")
## hyd traits
for (j in 1: length(tlplevels)){
  for (i in 1: length(size.class)){
    for (m in 1: length(all.traits)){
      if (q == 1) {
        hyd <- hyd %>%
          select(sp, Deciduousness, DeciLvl, AlAsw, TLPBrett, p80S, p88S, Fcap, CWR, Felbow,
                 WDBrett, LMABrett, LDMC, BarkThick, HSM88S,  HSMTLPBrett, HSMFelbow, HSMTLP.88S,
                 grate.adult, grate.juve, mrate.adult, mrate.juve, Panama.Moist, Panama.Moist2, Plot.swp, Plot.swp.ENSO)
          }
      ds.sp <- ds %>% select(sp, size, udi.best, root.95, tlplevel) %>% subset(tlplevel == tlplevels[j] & size == size.class[i])
      if (i == 4) {
        # species for which ds.sp did not have a udi.best
        sp.needed <- hyd$sp[!hyd$sp %in% ds.sp$sp]; sp.needed <- sp.needed[!is.na(sp.needed)] # 16 sp for i = 4
        ds.sp.smaller <- ds %>%
          subset(tlplevel == tlplevels[j] & size == size.class[i - 1] &
                   ## only retaining species for which ds.sp did not have a udi.best
                   sp %in% sp.needed) %>%
          select(sp, size, udi.best, root.95, tlplevel)
        ds.sp <- ds.sp %>% rbind(ds.sp.smaller) ## 18 sp for i = 4
      }
      dst <- left_join(hyd, ds.sp, by = "sp") %>%
        subset(!is.na(udi.best)) # 14 sp for i = 4

      if(all.traits[m] == "on") {
        dst <- dst %>%
          left_join(traits %>% select(-c(grate.adult, grate.juve, mrate.adult,
                                         mrate.juve, Panama.Moist, Panama.Moist2, Plot.swp, Plot.swp.ENSO, DeciLvl, Deciduousness)), by = "sp")
      } ## 14 sp for i = 4, but 9 with complete cases
      write.csv(dst, "results/dst.csv")
      ### plot of udi vs other traits-----
      if(all.traits[m] == "on") { title.piece = "Brett & Nobby's Data"} else { title.piece = "Brett's Data"}
      plot.dst <- dst %>%
        gather(trait, value, -sp, -Deciduousness, -size, -tlplevel, -udi.best) %>%
        group_by(trait) %>%
        mutate(pval = summary(lm(udi.best ~ value))$coefficients[2,4], # round(cor.test(value, udi.best, use = "complete.obs")$p.value, 2),
               stars = cut(pval, breaks = c(0, 0.001, 0.01, 0.05, 0.1), labels = c("***", "**", "*", "`"),
                           include.lowest = FALSE, right = TRUE))
      formula <- y ~ x
      p0 <- ggplot(plot.dst, aes(y = udi.best, x = value)) +
        geom_point() + ylab("Water Uptake Depth (m)") +
        geom_smooth(method = "lm", formula = formula) +
        facet_wrap( ~ trait, scales = "free_x") +
        ggtitle(paste0(title.piece, " | No. of Species = ", length(unique(plot.dst$sp)))) +
        stat_poly_eq(aes(label = paste(..rr.label..)),
                     npcx = 0.2, npcy = 0.1, rr.digits = 2,
                     formula = formula, parse = TRUE, size = 3) +
        stat_fit_glance(method = 'lm',
                        method.args = list(formula = formula),
                        geom = 'text_npc',
                        aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")),
                        npcx = 0.7, npcy = 0.1, size = 3) +
        geom_text_npc(inherit.aes = FALSE,
                      aes(npcx = 0.8, npcy = 0.07, label = stars), size = 6, color = "red", show.legend = FALSE) +
        ylim(c(max(plot.dst$udi.best, na.rm = TRUE) + 1), 0) + scale_y_reverse()
      ggsave(file.path(paste0(figures.folder.hyd, "udi_vs_traits_",
                              tlplevels[j], "_", size.class[i],"_", data.subset[q], "_all.traits_",
                              all.traits[m],"_demo.jpeg")), plot = p0,
             # height = 9, width = 12, units ='in')
             height = 8*(1 + 0.5*(m-1)), width = 10.5*(1 + 0.5*(m-1)), units ='in')
      if (nrow(subset(plot.dst, !is.na(stars) & !trait %in% c("Chl", "root.95"))) != 0) {
        p1 <- p0 %+% subset(plot.dst, !is.na(stars) & !trait %in% c("Chl", "root.95"))
        ggsave(file.path(paste0(figures.folder.hyd, "udi_vs_traits_",
                                tlplevels[j], "_", size.class[i],"_", data.subset[q], "_all.traits_",
                                all.traits[m],"_demo_signif.jpeg")), plot = p1,
               height = 5, width = 7, units ='in')
      }

      ###---traits by deciduousness-----
      traits.long <- dst %>%
        gather(trait, value, -sp, -Deciduousness, -size, -tlplevel, -DeciLvl)
      kruskal.list <- list()
      for(n in unique(traits.long$trait)) {
        xx <- traits.long %>% subset(trait == n)
        kruskal.list[[n]] <- cbind(trait = n, kruskal(xx$value, xx$Deciduousness, alpha = 0.1, group=TRUE, p.adj="bonferroni")$groups,
                                   Deciduousness = rownames(kruskal(xx$value, xx$Deciduousness, alpha = 0.1, group=TRUE, p.adj="bonferroni")$groups))
      }
      traits.kruskal.labels <- do.call(rbind.data.frame, kruskal.list)
      traits.labels.data <- traits.kruskal.labels %>%
        left_join(traits.long %>% group_by(trait) %>%
                    summarise(value = max(value, na.rm = TRUE)), by = c("trait"))
      ggplot(traits.labels.data, aes(x = Deciduousness, y = value)) +
        facet_wrap(. ~  trait, scales = "free_y") +
        geom_text(aes(label = groups, color = Deciduousness), vjust = 1, hjust = 0) +
        geom_boxplot(data = traits.long, aes(fill = Deciduousness), stat = "boxplot")
      ggsave(file.path(figures.folder.hyd, paste0(data.subset[q], "/deciduousness/traits_vs_deciduousness_notch_",
                              tlplevels[j], "_", size.class[i],"_", data.subset[q], "_all.traits_", all.traits[m],"_demo.jpeg")), height = 7, width = 9, units ='in')
      ### PCA ---------
      df.pca <- dst %>% remove_rownames %>% column_to_rownames(var = "sp") %>%
        select_if(is.numeric) %>% subset(complete.cases(dst)) %>% select(-root.95)## more species without Kmax data %>% select(-SafetyMargin.p50, -p50_Kmax, -Kmax)
      result.pca <- prcomp(df.pca, center = TRUE, scale = TRUE)
      pdf(file.path(paste0(figures.folder.hyd, "pca/sp_udi.best_traits_demo",
                           tlplevels[j], "_", size.class[i], "_", data.subset[q], "_all.traits_", all.traits[m], "_axes_1_2.pdf")), height = 8, width = 8)
      biplot(result.pca, choices = 1:2, pc.biplot = TRUE, main = paste0(title.piece, " | No. of Species = ", nrow(df.pca)))
      dev.off()
      pdf(file.path(paste0(figures.folder.hyd, "sp_udi.best_traits_demo",
                           tlplevels[j], "_", size.class[i], "_", data.subset[q], "_all.traits_", all.traits[m], "_axes_3_4.pdf")), height = 8, width = 8)
      biplot(result.pca, choices = 3:4, pc.biplot = TRUE, main = paste0(title.piece, " | No. of Species = ", nrow(df.pca)))
      dev.off()
      ### Correlation pair plots ---------
      df <- dplyr::select_if(dst, is.numeric) %>% subset(complete.cases(dst)) %>% select(-root.95)## more species without Kmax data %>% select(-SafetyMargin.p50, -p50_Kmax, -Kmax)
      cor.mat <- cor(df, use = "complete.obs")
      g.deci <- ggcorrplot(cor.mat,
                 hc.order = TRUE,
                 type = "full",
                 lab = TRUE,
                 title = paste0("TLPlevel = ", tlplevels[j], "\n Tree Size Class = ", size.class[i],
                                "\nNo. of species = ", nrow(df[complete.cases(df),])))
      ggsave(plot = g.deci, file.path(paste0(figures.folder.hyd, data.subset[q] ,"/sp_udi.best_vs_traits_corrplot_tlp",
                              tlplevels[j], "_", size.class[i],"_", data.subset[q], "_all.traits_", all.traits[m],"_demo.jpeg")),
             height = 10*(1 + 0.5*(m-1)), width = 10*(1 + 0.5*(m-1)), units ='in')
      ## Correlation chart
      pdf(file.path(paste0(figures.folder.hyd, "chart/sp_udi.best_vs_traits_corrplot_tlp",
                           tlplevels[j], "_", size.class[i], "_", data.subset[q], "_all.traits_", all.traits[m], "_demo_chart.pdf")),
          height = 10*(1 + 0.5*(m-1)), width = 10*(1 + 0.5*(m-1)))
      chart.Correlation(df, histogram=TRUE, pch=19) # does not take title
      dev.off()
    }
  }
}

