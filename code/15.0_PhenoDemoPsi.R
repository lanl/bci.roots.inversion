#---------------------------------
# Title: Phenology, demography, water availability
# Author : Rutuja Chitra-Tarak
# Original date: April 25, 2020
#---------------------------------

## 1. How is species phenology located on the fast-slow continuum?
## 1.a Are deciduous species fast while evergreen species slow?
## 1.b demographically fast == high growth and mortality rates; slow == low growth and mortality rates
## 1.c fast traits high SLA, high Kmax, high TLP, deep roots
## 1.d fast species located on greater resource environments:
##        wetter sites within 50-ha plot, wet distributed along Panama gradient

rm(list=ls())

if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, readxl, forcats, agricolae, gridExtra, scales)
# graphics info
theme_set(theme_bw())
theme_update(text = element_text(size = 14),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)

rev_sqrt_trans <- function() {
  scales::trans_new(
    name = "rev_sqrt",
    transform = function(x) -sqrt(abs(x)),
    inverse = function(x) x^2);
}

n.threshold = 50
figures.folder <- paste0("figures/mortality/phenology")

## Deciduousness-----
deci <- read_excel(file.path("data-raw/traits/nomenclature_R_20190524_Rready_Osvaldo Calderon & JoeWright_expert_opinion.xlsx"))
deci.level_key <- c("Evg" = 1, "DF" = 2, "DB" = 3, "DO" = 4, "D" = "4") #c(a = "apple", b = "banana", c = "carrot")

deci <- deci %>% mutate(sp = tolower(sp6)) %>%
  select(sp4, sp, deciduous) %>%
  subset(deciduous %in% c("E", "DF", "DO", "DB")) %>%
  mutate(deciduousness = recode_factor(as.factor(deciduous), `E` = "Evergreen", `DO` = "Obligate Deciduous", `DF` = "Facultative Deciduous",
                                       `DB` = "Brevideciduous"), ordered = TRUE) %>%
  mutate(DeciLvl = as.numeric(recode_factor(deciduous, !!!deci.level_key))) %>%
  select(sp4, sp, deciduous, deciduousness, DeciLvl)

head(deci)
# deci.2 <- read.csv("data-raw/traits/HydraulicTraits_Kunert/deciduous_species_Meakem.csv")
# deci.2 <- deci.2 %>% mutate(sp = as.character(Species.code), deciduousness = as.character(Deciduousness)) %>%
#   select(sp, deciduousness)

### LWP -----

lwp <- read_excel("data-raw/traits/2016ENSO_Panama_LWP_20170407181307/2016ENSO_Panama_LWP.xlsx",
                  sheet = "Panama_LWP")
lwp <- lwp %>%
  mutate(month = format(Date, "%b"),
         date.chr = as.character(Date),
         date = as.Date(Date),
         LWP_date.time.chr = paste(date.chr, LWP_coll_time),
         LWP_date.time = as.POSIXct(LWP_date.time.chr, tz = "America/Panama", format = "%Y-%m-%d %H:%M")) %>%
  subset(LWP_bar != -9999) %>% select(-date.chr, -LWP_date.time.chr, -Date)
head(lwp)
## same sp is not measured at two different locations
lwp.all <- lwp %>%
  mutate(sp = tolower(Species),
         # Brett's e-mail: ALBIED should be ALBIAD for Albizia adinocephala. We were using ALBIED mistakenly at the beginning
         sp = ifelse(sp == "albied", "albiad", sp)) %>%
  select(-Species) %>%
  rename(location = Location) %>%
  group_by(location, date, sp, Gasex_typeOrSequence, LWP_coll_time) %>%
  summarise(lwp.min = -min(LWP_bar/10, na.rm = TRUE)) %>% ## in MPa
  mutate(LWP_coll_time = as.numeric(LWP_coll_time)) %>%
  unite("sp_date", sp, date, remove = FALSE) %>%
  data.frame() %>%
  transform(location = factor(location, levels = c("PA-PNM", "PA-BCI", "PA-SLZ"), ordered = TRUE)) %>%
  rename(time = Gasex_typeOrSequence) %>%
  unite("sp_date_time", sp, date, time, remove = FALSE)


str(lwp.all)

ggplot(lwp.all, aes(x = LWP_coll_time, y = lwp.min)) +
  facet_grid(location ~ .) +
  geom_point(aes(color = as.factor(date))) +
  geom_line(aes(group = c(sp_date), color = as.factor(date)), show.legend = FALSE) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE) +
  guides(color = guide_legend(title = "Date"))
ggsave(file.path(figures.folder, paste0("LWP_time_series_all_Data.jpeg")), height = 6, width = 5, units ='in')

# Across all the dryseason measurements which diurnal measurement was minimum by species and location
## also get the predawn for the same day

## Most of the trees have multiple records per time step
lwp.all.indi <- lwp.all %>%
  group_by(location, sp, date, time, sp_date) %>%
  summarise(lwp.min = mean(lwp.min, na.rm = TRUE))

# First getting those species-(location)-wise dates and date-times with min.LWP
lwp.min.diurnal.dates <- lwp.all.indi %>%
  subset(!time %in% "Pre-dawn LWP") %>%
  group_by(sp, location, sp_date) %>%
  # Location not really needed as this data set does not have same sp measured at two different locations
  filter(lwp.min == min(lwp.min, na.rm = TRUE)) %>%
  unite("sp_date_time", sp, date, time, remove = FALSE)

## getting predawn LWPs for the same day as that of those minimum diurnal LWPs
## Retain multiple records per tree to calculate sd
lwp.min.diurnal <- lwp.all %>%
  subset(sp_date_time %in% lwp.min.diurnal.dates$sp_date_time) %>%
  group_by(sp, location) %>%
  # Location (sp_date_location) not really needed as this dat asets does not have same sp measured at two different locations
  bind_rows(lwp.all %>% subset(time == "Pre-dawn LWP" &
                                 # Location (sp_date_location) not really needed as this dat asets does not have same sp measured at two different locations
                                 sp_date %in% lwp.min.diurnal.dates$sp_date)) %>%
  mutate(time = ifelse(time == "Pre-dawn LWP", "Predawn", "Diurnal")) %>%
  mutate(time = recode_factor(time, `predawn` = "Predawn", `diurnal` = "Diurnal"), ordered = TRUE) %>%
  group_by(sp, time, location) %>%
  summarise(lwp.se = sd(lwp.min, na.rm = TRUE)/sqrt(n()),
            ## don't change this sequence
            lwp.min = mean(lwp.min, na.rm = TRUE)) %>%
  ungroup(sp, time, location) %>%
  arrange(sp, location, time) %>%
  group_by(sp, location) %>%
  mutate(lwp.diff = lwp.min - lag(lwp.min)) %>%
  ungroup(sp, location)

## also, to calculate difference between preadawn and min.dirunal measurement across all measurements
## first calculate minimun dirunal measurement for each day
# first get which time had the minimum water potential on each measurement day by species
lwp.diff.filter <- lwp.all %>%
  subset(!time %in% "Pre-dawn LWP") %>%
  group_by(sp_date) %>%
  # Location not really needed as this data set does not have same sp measured at two different locations
  slice(which.min(lwp.min)) %>% # selects only the match for the first  min
  arrange(sp, date, time)

lwp.diff <- lwp.all %>%
  subset(sp_date_time %in% lwp.diff.filter$sp_date_time) %>%
  bind_rows(lwp.all %>%
              subset(time == "Pre-dawn LWP")) %>%
  mutate(time = ifelse(time == "Pre-dawn LWP", "Predawn", "Diurnal")) %>%
  mutate(time = recode_factor(time, `predawn` = "Predawn", `diurnal` = "Diurnal"), ordered = TRUE) %>%
  arrange(sp, date, time) %>%
  group_by(sp, date, time, location) %>%
  summarise(lwp.min = mean(lwp.min, na.rm = TRUE)) %>%
  ungroup(sp, date, time, location) %>%
  group_by(sp, date, location) %>%
  arrange(sp, date, time) %>%
  mutate(lwp.diff = lwp.min - lag(lwp.min)) %>%
  ungroup(sp, date, location) %>%
  arrange(sp, date, time) %>%
  left_join(deci, by = "sp") %>%
  transform(deciduousness = factor(deciduousness,
                                   levels = c("Evergreen", "Brevideciduous",
                                              "Facultative Deciduous", "Obligate Deciduous"), ordered = TRUE)) %>%
  unite("sp_date", sp, date, remove = FALSE) %>%
  unite("deci_sp", deciduous, sp, remove = FALSE)
View(lwp.diff)
lwp.min.wide <- lwp.min.diurnal %>% pivot_wider(names_from = time, values_from = c(lwp.min, lwp.se, lwp.diff)) %>%
  mutate(lwp.diff = lwp.diff_Diurnal) %>% select(-lwp.diff_Diurnal, -lwp.diff_Predawn)

lwp.min <- lwp.min.diurnal %>%
  left_join(deci, by = "sp") %>%
  transform(deciduousness = factor(deciduousness,
                                   levels = c("Evergreen", "Brevideciduous",
                                              "Facultative Deciduous", "Obligate Deciduous"), ordered = TRUE))

# https://stackoverflow.com/questions/46971945/how-can-i-have-a-greek-symbol-within-one-of-my-facet-labels
plot.lwp.base <- ggplot(lwp.min %>% subset(!is.na(deciduousness)), aes(x = time, y = lwp.min)) +
  facet_grid(. ~ location) +
  geom_point(aes(color = deciduousness)) +
  ylab(expression(psi[min])) + xlab("Time") +
  guides(colour = guide_legend(order = 1, title = "Deciduousness")) +
  scale_x_discrete(name="",
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

## Panama rainfall gradient preference------
moist.pref <- read.csv("data-raw/Condit_et_al_2013/TreeCommunityDrySeasonSpeciesResponse.csv")
moist.pref <- moist.pref %>% mutate(sp = paste0(tolower(str_sub(species, 1, 4)), str_sub(genus, 1, 2))) %>%
  rename(moist.pref = Moist, moist.pref.2 = Moist.2) %>% select(sp, moist.pref, moist.pref.2)
which(moist.pref$sp == "pipeco") ## appears twice, so removing one.
moist.pref <- moist.pref[-372,]
hab.swp <- read.csv(file.path("data-raw/sp.plot.hab.swp.csv"))
sp.hab <- moist.pref %>% full_join(hab.swp, by = "sp") %>%
  rename(Panama.moist.pref = moist.pref, Panama.moist.pref.2 = moist.pref.2) %>%
  mutate(Plot.swp.pref = med.swp.reg - mean(med.swp.reg, na.rm = TRUE),
         Plot.swp.ENSO = med.swp.dry - mean(med.swp.dry, na.rm = TRUE)) %>%
  select(-med.swp.reg, -med.swp.dry, -sd.swp.reg, -sd.swp.dry, -Panama.moist.pref.2, -Plot.swp.ENSO)

### Hydraulic traits by Brett Wolfe ---------
hyd <- read.csv("data-raw/traits/HydraulicTraits_BrettWolfe/ht1_20200103.csv") #  # Brett's data
hyd <- hyd %>% select(-genus, -species, -deciduousness, -site) %>%
  rename(LMA = lma_gm2_m, WD = xylem_den_m, TLP = tlp_m,
         p50S = p50, p80S = p80, p88S = p88,
         Fcap_Xylem = cwr_xylem_elbow_origin_slope, CWR_Xylem = cwr_xylem_cwr_at_elbow, Felbow_Xylem = cwr_xylem_elbow,
         Fcap_Bark = cwr_bark_elbow_origin_slope, CWR_Bark = cwr_bark_cwr_at_elbow, Felbow_Bark = cwr_bark_elbow,
         BarkThick = barkthickness10mm) %>%
  ## Given Fcap slope values are -ve but they should be positive
  mutate(sp = tolower(sp), Fcap_Xylem = -Fcap_Xylem, Fcap_Bark = -Fcap_Bark) %>%
  left_join(lwp.min.wide %>% subset(location = "PA-BCI"), by = "sp") %>%
  mutate(HSMTLP.50S =  TLP - p50S, HSM50S = lwp.min_Diurnal - p50S,
         HSMTLP.88S =  TLP - p88S, HSM88S = lwp.min_Diurnal - p88S,
         HSMTLP = lwp.min_Diurnal - TLP,
         HSMFelbow_Xylem = lwp.min_Diurnal - Felbow_Xylem,
         HSMFelbow_Bark = lwp.min_Diurnal - Felbow_Bark,
         CWR_Total = CWR_Xylem + CWR_Bark) %>%
  left_join(deci %>% select(-sp4), by = "sp") %>%
  left_join(sp.hab, by = "sp")
length(unique(hyd$sp)) # 27 sp across BCI, PNM, San Lorenzo


####----Phenology by Wolfe hydraulic traits-----
hyd.long <- hyd %>% select(-DeciLvl) %>%
  select(sp, deciduousness, deciduous, location, TLP, p50S, p88S,
         CWR_Total, Fcap_Xylem, CWR_Xylem, Felbow_Xylem, Fcap_Bark, CWR_Bark, Felbow_Bark, WD, LMA,
         HSM50S, HSM88S, HSMTLP, HSMFelbow_Xylem, HSMFelbow_Bark, HSMTLP.50S, HSMTLP.88S,
         Panama.moist.pref, Plot.swp.pref, lwp.min_Diurnal, lwp.min_Predawn) %>%
  gather(trait, value, -sp, -deciduousness, -deciduous,  -location) %>%
  unite("deci_sp", deciduous, sp, remove = FALSE) %>%
  subset(deciduousness != "NA") %>%
  droplevels() %>%
  transform(deciduousness = factor(deciduousness,
                                   levels = c("Evergreen", "Brevideciduous",
                                              "Facultative Deciduous", "Obligate Deciduous"), ordered = TRUE))
kruskal.list <- list()
for(i in unique(hyd.long$trait)) {
  xx <- hyd.long %>% subset(trait == i)
  kruskal.list[[i]] <- cbind(trait = i, kruskal(xx$value, xx$deciduousness, alpha = 0.1, group=TRUE, p.adj="bonferroni")$groups,
                             deciduousness = rownames(kruskal(xx$value, xx$deciduousness, alpha = 0.1, group=TRUE, p.adj="bonferroni")$groups))
}
hyd.kruskal.labels <- do.call(rbind.data.frame, kruskal.list)
head(hyd.kruskal.labels)
hyd.labels <- hyd.kruskal.labels
unique(hyd.labels$trait)
hyd.labels.data <- hyd.labels %>%
  left_join(hyd.long %>% group_by(trait) %>%
              summarise(value = max(value, na.rm = TRUE)), by = c("trait")) %>%
  subset(deciduousness != "NA") %>%
  droplevels() %>%
  transform(deciduousness = factor(deciduousness,
                                   levels = c("Evergreen", "Brevideciduous",
                                              "Facultative Deciduous", "Obligate Deciduous"), ordered = TRUE),
            trait.plot = factor(trait, levels = c("lwp.min_Predawn", "lwp.min_Diurnal", "TLP",
                                                  "p50S", "p88S",
                                                  "HSMTLP", "HSM50S","HSM88S",
                                                  "HSMTLP.50S", "HSMTLP.88S",
                                                  "CWR_Total", "CWR_Xylem", "CWR_Bark",
                                                  "Felbow_Xylem", "Felbow_Bark", "HSMFelbow_Xylem", "HSMFelbow_Bark",
                                                  "Fcap_Xylem", "Fcap_Bark","WD",
                                                  "Panama.moist.pref", "Plot.swp.pref", "LMA"), ordered = TRUE,
                                labels = c(expression(Psi[Predawn]), expression(Psi[min]), expression(Psi[TLP]),
                                           expression(italic('P')['50, Stem']),  expression(italic('P')['88, Stem']),
                                           expression(Psi[min]*' - '*Psi[TLP]),
                                           expression(Psi[min]*' - '*italic('P')['50, Stem']),
                                           expression(Psi[min]*' - '*italic('P')['88, Stem']),
                                           expression(Psi[TLP]*' - '*italic('P')['50, Stem']),
                                           expression(Psi[TLP]*' - '*italic('P')['88, Stem']),
                                           expression('CWR'['Total']), expression('CWR'['Xylem']), expression('CWR'['Bark']),
                                           expression(italic('F')['Elbow, Xylem']), expression(italic('F')['Elbow, Bark']),
                                           expression(Psi[min]*' - '*italic('F')['Elbow, Xylem']),
                                           expression(Psi[min]*' - '*italic('F')['Elbow, Bark']),
                                           expression(italic('F')['Cap, Xylem']), expression(italic('F')['Cap, Bark']),
                                           expression('WD'[stem]),
                                           "Panama.moist.pref", "Plot.swp.pref", "LMA"))) %>% droplevels()

hyd.long <- hyd.long %>% transform(deciduousness = factor(deciduousness,
                                                          levels = c("Evergreen", "Brevideciduous",
                                                                     "Facultative Deciduous", "Obligate Deciduous"), ordered = TRUE),
                                   trait.plot = factor(trait, levels = c("lwp.min_Predawn", "lwp.min_Diurnal", "TLP",
                                                                         "p50S", "p88S",
                                                                         "HSMTLP", "HSM50S","HSM88S",
                                                                         "HSMTLP.50S", "HSMTLP.88S",
                                                                         "CWR_Total", "CWR_Xylem", "CWR_Bark",
                                                                         "Felbow_Xylem", "Felbow_Bark", "HSMFelbow_Xylem", "HSMFelbow_Bark",
                                                                         "Fcap_Xylem", "Fcap_Bark","WD",
                                                                         "Panama.moist.pref", "Plot.swp.pref", "LMA"), ordered = TRUE,
                                                       labels = c(expression(Psi[Predawn]), expression(Psi[min]), expression(Psi[TLP]),
                                                                  expression(italic('P')['50, Stem']),  expression(italic('P')['88, Stem']),
                                                                  expression(Psi[min]*' - '*Psi[TLP]),
                                                                  expression(Psi[min]*' - '*italic('P')['50, Stem']),
                                                                  expression(Psi[min]*' - '*italic('P')['88, Stem']),
                                                                  expression(Psi[TLP]*' - '*italic('P')['50, Stem']),
                                                                  expression(Psi[TLP]*' - '*italic('P')['88, Stem']),
                                                                  expression('CWR'['Total']), expression('CWR'['Xylem']), expression('CWR'['Bark']),
                                                                  expression(italic('F')['Elbow, Xylem']), expression(italic('F')['Elbow, Bark']),
                                                                  expression(Psi[min]*' - '*italic('F')['Elbow, Xylem']),
                                                                  expression(Psi[min]*' - '*italic('F')['Elbow, Bark']),
                                                                  expression(italic('F')['Cap, Xylem']), expression(italic('F')['Cap, Bark']),
                                                                  expression('WD'[stem]),
                                                                  "Panama.moist.pref", "Plot.swp.pref", "LMA"))) %>% droplevels()

h1 <- ggplot(hyd.labels.data, aes(x = deciduousness, y = value)) +
  geom_boxplot(data = hyd.long, aes(fill = deciduousness), stat = "boxplot", notch = TRUE) +
  geom_jitter(data = hyd.long, width = 0.05, shape = 21, fill = "darkgray", color = "black", show.legend = FALSE, alpha = 0.7) +
  facet_wrap(. ~  trait.plot, scales = "free_y", labeller = label_parsed) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 1, hjust = 1)) +
  scale_color_brewer(palette = "Greens", direction = -1) +
  scale_fill_brewer(name = "Deciduousness", palette = "Greens", direction = -1) +
  xlab("Deciduousness") + ylab("Value")
ggsave(file.path(figures.folder, paste0("BrettWolfe_traits_vs_deciduousness.jpeg")),
       plot = h1, height = 7, width = 10, units ='in')
h1.1 <- h1 + geom_text(aes(label = groups), vjust = 1, hjust = 0, show.legend = FALSE)
ggsave(file.path(figures.folder, paste0("BrettWolfe_traits_vs_deciduousness_kruskal.labels.jpeg")),
       plot = h1.1, height = 7, width = 10, units ='in')

select.traits <- c("lwp.min_Predawn", "lwp.min_Diurnal", "TLP", "p50S", "p88S",
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

###

## BCI traits
bci.traits <- read.csv("data-raw/traits/BCITRAITS_20101220.csv") %>%
  rename(form1 = GRWFRM1., sp = SP., WSG_Chave = WSG_CHAVE) %>% mutate(sp = tolower(sp))

## Nobert Kunert traits --------
traits.indi <- read.csv("data-raw/traits/HydraulicTraits_Kunert/hydraulic_traits_panama_kunert.csv") # Nobby's data
traits <- traits.indi %>%
  rename(TLP = mean_TLP_Mpa, Chl = Chl_m2_per_g) %>%
  select(sp, TLP, Chl) %>% group_by(sp) %>%
  summarise_all(mean, na.rm = TRUE)
traits <- traits %>%
  left_join(lwp.min.wide %>% select(sp, lwp.min_Predawn, lwp.min_Diurnal) %>%
              subset(location = "PA-BCI"), by = "sp")

leaf_cond.models <- read.csv("data-raw/traits/HydraulicTraits_Kunert/Panama_fits_leaf_K_p50_Kunert.csv")
leaf.k.p80 <- leaf_cond.models %>% subset(model_type == "Exponential") %>%
  mutate(sp = data.type, p50L = -psi_kl50, p80L = -psi_kl80, KmaxL = Kmax) %>% # these are Kmax that are extrapolated from the exponential curve
  select(sp, KmaxL, p50L, p80L) # 21 species
traits <- traits %>%
  full_join(leaf.k.p80 %>% mutate(sp = as.character(sp)), by = "sp") %>%
  mutate(HSMTLP.50L = TLP - p50L, HSMLWP.50L = lwp.min_Diurnal - p50L,
         HSMTLP.80L = TLP - p80L, HSMLWP.80L = lwp.min_Diurnal - p80L,
         HSMLWP.TLP = lwp.min_Diurnal - TLP) %>%
  left_join(sp.hab, by = "sp") %>%
  left_join(deci %>% select(-sp4), by = "sp") %>%
  left_join(bci.traits %>% select(sp, form1, WSG_Chave), by = "sp")

# > with(traits, table(deciduousness))
# Evergreen    Obligate Deciduous Facultative Deciduous        Brevideciduous
# 26                     3                    11                     8
####----Phenology by Kunert hydraulic traits-----

traits.long <- traits %>% select(-DeciLvl) %>%
  gather(trait, value, -sp, -deciduousness, -deciduous, -form1) %>%
  subset(deciduousness != "NA") %>%
  droplevels() %>%
  transform(deciduousness = factor(deciduousness,
                                   levels = c("Evergreen", "Brevideciduous",
                                              "Facultative Deciduous", "Obligate Deciduous"), ordered = TRUE))
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
              summarise(value = max(value, na.rm = TRUE)), by = c("trait")) %>%
  subset(deciduousness != "NA") %>%
  droplevels() %>%
  transform(deciduousness = factor(deciduousness,
                                   levels = c("Evergreen", "Brevideciduous",
                                              "Facultative Deciduous", "Obligate Deciduous"), ordered = TRUE),
            trait.plot = factor(trait, levels = c("KmaxL", "lwp.min_Predawn", "lwp.min_Diurnal", "TLP", "p50L", "p80L",
                                                  "HSMLWP.TLP", "HSMLWP.50L", "HSMTLP.50L",
                                                  "HSMLWP.80L", "HSMTLP.80L",
                                                  "Panama.moist.pref", "Plot.swp.pref", "WSG_Chave", "Chl"), ordered = TRUE,
                                labels = c(expression(italic(K)[max]), expression(Psi[predawn]), expression(Psi[min]),
                                           expression(Psi[TLP]), expression(italic('P')['50, Leaf']), expression(italic('P')['80, Leaf']),
                                           expression(Psi[min]*' - '*Psi[TLP]),
                                           expression(Psi[min]*' - '*italic('P')['50, Leaf']),
                                           expression(Psi[TLP]*' - '*italic('P')['50, Leaf']),
                                           expression(Psi[min]*' - '*italic('P')['80, Leaf']),
                                           expression(Psi[TLP]*' - '*italic('P')['80, Leaf']),
                                           "Panama.moist.pref", "Plot.swp.pref", "WSG_Chave", "Chl"))) %>% droplevels()
traits.long <- traits.long %>% transform(
  trait.plot = factor(trait, levels = c("KmaxL", "lwp.min_Predawn", "lwp.min_Diurnal", "TLP", "p50L", "p80L",
                                        "HSMLWP.TLP", "HSMLWP.50L", "HSMTLP.50L",
                                        "HSMLWP.80L", "HSMTLP.80L",
                                        "Panama.moist.pref", "Plot.swp.pref", "WSG_Chave", "Chl"), ordered = TRUE,
                      labels = c(expression(italic(K)[max]), expression(Psi[predawn]), expression(Psi[min]),
                                 expression(Psi[TLP]), expression(italic('P')['50, Leaf']), expression(italic('P')['80, Leaf']),
                                 expression(Psi[min]*' - '*Psi[TLP]),
                                 expression(Psi[min]*' - '*italic('P')['50, Leaf']),
                                 expression(Psi[TLP]*' - '*italic('P')['50, Leaf']),
                                 expression(Psi[min]*' - '*italic('P')['80, Leaf']),
                                 expression(Psi[TLP]*' - '*italic('P')['80, Leaf']),
                                 "Panama.moist.pref", "Plot.swp.pref", "WSG_Chave", "Chl"))) %>% droplevels()

ggplot(traits.labels.data, aes(x = deciduousness, y = value)) +
  facet_wrap(. ~  trait.plot, scales = "free_y", labeller = label_parsed) +
  geom_text(aes(label = groups), vjust = 1, hjust = 0, show.legend = FALSE) +
  geom_boxplot(data = traits.long, aes(fill = deciduousness), stat = "boxplot", notch = TRUE) +
  geom_jitter(data = traits.long, width = 0.05, shape = 21, fill = "darkgray",
              color = "black", show.legend = FALSE, alpha = 0.7) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 1, hjust = 1)) +
  scale_color_brewer(palette = "Greens", direction = -1) +
  scale_fill_brewer(name = "Deciduousness", palette = "Greens", direction = -1)
ggsave(file.path(figures.folder, paste0("kunert_traits_vs_deciduousness.jpeg")), height = 9, width = 12, units ='in')

## species wise for sp in hyd.traits

select.traits <- c("lwp.min_Predawn", "lwp.min_Diurnal", "TLP", "p50S", "p88S",
                   "HSMTLP", "HSM50S","HSM88S", "HSMTLP.50S", "HSMTLP.88S")
traits.long <- traits.long %>%
  unite("deci_sp", deciduous, sp, remove = FALSE) %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(deciduousness)]), ordered=TRUE),
         deci_sp.plot = factor(deci_sp, levels=unique(deci_sp[order(deciduousness)]), ordered=TRUE))
# just for sp with hyd.traits, but traits.long does not have all those sp, and hab preference and WSG traits will be missed
## so beginning with those other traits
traits.wide <- traits.long %>% pivot_wider(names_from = trait, values_from = value, -trait.plot)
traits.long.hyd <- sp.hab %>%
  full_join(bci.traits %>% select(sp, form1, WSG_Chave), by = "sp") %>%
  full_join(deci %>% select(-sp4), by = "sp") %>%
  subset(sp %in% hyd$sp) %>%
  left_join(traits.wide %>% select(-form1, -deciduous, -deciduousness,
                                   -WSG_Chave, -Panama.moist.pref, -Plot.swp.pref,
                                   -deci_sp, -sp.plot, -deci_sp.plot), by = "sp") %>%
  pivot_longer(cols = c(-sp, -form1, -deciduous, - deciduousness, - DeciLvl),
               names_to = "trait", values_to = "value") %>%
  unite("deci_sp", deciduous, sp, remove = FALSE) %>%
  transform(deciduousness = factor(deciduousness,
                                   levels = c("Evergreen", "Brevideciduous",
                                              "Facultative Deciduous", "Obligate Deciduous"), ordered = TRUE),
            trait.plot = factor(trait, levels = c("KmaxL", "lwp.min_Predawn", "lwp.min_Diurnal", "TLP", "p50L", "p80L",
                                                  "HSMLWP.TLP", "HSMLWP.50L", "HSMTLP.50L",
                                                  "HSMLWP.80L", "HSMTLP.80L",
                                                  "Panama.moist.pref", "Plot.swp.pref", "WSG_Chave", "Chl"), ordered = TRUE,
                                labels = c(expression(italic(K)[max]), expression(Psi[predawn]), expression(Psi[min]),
                                           expression(Psi[TLP]), expression(italic('P')['50, Leaf']), expression(italic('P')['80, Leaf']),
                                           expression(Psi[min]*' - '*Psi[TLP]),
                                           expression(Psi[min]*' - '*italic('P')['50, Leaf']),
                                           expression(Psi[TLP]*' - '*italic('P')['50, Leaf']),
                                           expression(Psi[min]*' - '*italic('P')['80, Leaf']),
                                           expression(Psi[TLP]*' - '*italic('P')['80, Leaf']),
                                           "Panama.moist.pref", "Plot.swp.pref", "WSG_Chave", "Chl"))) %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(deciduousness)]), ordered=TRUE),
         deci_sp.plot = factor(deci_sp, levels=unique(deci_sp[order(deciduousness)]), ordered=TRUE)) %>%
  droplevels()

t2 <- ggplot(traits.long.hyd) +
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

## Demographic data----
load("results/demo.sp.RData")
load("results/demo.sp_size.RData")
load("results/mrate.long.RData")

census.years <- c(1982, 1985, 1990, 1995, 2000, 2005, 2010, 2015)

size.level_key <- c(tiny = "Tiny (1-5cm)", small = "Small (5-10cm)",
                    medium = "Medium (10-30cm)", large = "Large (>= 30cm)")
census.years.short <- format(strptime(census.years, "%Y"), "%y")

mrate.long <- mrate.long %>%
  separate(sp_size, c("sp", "size", sep = "_"), remove = FALSE, extra = "drop", fill = "right") %>%
  select(-"_") %>%
  left_join(deci, by = "sp") %>%
  mutate(censusint.m = recode(census, `1985` = "1982-85", `1990` = "1985-90", `1995` = "1990-95", `2000` = "1995-00", `2005` = "2000-05", `2010` = "2005-10", `2015` = "2010-15"),
         size.num = recode_factor(size, !!!size.level_key),
         size = factor(size, levels = c("tiny", "small", "medium", "large"))) %>%
  left_join(demo.sp_size %>% mutate(mean.mrate = mrate, mean.grate = grate) %>%
              select(sp_size, mean.mrate, mean.grate), by = "sp_size") %>%
  mutate(diff.mrate = mrate - mean.mrate) %>%
  transform(deciduousness = factor(deciduousness,
                                   levels = c("Evergreen", "Brevideciduous",
                                              "Facultative Deciduous", "Obligate Deciduous"), ordered = TRUE))

## Mortality vs growth rate by phenology-----

growth.type <- "med"
formula = y ~ x
demo.sp_size <- demo.sp_size %>%
  left_join(deci, by = "sp")
p.d0 <- ggplot(demo.sp_size %>% subset(size != "NA"),
               aes(y = mrate, x = grate)) +
  geom_point() +
  facet_wrap(. ~ size, scales = "free_y") +
  geom_smooth(method = "lm") +
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

### Is leaf phenology linked to vulnerability to different drought intensity and duration?------

m1 <- ggplot(mrate.long %>%
         subset(!is.na(size) & avg.abund >= n.threshold & !is.na(deciduousness)),
       aes(x = deciduousness, y = diff.mrate, color = avg.abund)) +
  scale_color_gradient(name = "Mean\nAbundance", trans = "rev_sqrt",
                       low = "red", high = "blue", breaks = c(100, 1000, 5000, 10000, 20000, 30000)) +
  facet_grid(size.num ~ censusint.m, scales = "free_y") +
  geom_hline(aes(yintercept = 0), color = "blue", size = 0.5) +
  geom_boxplot(aes(fill = deciduousness), stat = "boxplot", notch = TRUE) +
  scale_fill_brewer(name = "", palette = "Greens", direction = -1) +
  theme(legend.position = "top") +
  ylab(expression("Mortality Rate - Mean (% per year)")) + xlab("Deciduousness") +
  ggtitle("Drought Mortality by Leaf Phenology") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 1, hjust = 1))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_deci_by_size_aboveN",
                        n.threshold, ".jpeg")), plot = m1, height = 8, width = 9, units='in')
m2 <- m1 %+% subset(mrate.long, size == "large" & avg.abund >= n.threshold & !is.na(deciduousness))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_deci_by_size_aboveN",
                        n.threshold, "_large.jpeg")), plot = m2, height = 5, width = 9, units='in')
m2 <- m1 %+% subset(mrate.long, size == "large" & avg.abund >= n.threshold & !is.na(deciduousness)) +
  geom_jitter(width = 0.05, shape = 21, fill = "darkgray", color = "black", show.legend = FALSE, alpha = 0.7)
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_deci_by_size_aboveN",
                        n.threshold, "_large_points.jpeg")), plot = m2, height = 5, width = 9, units='in')

hyd.wide <- hyd.long %>% pivot_wider(names_from = trait, values_from = value, -trait.plot)
mrate.long.hyd <- subset(mrate.long, sp %in% hyd$sp) %>%
  left_join(hyd.wide %>% select(sp, p88S, HSMTLP.88S, HSM88S), by = "sp") %>%
  mutate(sp.plot = factor(sp, levels=unique(sp[order(HSMTLP.88S)]), ordered=TRUE))
# show_col(viridis_pal()(4))
m3.base <- ggplot(subset(mrate.long.hyd, size == "large" & sp %in% c("cordal", "luehse", "tab1ro"))) +
  facet_grid(size.num ~ censusint.m, scales = "free_y") +
  theme(legend.position = "top") +
  ylab(expression("Mortality Rate - Mean (% per year)")) + xlab("Deciduousness") +
  guides(fill = guide_legend(title = "Deciduousness")) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 0.5, hjust = 1))
m3.1 <- m3.base + geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Species leafless in early wet season, with increasing HSM '*Psi[TLP]*' - '*italic('P')['88,Stem']))
  #  does not work: fill = "Facultative Deciduous"
  # guides(fill = guide_legend(title = "Deciduousness",
  #                            override.aes = list(fill = c("Facultative Deciduous" = "#35B779FF"))))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_deci_HSMTLP.88S_increasing_spp_leafless in early wet season.jpeg")),
       plot = m3.1, height = 4, width = 9, units='in')

m4.1 <- m3.base %+% subset(mrate.long.hyd, sp %in% hyd$sp & size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Mortality for Evergreen Species with increasing HSM '*Psi[TLP]*' - '*italic('P')['88,Stem']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_HSMTLP.88S_increasing_evergreens.jpeg")),
       plot = m4.1, height = 4, width = 9, units='in')

mrate.long.hyd <- mrate.long.hyd %>% mutate(sp.plot = factor(sp, levels=unique(sp[order(-p88S)]), ordered=TRUE))
m3.2 <- m3.base + geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Species leafless in early wet season, with increasingly more negative '*italic('P')['88,Stem']))
ggsave(file.path(paste0(figures.folder, "/sp_Mortality_rate_by_period_deci_p88S_increasing_spp_leafless in early wet season.jpeg")),
       plot = m3.2, height = 4, width = 9, units='in')

m4.2 <- m3.base %+% subset(mrate.long.hyd, sp %in% hyd$sp & size == "large" & deciduous == "E") +
  geom_col(aes(x = sp.plot, y = diff.mrate, fill = deciduousness)) +
  ggtitle(expression('Mortality for Evergreen Species with increasingly more negative '*italic('P')['88,Stem']))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_p88S_increasing_evergreens.jpeg")),
       plot = m4.2, height = 4, width = 9, units='in')

## all evergreens
mrate.long <- mrate.long %>% mutate(sp.plot = factor(sp, levels=unique(sp[order(mean.mrate)]), ordered=TRUE))
m5 <- m3.base %+% subset(mrate.long, size == "large" & deciduous == "E") +
  geom_col(aes(x = sp, y = diff.mrate, fill = deciduousness)) +
  ggtitle("Mortality for Evergreen Species") +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 0.5, hjust = 1, size = 4))
ggsave(file.path(paste0(figures.folder,
                        "/sp_Mortality_rate_by_period_evergreens.jpeg")),
       plot = m5, height = 4, width = 15, units='in')

## Plot addtional mortality by period mean swp-------

### psi

census.meds <- readr::read_rds("results/census.mediandates.rds")
census.beg <- census.meds[3: length(census.meds)]
cut.breaks <- census.beg
cut.breaks.2 <- as.Date(paste0(seq(1990, 2015, by = 5), "-01-01"))
cut.labels.2 <- paste0(seq(1990, 2010, by = 5), "-", seq(1995, 2015, by = 5))

load(file.path("data-raw/psi.rda"))

## by depth panels

psi <- psi %>%
  mutate(interval.yrs = cut(date, include.lowest = TRUE, breaks = cut.breaks,
                            labels = cut.labels.2, right = TRUE))
psi.stat.1 <- psi %>%
  group_by(interval.yrs, date, depth) %>%
  summarise(mean = -mean(psi, na.rm = TRUE),
            upper.CI = -quantile(psi, probs = 0.975),
            lower.CI = -quantile(psi, probs = 0.025))
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
  geom_hline(yintercept = c(0.5, 1, 1.5, 2), size = 0.2) +
  geom_line(aes(x = date, y = mean, group = depth, color = depth), size = 0.3) +
  scale_color_continuous(trans = "reverse", guide = guide_colorbar(title = "Depth\n(cm)")) +
  geom_vline(xintercept = census.beg) + # as.Date(paste0(c(1990:2015), "-01-01")) +
  scale_y_reverse(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 10, 15), limits = c(2.5, 0)) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%Y")) +
  coord_cartesian(xlim = c(as.Date("1990-01-01"), as.Date("2018-12-31"))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  # facet_grid(interval.yrs ~ .) +
  ylab("-Soil Water Potential [MPa]") +
  theme(text = element_text(size = 12)) + xlab("")
  # ggtitle("PSI:Violins of interval means of best-fit ensembles")
ggsave("psi_model_daily_bestfit_params.top.few_CI_full.jpeg", plot = plot.psi.stat.1,
       file.path(figures.folder), device = "jpeg", height = 3, width = 20, units='in')

psi.2 <- psi %>%
  mutate(interval.yrs.to.plot = cut(date, include.lowest = TRUE, breaks = cut.breaks.2,
                                    labels = cut.labels.2, right = TRUE))

psi.stat.2 <- psi.2 %>%
  group_by(interval.yrs.to.plot, date, depth) %>%
  summarise(mean = mean(psi, na.rm = TRUE),
            upper.CI = quantile(psi, probs = 0.975),
            lower.CI = quantile(psi, probs = 0.025))
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
  geom_line(aes(x = days, y = mean, group = depth, color = depth), size = 0.3) +
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

psi.stat.3 <- psi %>%
  mutate(doy = format(date, "%j")) %>%
  group_by(doy, depth) %>%
  summarise(mean = mean(psi, na.rm = TRUE),
            sd = sd(psi, na.rm = TRUE),
            upper.CI = quantile(psi, probs = 0.975),
            lower.CI = quantile(psi, probs = 0.025)) %>%
  ungroup(doy, depth) %>%
  mutate(doy = as.numeric(doy))
psi.stat.4 <- psi %>%
  mutate(doy = format(date, "%j"),
         year = format(date, "%Y")) %>%
  group_by(doy, year, depth) %>%
  summarise(mean = mean(psi, na.rm = TRUE),
            median = median(psi, na.rm = TRUE),
            sd = sd(psi, na.rm = TRUE),
            upper.CI = quantile(psi, probs = 0.975),
            lower.CI = quantile(psi, probs = 0.025)) %>%
  ungroup(doy, year, depth) %>%
  mutate(doy = as.numeric(doy))
psi.stat.5 <- psi.stat.4 %>%
  group_by(doy, depth) %>%
  summarise(mean.clim = mean(mean, na.rm = TRUE),
            sd = sd(mean, na.rm = TRUE),
            upper.CI = quantile(mean, probs = 0.975),
            lower.CI = quantile(mean, probs = 0.025)) %>%
  ungroup(doy, depth) %>%
  mutate(doy = as.numeric(doy))
rectangles.3 <- data.frame(
  xmin = 120,
  xmax = 335,
  ymin = 0,
  ymax = -2.5
)
source("code/Utilities/plot.ticks.R")
plot.psi.stat.5.base <- ggplot(psi.stat.5 %>% droplevels()) +
  # geom_rect(data=rectangles.3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
  #           fill='gray80', alpha=0.8) +
  scale_x_continuous(breaks = c(seq(0, 360, by = 60))) +
  coord_cartesian(ylim = c(-2.5, 0)) +
  theme(panel.grid.major.y = element_line()) +
  ylab("Soil Water Potential (MPa)") + xlab("Day of the Year") +
  theme(text = element_text(size = 12))
plot.psi.stat.5 <- plot.psi.stat.5.base +
  geom_ribbon(aes(x = doy, ymin = lower.CI, ymax = upper.CI), alpha = 0.3, fill = "grey20") +
  geom_line(aes(x = doy, y = mean.clim, group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
  guides(color = guide_legend(title = "Depth(m)", order = 2, override.aes = list(size = 3)))
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_censuspanels_climatology.jpeg",
       plot = plot.ticks(plot.psi.stat.5),
       file.path(figures.folder), device = "jpeg", height = 4, width = 7, units='in')

plot.psi.stat.5.over <- plot.psi.stat.5.base %+% subset(psi.stat.5, depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)) +
  geom_line(aes(x = doy, y = mean.clim, linetype = "Climatalogy", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
  geom_line(data = subset(psi.stat.4, year == 2016 & depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)),
            aes(x = doy, y = mean, linetype = "2016", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
  guides(color = guide_legend(title = "Depth(m)", legend.position = "right", order = 1, override.aes = list(size = 3)),
         linetype = guide_legend(order = 2, title = NULL, legend.position = "top", override.aes =
                                   list(linetype = c("Climatalogy" = "dashed", "2016" = "solid")))) +
  coord_cartesian(ylim = c(-3, 0)) + ggtitle("2016")
ggsave("psi_model_daily_bestfit_params.top.few_CI_full_censuspanels_climatology_over.jpeg",
       plot = plot.ticks(plot.psi.stat.5.over),
       file.path(figures.folder), device = "jpeg", height = 4, width = 7, units='in')

pdf(paste0(figures.folder, "/psi_model_daily_bestfit_params.top.few_CI_full_censuspanels_climatology_over_by_year.pdf"), height = 4, width = 7)
for (i in unique(psi.stat.4$year)) {
  plot.psi.stat.5.yr <- plot.psi.stat.5.base %+% subset(psi.stat.5, depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)) +
    geom_line(aes(x = doy, y = mean.clim, linetype = "Climatalogy", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
    geom_line(data = subset(psi.stat.4, year == i & depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)),
              aes(x = doy, y = mean, linetype = "Year", group = as.factor(depth), color = as.factor(depth)), size = 0.5) +
    guides(color = guide_legend(title = "Depth(m)", order = 1, override.aes = list(size = 3)),
           linetype = guide_legend(order = 2, title = NULL, override.aes =
                                     list(linetype = c("Climatalogy" = "solid", "Year" = "dashed")))) +
    coord_cartesian(ylim = c(-3, 0), xlim = c(0, 200)) + ggtitle(i)
  print(plot.psi.stat.5.yr)
}
dev.off()

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
    geom_line(aes(x = doy, y = mean.clim, linetype = "Climatalogy", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
    geom_line(data = subset(psi.stat.4, year == year.on & depth %in% c(0.06, 0.12, 0.37, 0.62, 1, 1.7, 2.9)),
              aes(x = doy, y = mean, linetype = "Year", group = as.factor(depth), color = as.factor(depth)), size = 0.3) +
    guides(color = guide_legend(title = "Depth(m)", order = 1, override.aes = list(size = 3)),
           linetype = guide_legend(order = 2, title = NULL, override.aes =
                                     list(linetype = c("Climatalogy" = "solid", "Year" = "dashed")))) +
    coord_cartesian(ylim = c(-3, 0), xlim = c(0, 200)) +
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

ggsave("arrange.pdf", arrangeGrob(grobs = plot.list.grid, ncol = 5), width = 15, height = 10, units = "in")
ggsave("arrange.pdf", arrangeGrob(grobs = plot.list.grid, ncol = 5), width = 15, height = 10, units = "in")

