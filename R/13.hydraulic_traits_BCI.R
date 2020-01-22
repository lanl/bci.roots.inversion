##-------------------------
## Date : 22 Jul 2018
## Author : Rutuja
## Title : Uptake depth and hydrualic traits
##-------------------------
## Do species with shallow uptake depth have greater drought tolerance? (more -ve TLP)

rm(list = ls())
gc()
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(tidyverse, ggcorrplot, Hmisc, PerformanceAnalytics, corrplot, devtools)
install_github("vqv/ggbiplot")
  # graphics info
theme_set(theme_bw())
theme_update(text = element_text(size=14),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)
## adding uptake depth index
load("results/GLUEsetup_part1_BCI.RData") # has model info and data on obs
load("results/4.1GLUEsetup_part2_BCI.RData") # has n.ensembles and growth and si matrix
load("results/demo.sp_size.RData")
intervals <- info$intervals
n.ensembles <- growth_by_si.info$n.ensembles
growth.type <- growth_by_si.info$growth.type
si.type <- growth_by_si.info$si.type
goodness.fit <- 0.3
soil.depths <- unique(info$root.param.long$depth)
##
load(file = paste("results/splevel/ds.bestfit_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", intervals, "_id.Rdata", sep = ""))
ds <- ds.bestfit
load(file = paste("results/commlevel/ds.bestfit_cor", goodness.fit, "_", si.type, "_", n.ensembles, "_", growth.type, "_", intervals, "_id.Rdata", sep = ""))
ds <- rbind(ds, ds.bestfit) %>% mutate(tlplevel = as.factor(tlplevel)) %>%
  subset(!is.na(udi)) %>%
  mutate(size = factor(size, levels = c("tiny", "small", "medium", "large"))) %>%
  droplevels()
hyd <- read.csv("data-raw/traits/ht1_20200103.csv") #  # Brett's data
hyd <- hyd %>% select(-genus, -species, -deciduousness, -site) %>% mutate(sp = tolower(sp))
traits.indi <- read.csv("data-raw/traits/hydraulic_traits_panama_kunert.csv") # Nobby's data
traits <- traits.indi %>% group_by(sp) %>%
  select(-idividual, -ind_ID,
         -leaf_area_dry_cm2, -leaf_area_fresh_cm2,
         -sum_dry_mass_leaf_blade_g, -sum_fresh_mass_leaf_blade_g, -PLA_dry_percent) %>%
  summarise_all(mean, na.rm = TRUE)
colnames(traits) <- c("sp", "LA", "LMA", "LT", "LD", "WMA", "SPAD", "Chl", "TLP", "WD")
tlplevels <- c("sp", "comm")
size.class <- levels(ds$size)

for (j in 1: length(tlplevels)){
  for (i in 1: length(size.class)){
    ds.sp <- ds %>% select(sp, size, udi, tlplevel) %>% subset(tlplevel == tlplevels[j] & size == size.class[i])
    dst <- left_join(traits, ds.sp, by = "sp") %>% rename(UDI = udi)
    df <- dplyr::select_if(dst, is.numeric)
    cor.mat <- cor(df, use = "complete.obs")
    ggcorrplot(cor.mat,
               hc.order = FALSE,
               type = "lower",
               lab = TRUE,
               title = paste0("TLPlevel = ", tlplevels[j], "\n Tree Size Class = ", size.class[i]))
    ggsave(file.path(paste0("figures/traits/", growth.type, "/sp_UDI_vs_traits_corrplot_tlp", tlplevels[j], "_", size.class[i], ".jpeg")), height = 10, width = 10, units ='in')
    }
}


for (j in 1: length(tlplevels)){
  for (i in 1: length(size.class)){
    ds.sp <- ds %>% select(sp, size, udi, tlplevel) %>% subset(tlplevel == tlplevels[j] & size == size.class[i])
    dst <- left_join(traits, ds.sp, by = "sp") %>%
      left_join(demo.sp_size %>%
                  subset(size == size.class[i]) %>%
                  select(sp, growth.mean, mrate), by = "sp") %>%
      rename(UDI = udi, MortRate = mrate, GrowthRate = growth.mean)
    df <- dplyr::select_if(dst, is.numeric)
    cor.mat <- cor(df, use = "complete.obs")
    ggcorrplot(cor.mat,
               hc.order = FALSE,
               type = "lower",
               lab = TRUE,
               title = paste0("TLPlevel = ", tlplevels[j], "\n Tree Size Class = ", size.class[i]))
    ggsave(file.path(paste0("figures/traits/", growth.type, "/sp_UDI_vs_traits_corrplot_tlp", tlplevels[j], "_", size.class[i], "_demo.jpeg")), height = 10, width = 10, units ='in')
  }
}


for (j in 1: length(tlplevels)){
  for (i in 1: length(size.class)){
    ds.sp <- ds %>% select(sp, size, udi, tlplevel) %>% subset(tlplevel == tlplevels[j] & size == size.class[i])
    dst <- left_join(traits, ds.sp, by = "sp") %>%
      left_join(demo.sp_size %>%
                  subset(size == size.class[i]) %>%
                  select(sp, growth.mean, mrate), by = "sp") %>%
      rename(UDI = udi, MortRate = mrate, GrowthRate = growth.mean)
    df <- dplyr::select_if(dst, is.numeric)
    cor.mat <- cor(df, use = "complete.obs")
    pdf(file.path(paste0("figures/traits/", growth.type, "/chart/sp_UDI_vs_traits_corrplot_tlp",
                                     tlplevels[j], "_", size.class[i], "_demo_chart.pdf")))
    chart.Correlation(cor.mat, histogram=TRUE, pch=19)# does not take title
    dev.off()
  }
}

data.subset <- c("sub", "full")
k <- 1 #for data.subset[k]
## hyd traits
for (j in 1: length(tlplevels)){
  for (i in 1: length(size.class)){
    if (k == 1) {
      hydt <- hyd %>% select(sp, p50, tlp_m, barkthickness10mm, ldmc_m, cwr_xylem_cwr_at_elbow)
    }
    ds.sp <- ds %>% select(sp, size, udi, tlplevel) %>% subset(tlplevel == tlplevels[j] & size == size.class[i])
    dst <- left_join(hydt, ds.sp, by = "sp") %>%
      left_join(demo.sp_size %>%
                  subset(size == size.class[i]) %>%
                  select(sp, growth.mean, mrate), by = "sp") %>%
      rename(UDI = udi, MortRate = mrate, GrowthRate = growth.mean)
    df <- dplyr::select_if(dst, is.numeric)
    cor.mat <- cor(df, use = "complete.obs")
    ggcorrplot(cor.mat,
               hc.order = TRUE,
               type = "lower",
               lab = TRUE,
               title = paste0("TLPlevel = ", tlplevels[j], "\n Tree Size Class = ", size.class[i]))
    ggsave(file.path(paste0("figures/traits/", growth.type, "/hyd/", data.subset[k] ,"/sp_UDI_vs_traits_corrplot_tlp",
                            tlplevels[j], "_", size.class[i],"_", data.subset[k] ,"_demo.jpeg")),
           height = 10*k, width = 10*k, units ='in')
  }
}
### PCA

## knobby's data
for (j in 1: length(tlplevels)){
  for (i in 1: length(size.class)){
    ds.sp <- ds %>% select(sp, size, udi, tlplevel) %>%
      subset(tlplevel == tlplevels[j] & size == size.class[i])
    dst <- left_join(traits %>% mutate(sp = as.character(sp)), ds.sp, by = "sp") %>%
      left_join(demo.sp_size %>%
                  subset(size == size.class[i] & tlplevel == tlplevels[j]) %>%
                  select(sp, growth.mean, mrate), by = "sp") %>%
      dplyr::rename(UDI = udi, MortRate = mrate, GrowthRate = growth.mean)

    df.pca <- dst %>% as_tibble() %>% column_to_rownames(var = "sp") %>% select(-tlplevel, -size)
    df.pca <- df.pca[complete.cases(df.pca),]
    pca <- prcomp(df.pca, center = TRUE, scale. = TRUE)
    tiff(file.path(paste0("figures/traits/", growth.type, "/pca/sp_UDI_traits_demo",
                         tlplevels[j], "_", size.class[i], ".tiff")))
    ggbiplot(pca)
    dev.off()
    tiff(file.path(paste0("figures/traits/", growth.type, "/pca/sp_UDI_traits_demo",
                          tlplevels[j], "_", size.class[i], "_names.tiff")))
    ggbiplot(pca, labels = rownames(df.pca))
    dev.off()
  }
}

## brett's data

for (j in 1: length(tlplevels)){
  for (i in 1: length(size.class)){
    ds.sp <- ds %>% select(sp, size, udi, tlplevel) %>%
      subset(tlplevel == tlplevels[j] & size == size.class[i])
    dst <- left_join(hyd %>% mutate(sp = as.character(sp)), ds.sp, by = "sp") %>%
      left_join(demo.sp_size %>%
                  subset(size == size.class[i] & tlplevel == tlplevels[j]) %>%
                  select(sp, growth.mean, mrate), by = "sp") %>%
      dplyr::rename(UDI = udi, MortRate = mrate, GrowthRate = growth.mean)

    df.pca <- dst %>% as_tibble() %>% column_to_rownames(var = "sp") %>% select(-tlplevel, -size)
    df.pca <- df.pca[complete.cases(df.pca),]
    pca <- prcomp(df.pca, center = TRUE, scale. = TRUE)
    tiff(file.path(paste0("figures/traits/", growth.type, "/hyd/pca/sp_UDI_traits_demo",
                          tlplevels[j], "_", size.class[i], ".tiff")))
    ggbiplot(pca)
    dev.off()
    tiff(file.path(paste0("figures/traits/", growth.type, "hyd//pca/sp_UDI_traits_demo",
                          tlplevels[j], "_", size.class[i], "_names.tiff")))
    ggbiplot(pca, labels = rownames(df.pca))
    dev.off()
  }
}
