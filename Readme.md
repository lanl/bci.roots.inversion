Workflow
================
Rutuja Chitra-Tarak
10 November, 2020

Pre-requisite data
==================

**Growth data**

Calculate species-size wise growth time series from forest inventory data. Calculate growth residuals after removing effect of size.

``` r
source("code/01.0_Preparing_stem_growth_time_series.R")
```

**Mortality data**

``` r
source("code/02.0_demographic_data_BCI.R")
```

#### All pre-requisites

Gather all pre-requisites for rooting depth inverse model (effective rooting depth ) in a place and traits data for post-processing

-   VPD-GPP relationship
-   Leaf vulnerability curves
-   LAI seasonality
-   Hydraulic traits

``` r
source("code/05.0_Prepare_data_for_correlative analyses.R")
```

Script to load data generated from the data preparation script above for use in inverse model as well as manuscript

``` r
source("code/06.0_load.R")
```

Model Setup
===========

Set-up and run effective rooting depth models, evaluate model, and join ERD from the predicted model with hydraulic traits and mortality data, and calcualted realized hydraulic risk

``` r
source("code/07.0_PhenoDemoTraitsPsi.R")
```

Manuscript
==========

**Figures**

``` r
source("code/08.0_manuscript_data_figs.R")
```

**Main text**

``` r
source("code/09.0_manuscript.rmd")
```

**Supporting Information**

``` r
source("code/10.0_Supporting_Information.Rmd")
```
