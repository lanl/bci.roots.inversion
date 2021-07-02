Scripts in this repo are a part of the workflow associated with the
following manuscript:

Chitra-Tarak, R, C Xu, S Aguilar, K Anderson-Teixeira, J Chambers, M
Detto, B Faybishenko, RA Fisher, R Knox, C Koven, L Kueppers, N Kunert,
SJ Kupers, NG McDowell, BD Newman, SR Paton, R Pérez, L Ruiz, L Sack, JM
Warren, BT Wolfe, C Wright, SJ Wright, J Zailaa, SM McMahon (2021)
Hydraulically vulnerable trees survive on deep-water access during
droughts in a tropical forest. New Phytologist. https:
//doi.org/10.1111/nph.17464

For the full work-flow and output datasets associated with the above
manuscript see the following archived dataset:

Chitra-Tarak R, Xu C, Aguilar S, Anderson-Teixeira K, Chambers J, Detto
M, Faybishenko B, Fisher R, Knox R, Koven C et al. 2020. Soil water
potentials (1990–2018) from a calibrated ELM-FATES, and rooting depth
analyses scripts, PA-BCI, Panama. 2.0. NGEE Tropics Data Collection.
doi: 10.15486/ngt/1696806

# Workflow Part II: Inverse root modeling and post-processing

# Pre-requisite data

**Growth data**

Calculate species-specific growth time series from forest inventory
data.

    source("code/01.0_Preparing_stem_growth_time_series.R")

**Mortality data**

    source("code/02.0_demographic_data_BCI.R")

#### All pre-requisites

Gather all pre-requisites for rooting depth inverse model (effective
rooting depth ) in a place and traits data for post-processing

-   Soil Water Potentials (outputs from ELM-FATES)
-   VPD-GPP relationship
-   Leaf vulnerability curves
-   LAI seasonality
-   Hydraulic traits

<!-- -->

    source("code/05.0_Prepare_data_for_correlative analyses.R")

Script to load data generated from the data preparation script above for
use in inverse model as well as manuscript. It’s used within the rest of
the scripts.

<!-- ```{r Chunk 06.0, eval = FALSE} -->
<!-- source("code/06.0_load.R") -->
<!-- ``` -->

# Model Setup

Set-up and run effective rooting depth models, evaluate model, and join
ERD from the predicted model with hydraulic traits and mortality data,
and calculate realized hydraulic risk

    source("code/07.0_PhenoDemoTraitsPsi.R")

# Manuscript

**Main text**

    # Two files are used within the report
    # For data, 06.0_load.R and 
    # For figures and data, code/08.0_manuscript_data_figs.R
    if (!require("groundhog")) install.packages("groundhog")
    library(groundhog)
    groundhog.folder <- paste0("groundhog.library")
    if(!dir.exists(file.path(groundhog.folder))) {dir.create(file.path(groundhog.folder))}
    set.groundhog.folder(groundhog.folder)

    groundhog.day = "2021-01-01"
    groundhog.library('rmarkdown', groundhog.day)
    # .rs.restartR() to restart R session
    groundhog.library('bookdown', groundhog.day)
    # groundhog.library('knitr', groundhog.day)
    # rmarkdown::render("09.0_manuscript.rmd", output_format = "bookdown::html_document2")
    rmarkdown::render("09.0_manuscript.rmd", output_format = "bookdown::word_document2")

**Supporting Information**

    # rmarkdown::render("10.0_Supporting_Information.Rmd", output_format = "bookdown::html_document2")
    rmarkdown::render("10.0_Supporting_Information.Rmd", output_format = "bookdown::word_document2")

    remotes::install_github("ropenscilabs/reviewer")
    modified_file <- system.file("extdata/CheatSheet-modified.Rmd", package = "reviewer")
    reference_file <- system.file("extdata/CheatSheet.Rmd", package = "reviewer")
    library(reviewer)
    result <- diff_rmd(modified_file, reference_file)
    groundhog.library('diffobj', groundhog.day)
    browseVignettes("diffobj")
    remotes::install_github("ropenscilabs/trackmd")
    trackmd::trackChanges("09.0_manuscript.rmd")

## Copyrights

© 2021. Triad National Security, LLC. All rights reserved. LANL C20130.

This program was produced under U.S. Government contract
89233218CNA000001 for Los Alamos National Laboratory (LANL), which is
operated by Triad National Security, LLC for the U.S. Department of
Energy/National Nuclear Security Administration. All rights in the
program are reserved by Triad National Security, LLC, and the U.S.
Department of Energy/National Nuclear Security Administration. The
Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to
reproduce, prepare derivative works, distribute copies to the public,
perform publicly and display publicly, and to permit others to do so.

## Open Source Redistribution License

This program is open source under the BSD-3 License.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

    3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS
IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
