################################################################################
# script : 00_load_and_tidy.R
#
# Commentary : This script loads the necessary data to start the analyses
#
################################################################################

library(tidyverse)

source(here::here("R", "load_data.R"))
source(here::here("R", "tidy_sites.R"))

load_rdata()
covariates_corrected = tidy_sites()
all_cov = all_cov_matrix()
