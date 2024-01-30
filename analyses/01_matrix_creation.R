################################################################################
# script : 01_matrix_creation.R
#
# Commentary : This script allows the creation of the 3 matrices (occurrence, 
#              abundance, and biomass), and also saves them in the "data" 
#              folder, allowing to start the analyses directly from script 02.
#
################################################################################

library(tidyverse)

# Source R code

source(here::here("R", "occu_matrix.R"))
source(here::here("R", "abun_matrix.R"))
source(here::here("R", "biom_matrix.R"))

# Occurrence 

melted = make_ssmatrix_melted(30) # Molten matrix with species w/ n minimum occu
melted_filt = rid_rare_outside() # Filters out species with < 5 obs outside
melted_cov = add_covariates() # Adds covariates to each observations

save(melted_cov, file = here::here("data", "data", "melted_cov.Rdata"))

# Abundance

melted_ab = make_ssmatrix_melted_ab(30) # Molten matrix with species w/ n mini
melted_filt_ab = rid_rare_outside_ab() # Filters out species with < 5 obs out
melted_cov_ab = add_covariates_ab() # Adds covariates to each observations

save(melted_cov_ab, file = here::here("data", "data", "melted_cov_ab.Rdata"))

# Biomass

melted_bm = make_ssmatrix_melted_bm(30) # Molten matrix with species w/ n mini
melted_filt_bm = rid_rare_outside_bm() # Filters out species with < 5 obs out
melted_cov_bm = add_covariates_bm() # Adds covariates to each observations

save(melted_cov_bm, file = here::here("data", "data", "melted_cov_bm.Rdata"))