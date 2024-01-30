################################################################################
# script : 02_effect_sizes.R
#
# Commentary : This script computes the effect sizes of protection on all 3 
#              facets (occu/abun/biom). Analyses can be started directly from 
#              here if needed thanks to the "data" folder, containing data that
#              have been cleaned and set up.
#
#              !! This part should be used with mclapply, make sure your system 
#                 allows parallel computations or it would take too long !!
#
################################################################################

load(here::here("data", "data", "melted_cov_bm.RData"))
load(here::here("data", "data", "melted_cov_ab.RData"))
load(here::here("data", "data", "melted_cov.RData"))
load(here::here("data", "data", "all_cov.RData"))

source(here::here("R", "effect_sizes.R"))

OS_occu = compute_ES_occu(x) # Set x = number of cores you're able to use
df_ab = compute_ES_abun(x)
df_bm_product = compute_ES_bm(x) # This one needs to have OS_occu computed already !

# Plots

histograms()
