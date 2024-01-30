################################################################################
# script : 03_traits_models.R
#
# Commentary : This script computes models for each level of protection, and 
#              and each facet of biodiversity, then plots the results.
#
################################################################################
library(tidyverse)

load(here::here("data", "raw_data", "traits.RData"))
load(here::here("data", "data", "melted_cov.RData"))
load(here::here("outputs", "OS_occu.RData"))
load(here::here("outputs", "df_ab.RData"))
load(here::here("outputs", "df_bm.RData"))

source(here::here("R", "prep_models.R"))
source(here::here("R", "models.R"))

melted_rares = add_rarity()
data_occu = prep_data_occu()
data_abun = prep_data_abun()
data_biom = prep_data_bm()

occu_full = occu_full()
occu_high = occu_high()
occu_light = occu_light()

abun_full = abun_full()
abun_high = abun_high()
abun_light = abun_light()

biom_full = biom_full()
biom_high = biom_high()
biom_light = biom_light()

#  FX plots 

all_rarity_trophic()
threeway_interaction()
