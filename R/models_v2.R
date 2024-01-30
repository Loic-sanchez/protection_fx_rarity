# Separate all DF

library(tidyverse)

prep_data_occu = function(){

  # Get rarity 
  
  traits$Species = traits$CURRENT_SPECIES_NAME
  
  melted_rares = melted_cov |> 
    filter(Presence > 0) |>  
    group_by(Species) |>  
    summarise(n_occu = n()) |>  
    mutate(decile_rank = ntile(n_occu, 10))
  
  df_occu_traits = ES_occu |>  
    left_join(traits, by = "Species") |>    
    left_join(melted_rares |>  dplyr::select(Species, decile_rank), by = "Species") |>  
    mutate(log_length = log10(MaxLength)) 

  df_occu_full = df_occu_traits |>  
    drop_na(RR_full, Trophic.Level, log_length) |> 
    mutate(winner = ifelse(RR_full > 1, 1, 0))
  df_occu_high = df_occu_traits |> 
    drop_na(RR_high, Trophic.Level, log_length)|> 
    mutate(winner = ifelse(RR_high > 1, 1, 0))
  df_occu_light = df_occu_traits |> 
    drop_na(RR_light, Trophic.Level, log_length)|> 
    mutate(winner = ifelse(RR_light > 1, 1, 0))
  
  list_occu = list(df_occu_full, df_occu_high, df_occu_light)
  return(list_occu)
}

prep_data_ab = function(){
  
  # Get rarity 
  
  traits$Species = traits$CURRENT_SPECIES_NAME
  
  melted_rares = melted_cov |> 
    filter(Presence > 0) |>  
    group_by(Species) |>  
    summarise(n_occu = n()) |>  
    mutate(decile_rank = ntile(n_occu, 10))
  
  df_ab_traits = ES_ab |> 
    left_join(traits, by = "Species") |> 
    left_join(melted_rares %>% dplyr::select(Species, decile_rank), by = "Species") |> 
    mutate(log_length = log10(MaxLength))
  
  df_ab_full = df_ab_traits |>  
    drop_na(IRR_full, Trophic.Level, log_length) |> 
    mutate(winner = ifelse(IRR_full > 1, 1, 0))
  df_ab_high = df_ab_traits |>  
    drop_na(IRR_high, Trophic.Level, log_length) |> 
    mutate(winner = ifelse(IRR_high > 1, 1, 0))
  df_ab_light = df_ab_traits |>  
    drop_na(IRR_light, Trophic.Level, log_length) |> 
    mutate(winner = ifelse(IRR_light > 1, 1, 0))

  list_abun = list(df_ab_full, df_ab_high, df_ab_light)
  return(list_abun)

}

prep_data_bm = function(){
  
  # Get rarity 
  
  traits$Species = traits$CURRENT_SPECIES_NAME
  
  melted_rares = melted_cov |> 
    filter(Presence > 0) |>  
    group_by(Species) |>  
    summarise(n_occu = n()) |>  
    mutate(decile_rank = ntile(n_occu, 10))
  
  df_bm_product = ES_bm |> 
    left_join(ES_occu |> 
                dplyr::select(Species, 
                              estimate_outside, 
                              estimate_full, 
                              estimate_high,
                              estimate_light, 
                              AUC), 
              by = "Species") |> 
    rowwise() |> 
    mutate(product_full = estimate_full.x*estimate_full.y,
           product_high = estimate_high.x*estimate_high.y,
           product_light = estimate_light.x*estimate_light.y,
           product_outside = estimate_outside.x*estimate_outside.y) |> 
    mutate(IRR_product_full = product_full/product_outside,
           IRR_product_high = product_high/product_outside,
           IRR_product_light = product_light/product_outside) |> 
    drop_na(estimate_outside.y)
  
  df_bm_traits = df_bm_product |> 
    left_join(traits, by = "Species") |>  
    left_join(melted_rares %>% dplyr::select(Species, decile_rank), by = "Species") |> 
    mutate(log_length = log10(MaxLength))
  
  df_bm_full = df_bm_traits |> 
    drop_na(IRR_product_full, Trophic.Level, log_length) |> 
    mutate(winner = ifelse(IRR_product_full > 1, 1, 0))
  df_bm_high = df_bm_traits |>  
    drop_na(IRR_product_high, Trophic.Level, log_length) |> 
    mutate(winner = ifelse(IRR_product_high > 1, 1, 0))
  df_bm_light = df_bm_traits |>  
    drop_na(IRR_product_light, Trophic.Level, log_length) |> 
    mutate(winner = ifelse(IRR_product_high > 1, 1, 0))
  
  list_biom = list(df_bm_full, df_bm_high, df_bm_light)

}

gamma_occurrence = function(){

  # Compute 9 models 
  
  data_occu = prep_data_occu()
  data_of = data_occu[[1]]
  data_oh = data_occu[[2]]
  data_ol = data_occu[[3]]
  
  gamma_occu_full = glm(RR_full ~ decile_rank * MaxLength * poly(Trophic.Level, 2),
                        data = data_of,
                        family = Gamma(link = "log"))
  gamma_of_aic = MASS::stepAIC(gamma_occu_full, direction = "both")
  RSQ_of = pscl::pR2(gamma_of_aic)["r2CU"]
  
  gamma_occu_high = glm(RR_high ~ decile_rank * MaxLength * poly(Trophic.Level, 2),
                        data = data_oh,
                        family = Gamma(link = "log"))
  gamma_oh_aic = MASS::stepAIC(gamma_occu_high, direction = "both")
  RSQ_oh = pscl::pR2(gamma_oh_aic)["r2CU"]

  gamma_occu_light = glm(RR_light ~ decile_rank * MaxLength * poly(Trophic.Level, 2),
                        data = data_ol,
                        family = Gamma(link = "log"))
  gamma_ol_aic = MASS::stepAIC(gamma_occu_light, direction = "both")
  RSQ_ol = pscl::pR2(gamma_ol_aic)["r2CU"]
  
  list_gamma_occu = list(gamma_of_aic, gamma_oh_aic, gamma_ol_aic)
  
}

gamma_abundance = function(){
  
  data_ab = prep_data_ab()
  data_af = data_ab[[1]]
  data_ah = data_ab[[2]]
  data_al = data_ab[[3]]
  
  gamma_ab_full = glm(IRR_full ~ decile_rank * MaxLength * I(Trophic.Level^2),
                      data = data_af,
                      family = Gamma(link = "log"))
  gamma_af_aic = MASS::stepAIC(gamma_ab_full, direction = "both")
  RSQ_af = pscl::pR2(gamma_af_aic)["r2CU"]
  
  gamma_ab_high = glm(IRR_high ~ decile_rank * MaxLength * poly(Trophic.Level, 2),
                      data = data_ah,
                      family = Gamma(link = "log"))
  gamma_ah_aic = MASS::stepAIC(gamma_ab_high, direction = "both")
  RSQ_ah = pscl::pR2(gamma_ah_aic)["r2CU"]

  gamma_ab_light = glm(IRR_light ~ decile_rank * MaxLength * poly(Trophic.Level, 2),
                      data = data_al,
                      family = Gamma(link = "log"))
  gamma_al_aic = MASS::stepAIC(gamma_ab_light, direction = "both")
  RSQ_al = pscl::pR2(gamma_al_aic)["r2CU"]
  
  list_gamma_abun = list(gamma_af_aic, gamma_ah_aic, gamma_al_aic)
  
}

gamma_biomass = function(){
  
  data_bm = prep_data_bm()
  data_bf = data_bm[[1]]
  data_bh = data_bm[[2]]
  data_bl = data_bm[[3]]
  
  gamma_bm_full = glm(IRR_product_full ~ decile_rank * MaxLength * poly(Trophic.Level, 2),
                      data = data_bf,
                      family = Gamma(link = "log"))
  gamma_bf_aic = MASS::stepAIC(gamma_bm_full, direction = "both")
  RSQ_bf = pscl::pR2(gamma_bf_aic)["r2CU"]
  
  gamma_bm_high = glm(IRR_product_high ~ decile_rank * MaxLength * poly(Trophic.Level, 2),
                      data = data_bh,
                      family = Gamma(link = "log"))
  gamma_bh_aic = MASS::stepAIC(gamma_bm_high, direction = "both")
  RSQ_bh = pscl::pR2(gamma_bh_aic)["r2CU"]
  
  gamma_bm_light = glm(IRR_product_light ~ decile_rank * MaxLength * poly(Trophic.Level, 2),
                      data = data_bl,
                      family = Gamma(link = "log"))
  gamma_bl_aic = MASS::stepAIC(gamma_bm_light, direction = "both")
  RSQ_bl = pscl::pR2(gamma_bl_aic)["r2CU"]
  
  list_gamma_biom = list(gamma_bf_aic, gamma_bh_aic, gamma_bl_aic)

}

# Test if pR2 are random for full protection
#   Returns an alpha error vector, 0 indicating that the model is not random
#   mean = 0.0133; sd = 0.011 --> The pR2 are not random

randomize_traits = function(){
  
  df_occu_full = prep_data_occu()[[1]]
  df_occu_high = prep_data_occu()[[2]]
  df_occu_light = prep_data_occu()[[3]]
  
  df_ab_full = prep_data_ab()[[1]]
  df_ab_high = prep_data_ab()[[2]]
  df_ab_light = prep_data_ab()[[3]]
  
  df_bm_full = prep_data_bm()[[1]]
  df_bm_high = prep_data_bm()[[2]]
  df_bm_light = prep_data_bm()[[3]]
  
  original_occu_full = prep_data_occu()[[1]]
  original_occu_high = prep_data_occu()[[2]]
  original_occu_light = prep_data_occu()[[3]]
  
  original_ab_full = prep_data_ab()[[1]]
  original_ab_high = prep_data_ab()[[2]]
  original_ab_light = prep_data_ab()[[3]]
  
  original_bm_full = prep_data_bm()[[1]]
  original_bm_high = prep_data_bm()[[2]]
  original_bm_light = prep_data_bm()[[3]]
  
  randomized_rsq_of = vector()
  randomized_rsq_oh = vector()
  randomized_rsq_ol = vector()
  
  randomized_rsq_af = vector()
  randomized_rsq_ah = vector()
  randomized_rsq_al = vector()
  
  randomized_rsq_bf = vector()
  randomized_rsq_bh = vector()
  randomized_rsq_bl = vector()
  
  n_it = 999
  
  set.seed(001)
  
  # occu full 
  
  for (i in 1:n_it) {
    
    species_order = sample(unique(df_occu_full$Species))
    df_occu_full$RR_full <-
      original_occu_full$RR_full[match(df_occu_full$Species, species_order)]
    
    tryCatch({
    mod_occu_full = glm(
      RR_full ~ decile_rank * MaxLength * poly(Trophic.Level, 2),
      data = df_occu_full,
      family = Gamma(link = "log")
    )
    
    
    mod_occu_full_AIC = MASS::stepAIC(mod_occu_full, 
                                      direction = "both", 
                                      trace = 0)}, error = function(e) mod_occu_full_AIC <<- NULL)
    randomized_rsq_of[[i]] = ifelse(is.null(mod_occu_full_AIC) == F, pscl::pR2(mod_occu_full_AIC)["r2CU"], NA)
    
  }
  
  # occu high
  
  for (i in 1:n_it) {
    
    species_order = sample(unique(df_occu_high$Species))
    df_occu_high$RR_high <-
      original_occu_high$RR_high[match(df_occu_high$Species, species_order)]
    
    tryCatch({
      mod_occu_high = glm(
        RR_high ~ decile_rank * MaxLength * poly(Trophic.Level, 2),
        data = df_occu_high,
        family = Gamma(link = "log")
      )
      
      
      mod_occu_high_AIC = MASS::stepAIC(mod_occu_high, 
                                        direction = "both", 
                                        trace = 0)}, error = function(e) mod_occu_high_AIC <<- NULL)
    randomized_rsq_oh[[i]] = ifelse(is.null(mod_occu_high_AIC) == F, pscl::pR2(mod_occu_high_AIC)["r2CU"], NA)
    
  } 
  
  # occu light
  
  for (i in 1:n_it) {
    
    species_order = sample(unique(df_occu_light$Species))
    df_occu_light$RR_light <-
      original_occu_light$RR_light[match(df_occu_light$Species, species_order)]
    
    tryCatch({
      mod_occu_light = glm(
        RR_light ~ decile_rank * MaxLength * poly(Trophic.Level, 2),
        data = df_occu_light,
        family = Gamma(link = "log")
      )
      
      
      mod_occu_light_AIC = MASS::stepAIC(mod_occu_light, 
                                        direction = "both", 
                                        trace = 0)}, error = function(e) mod_occu_light_AIC <<- NULL)
    randomized_rsq_ol[[i]] = ifelse(is.null(mod_occu_light_AIC) == F, pscl::pR2(mod_occu_light_AIC)["r2CU"], NA)
    
  } 
  
  hist(randomized_rsq_of)
  alpha_of = 1 - ((sum(randomized_rsq_of < RSQ_of, na.rm = T) / length(na.omit(randomized_rsq_of))))
  
  hist(randomized_rsq_oh)
  alpha_oh = 1 - ((sum(randomized_rsq_oh < RSQ_oh, na.rm = T) / length(na.omit(randomized_rsq_oh))))
  
  hist(randomized_rsq_ol)
  alpha_ol = 1 - ((sum(randomized_rsq_ol < RSQ_ol, na.rm = T) / length(na.omit(randomized_rsq_ol))))
  
  # abun full
  
  for (i in 1:n_it) {
    species_order = sample(unique(df_ab_full$Species))
    df_ab_full$IRR_full <-
      original_ab_full$IRR_full[match(df_ab_full$Species, species_order)]
    
    tryCatch({
    mod_ab_full = glm(
      IRR_full ~ decile_rank * MaxLength * I(Trophic.Level ^ 2),
      data = df_ab_full,
      family = Gamma(link = "log")
    )
    
      mod_ab_full_AIC = MASS::stepAIC(mod_ab_full, 
                                      direction = "both", 
                                      trace = 0)}, error = function(e) mod_ab_full_AIC <<- NULL)
    
      randomized_rsq_af[[i]] = ifelse(is.null(mod_ab_full_AIC) == F, pscl::pR2(mod_ab_full_AIC)["r2CU"], NA)
    
  }
  
  # abun high

  for (i in 1:n_it) {
    species_order = sample(unique(df_ab_high$Species))
    df_ab_high$IRR_high <-
      original_ab_high$IRR_high[match(df_ab_high$Species, species_order)]
    
    tryCatch({
      mod_ab_high = glm(
        IRR_high ~ decile_rank * MaxLength * poly(Trophic.Level, 2),
        data = df_ab_high,
        family = Gamma(link = "log")
      )
      
      mod_ab_high_AIC = MASS::stepAIC(mod_ab_high, 
                                      direction = "both", 
                                      trace = 0)}, error = function(e) mod_ab_high_AIC <<- NULL)
    
    randomized_rsq_ah[[i]] = ifelse(is.null(mod_ab_high_AIC) == F, pscl::pR2(mod_ab_high_AIC)["r2CU"], NA)
    
  }
  
  # abun light 
  
  for (i in 1:n_it) {
    species_order = sample(unique(df_ab_light$Species))
    df_ab_light$IRR_light <-
      original_ab_light$IRR_light[match(df_ab_light$Species, species_order)]
    
    tryCatch({
      mod_ab_light = glm(
        IRR_light ~ decile_rank * MaxLength * poly(Trophic.Level, 2),
        data = df_ab_light,
        family = Gamma(link = "log")
      )
      
      mod_ab_light_AIC = MASS::stepAIC(mod_ab_light, 
                                      direction = "both", 
                                      trace = 0)}, error = function(e) mod_ab_light_AIC <<- NULL)
    
    randomized_rsq_al[[i]] = ifelse(is.null(mod_ab_light_AIC) == F, pscl::pR2(mod_ab_light_AIC)["r2CU"], NA)
    
  }
  
  hist(randomized_rsq_af)
  alpha_af = 1 - ((sum(randomized_rsq_af < RSQ_af, na.rm = T) / length(na.omit(randomized_rsq_af))))
  
  hist(randomized_rsq_ah)
  alpha_ah = 1 - ((sum(randomized_rsq_ah < RSQ_ah, na.rm = T) / length(na.omit(randomized_rsq_ah))))
  
  hist(randomized_rsq_al)
  alpha_al = 1 - ((sum(randomized_rsq_al < RSQ_al, na.rm = T) / length(na.omit(randomized_rsq_al))))
  
  # biom full
  
  for (i in 1:n_it) {
    
    species_order = sample(unique(df_bm_full$Species))
    df_bm_full$IRR_product_full <-
      original_bm_full$IRR_product_full[match(df_bm_full$Species, species_order)]
    
    tryCatch({
    mod_bm_full = glm(
      IRR_product_full ~ decile_rank * MaxLength * poly(Trophic.Level, 2),
      data = df_bm_full,
      family = Gamma(link = "log")
    )
    
    
    mod_bm_full_AIC = MASS::stepAIC(mod_bm_full, 
                                    direction = "both", 
                                    trace = 0)}, error = function(e) mod_bm_full_AIC <<- NULL)
    
    randomized_rsq_bf[[i]] = ifelse(is.null(mod_bm_full_AIC) == F, pscl::pR2(mod_bm_full_AIC)["r2CU"], NA)
  }
  
  # biom high
  
  for (i in 1:n_it) {
    
    species_order = sample(unique(df_bm_high$Species))
    df_bm_high$IRR_product_high <-
      original_bm_high$IRR_product_high[match(df_bm_high$Species, species_order)]
    
    tryCatch({
      mod_bm_high = glm(
        IRR_product_high ~ decile_rank * MaxLength * poly(Trophic.Level, 2),
        data = df_bm_high,
        family = Gamma(link = "log")
      )
      
      
      mod_bm_high_AIC = MASS::stepAIC(mod_bm_high, 
                                      direction = "both", 
                                      trace = 0)}, error = function(e) mod_bm_high_AIC <<- NULL)
    
    randomized_rsq_bh[[i]] = ifelse(is.null(mod_bm_high_AIC) == F, pscl::pR2(mod_bm_high_AIC)["r2CU"], NA)
  }
  
  # biom light
  
  for (i in 1:n_it) {
    
    species_order = sample(unique(df_bm_light$Species))
    df_bm_light$IRR_product_light <-
      original_bm_light$IRR_product_light[match(df_bm_light$Species, species_order)]
    
    tryCatch({
      mod_bm_light = glm(
        IRR_product_light ~ decile_rank * MaxLength * poly(Trophic.Level, 2),
        data = df_bm_light,
        family = Gamma(link = "log")
      )
      
      
      mod_bm_light_AIC = MASS::stepAIC(mod_bm_light, 
                                      direction = "both", 
                                      trace = 0)}, error = function(e) mod_bm_light_AIC <<- NULL)
    
    randomized_rsq_bl[[i]] = ifelse(is.null(mod_bm_light_AIC) == F, pscl::pR2(mod_bm_light_AIC)["r2CU"], NA)
  }
  
  hist(randomized_rsq_bf)
  alpha_bf = 1 - ((sum(randomized_rsq_bf < RSQ_bf, na.rm = T) / length(na.omit(randomized_rsq_bf))))
  
  hist(randomized_rsq_bh)
  alpha_bh = 1 - ((sum(randomized_rsq_bh < RSQ_bh, na.rm = T) / length(na.omit(randomized_rsq_bh))))
  
  hist(randomized_rsq_ol)
  alpha_bl = 1 - ((sum(randomized_rsq_bl < RSQ_bl, na.rm = T) / length(na.omit(randomized_rsq_bl))))

  alpha_vec = c(alpha_of, alpha_af, alpha_bf)
  return(alpha_vec)
  
} 

# Check model residuals for full protection and save them

check_residuals_gamma = function(){

  jpeg(file = here::here("plots_revised", "suppl_residuals_traits.jpeg"),
       units = "px", 
       width = 2540, 
       height = 2000, 
       res = 300)
  
  par(mfrow = c(3, 2))
  
  res_of1 = DHARMa::testUniformity(gamma_of_aic)
  res_of2 = DHARMa::testQuantiles(gamma_of_aic, quantiles = .5)

  res_af1 = DHARMa::testUniformity(gamma_af_aic)
  res_af2 = DHARMa::testQuantiles(gamma_af_aic, quantiles = .5)

  res_bf1 = DHARMa::testUniformity(gamma_bf_aic)
  res_bf2 = DHARMa::testQuantiles(gamma_bf_aic, quantiles = .5)

  dev.off()

}

# Save models in plots_revised folder
#
# - Main figure of both abundance and biomass plots: interactions at same scale without conf. intervals
# - 2 suppl. figures (abundance & biomass) with free y scales and conf. intervals

plot_af_and_bf = function(){
  
  # extract model output for Trophic X Rarity decile X small species
  
  int_af_small = visreg::visreg(gamma_abundance()[[1]],
                                type = "conditional",
                                xvar = "Trophic.Level",
                                scale = "response", 
                                by = "decile_rank",
                                nn = 505,
                                cond = list(MaxLength = 9),
                                gg = F)$fit
  
  # extract model output for Trophic X Rarity decile X medium species
  
  int_af_medium = visreg::visreg(gamma_abundance()[[1]],
                                 type = "conditional",
                                 xvar = "Trophic.Level",
                                 scale = "response", 
                                 by = "decile_rank",
                                 nn = 505,
                                 cond = list(MaxLength = 22),
                                 gg = F)$fit
  
  # extract model output for Trophic X Rarity decile X large species
  
  int_af_large = visreg::visreg(gamma_abundance()[[1]],
                                type = "conditional",
                                xvar = "Trophic.Level",
                                scale = "response", 
                                by = "decile_rank",
                                nn = 505,
                                cond = list(MaxLength = 60),
                                gg = F)$fit
  
  all_sizes_ab = rbind(int_af_small, int_af_medium, int_af_large) # Join 3 outputs
  
  all_sizes_ab = all_sizes_ab |> 
    dplyr::mutate(Rarity = ifelse(decile_rank == 5, "Rare", 
                                  ifelse(decile_rank == 8, "Common", "Very frequent")),
                  length = ifelse(MaxLength == 9, "Small", 
                                  ifelse(MaxLength == 22, "Medium", "Large")))
  
  all_sizes_ab$Rarity <- factor(all_sizes_ab$Rarity, levels = c("Rare", "Common", "Very frequent"))
  all_sizes_ab$length <- factor(all_sizes_ab$length, levels = c("Small", "Medium", "Large"))
  
  # extract model output for Trophic X Rarity decile X small species
  
  int_bf_small = visreg::visreg(gamma_biomass()[[1]],
                         type = "conditional",
                         xvar = "Trophic.Level",
                         scale = "response", 
                         by = "decile_rank",
                         cond = list(MaxLength = 9),
                         nn = 505,
                         gg = F)$fit
  int_bf_small$length = "Small"
  
  # extract model output for Trophic X Rarity decile X medium species
  
  int_bf_medium = visreg::visreg(gamma_biomass()[[1]],
                          type = "conditional",
                          xvar = "Trophic.Level",
                          scale = "response", 
                          by = "decile_rank",
                          nn = 505,
                          cond = list(MaxLength = 20),
                          gg = F)$fit
  int_bf_medium$length = "Medium"
  
  # extract model output for Trophic X Rarity decile X large species
  
  int_bf_large = visreg::visreg(gamma_biomass()[[1]],
                         type = "conditional",
                         xvar = "Trophic.Level",
                         scale = "response", 
                         by = "decile_rank",
                         cond = list(MaxLength = 60),               
                         nn = 505,
                         gg = F)$fit
  int_bf_large$length = "Large"
  
  all_sizes = rbind(int_bf_small, int_bf_medium, int_bf_large) # Join all 3
  
  all_sizes$length <- factor(all_sizes$length, levels = c("Small", "Medium", "Large"))
  all_sizes = all_sizes |> 
    dplyr::mutate(Rarity = ifelse(decile_rank == 5, "Rare", 
                                  ifelse(decile_rank == 8, "Common", "Very frequent")))
  all_sizes$Rarity <- factor(all_sizes$Rarity, levels = c("Rare", "Common", "Very frequent"))
  
  all_sizes$Type = "Biomass"
  all_sizes_ab$Type = "Abundance"
  
  ab_and_bm = rbind(all_sizes_ab |> dplyr::select(-IRR_full), 
                    all_sizes |> dplyr::select(-IRR_product_full))
  
  # Make ggplot without confidence intervals
  
  # pal = harrypotter::hp(3, option = "Gryffindor", begin = 0, end = 1, direction = -1)
  # pal = harrypotter::hp(3, option = "Always", begin = 0.1, end = 0.9, direction = -1)
  # scales::show_col(pal)
  scales::show_col(pal_custom)
  
  # pal_custom = c("#AD0000FF", "#428185FF", "#ABE2B2FF")
  
  pal_custom = c("#AD0000FF", "#F06400FF", "#FFa700FF")
  
  ggplot(ab_and_bm, aes(x = Trophic.Level, y = visregFit)) +
    geom_line(linewidth = 1, aes(color = Rarity)) +
    facet_grid(cols = vars(length), rows = vars(Type), scales = "free_y", ) + 
    ylab("Effect size of high protection") +
    xlab("Trophic level") +
    theme_bw() +
    coord_cartesian(y = c(0, 10)) +
    geom_hline(yintercept = 1, linetype = 2, linewidth = 0.5) +
    scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
    scale_color_manual(values = pal_custom) +
    # harrypotter::scale_color_hp_d("Gryffindor", direction  = -1) +
    theme(strip.text = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 10),
          strip.background = element_blank(),
          panel.spacing = unit(2, "lines"))
  
  ggsave(here::here("plots_revised", "int_af_bf.jpeg"), 
         units = "px", 
         width = 2540, 
         height = 2000, 
         dpi = 300)
  
  dev.off()
  
  # GGplot with ribbons (abundance)
  
  ggplot(all_sizes_ab, aes(x = Trophic.Level, 
                           y = visregFit)) +
    geom_ribbon(aes(ymin = visregLwr, 
                    ymax = visregUpr,
                    fill = Rarity),
                alpha = 0.5) +
    geom_line(linewidth = 1, aes(color = Rarity)) +
    facet_grid(cols = vars(length),
               rows = vars(Rarity), 
               scales = "free_y", ) + 
    ylab("Effect size of high protection") +
    xlab("Trophic level") +
    theme_bw() +
    # coord_cartesian(y = c(0, 10)) +
    geom_hline(yintercept = 1, 
               linetype = 2, 
               linewidth = 0.5) +
    # scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
    # harrypotter::scale_fill_hp_d("Gryffindor", direction  = -1) +
    # harrypotter::scale_color_hp_d("Gryffindor", direction  = -1) +
    scale_color_manual(values = pal_custom) +
    scale_fill_manual(values = pal_custom) +
    theme(strip.text = element_text(size = 10, 
                                    face = "bold"),
          axis.text = element_text(size = 10),
          strip.background = element_blank(),
          panel.spacing = unit(2, "lines"))
  
  ggsave(here::here("plots_revised", "suppl_int_af.jpeg"), 
         units = "px", 
         width = 2540, 
         height = 2000, 
         dpi = 300)
  dev.off()
  
  # GGplot with ribbons (biomass)
  
  ggplot(all_sizes, aes(x = Trophic.Level, 
                           y = visregFit)) +
    geom_line(linewidth = 1, aes(color = Rarity)) +
    geom_ribbon(aes(ymin = visregLwr, 
                    ymax = visregUpr,
                    fill = Rarity),
                alpha = 0.5) +
    facet_grid(cols = vars(length),
               rows = vars(Rarity), 
               scales = "free_y", ) + 
    ylab("Effect size of high protection") +
    xlab("Trophic level") +
    theme_bw() +
    # coord_cartesian(y = c(0, 10)) +
    geom_hline(yintercept = 1, 
               linetype = 2, 
               linewidth = 0.5) +
    # scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
    # harrypotter::scale_fill_hp_d("Gryffindor", direction = -1) +
    # harrypotter::scale_color_hp_d("Gryffindor", direction = -1) +
    scale_color_manual(values = pal_custom) +
    scale_fill_manual(values = pal_custom) +
    theme(strip.text = element_text(size = 10, 
                                    face = "bold"),
          axis.text = element_text(size = 10),
          strip.background = element_blank(),
          panel.spacing = unit(2, "lines"))
  
  ggsave(here::here("plots_revised", "suppl_int_bf.jpeg"), 
         units = "px", 
         width = 2540, 
         height = 2000, 
         dpi = 300)
  dev.off()
  
}

plot_of = function(){
  
  int_of = visreg::visreg(gamma_occurrence()[[1]],
                          type = "conditional",
                          xvar = "Trophic.Level",
                          scale = "response", 
                          by = "MaxLength",
                          gg = F)$fit
  
  int_of = int_of |> 
    dplyr::mutate(Length = ifelse(MaxLength == 9, "Small", 
                                  ifelse(MaxLength == 21, "Medium", "Large")))
  int_of$Length <- factor(int_of$Length, levels = c("Small", "Medium", "Large"))
  
  # With conf. int
  
  ggplot(int_of, aes(x = Trophic.Level, y = visregFit)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr),
                alpha = 0.5,
                fill = "#03045e") +
    facet_wrap(vars(Length)) + 
    ylab("Effect size on occurrence of high protection") +
    xlab("Trophic level") +
    theme_bw() +
    geom_hline(yintercept = 1, linetype = 2, linewidth = 0.5) +
    theme(strip.text = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 10),
          strip.background = element_blank(),
          panel.spacing = unit(2, "lines"))
  
  ggsave(here::here("plots_revised", "int_of.jpeg"), 
         units = "px", 
         width = 2540, 
         height = 2000, 
         dpi = 300)
  dev.off()

}

plot_bh = function(){
  
  int_bh_small = visreg::visreg(gamma_biomass()[[2]],
                                type = "conditional",
                                xvar = "Trophic.Level",
                                scale = "response", 
                                by = "decile_rank",
                                nn = 505,
                                cond = list(MaxLength = 9),
                                gg = F)$fit
  
  int_bh_medium = visreg::visreg(gamma_biomass()[[2]],
                                type = "conditional",
                                xvar = "Trophic.Level",
                                scale = "response", 
                                by = "decile_rank",
                                nn = 505,
                                cond = list(MaxLength = 24),
                                gg = F)$fit
  
  int_bh_large =  visreg::visreg(gamma_biomass()[[2]],
                                 type = "conditional",
                                 xvar = "Trophic.Level",
                                 scale = "response", 
                                 by = "decile_rank",
                                 nn = 505,
                                 cond = list(MaxLength = 61),
                                 gg = F)$fit
  
  all_sizes_bh = rbind(int_bh_large, 
                       int_bh_medium, 
                       int_bh_small)
  
  all_sizes_bh = all_sizes_bh |> 
    dplyr::mutate(Rarity = ifelse(decile_rank == 4, "Rare", 
                                  ifelse(decile_rank == 7, "Common", "Very frequent")),
                  length = ifelse(MaxLength == 9, "Small", 
                                  ifelse(MaxLength == 24, "Medium", "Large")))
  
  all_sizes_bh$Rarity <- factor(all_sizes_bh$Rarity, levels = c("Rare", "Common", "Very frequent"))
  all_sizes_bh$length <- factor(all_sizes_bh$length, levels = c("Small", "Medium", "Large"))
  
  ggplot(all_sizes_bh, aes(x = Trophic.Level, 
                           y = visregFit)) +
    geom_line(linewidth = 1, aes(color = Rarity)) +
    geom_ribbon(aes(ymin = visregLwr, 
                    ymax = visregUpr,
                    fill = Rarity),
                alpha = 0.5) +
    facet_grid(cols = vars(length),
               rows = vars(Rarity), 
               scales = "free_y", ) + 
    ylab("Effect size of medium protection") +
    xlab("Trophic level") +
    theme_bw() +
    geom_hline(yintercept = 1, 
               linetype = 2, 
               linewidth = 0.5) +
    # harrypotter::scale_fill_hp_d("Gryffindor", direction  = -1) +
    # harrypotter::scale_color_hp_d("Gryffindor", direction  = -1) +
    scale_color_manual(values = pal_custom) +
    scale_fill_manual(values = pal_custom) +
    theme(strip.text = element_text(size = 10, 
                                    face = "bold"),
          axis.text = element_text(size = 10),
          strip.background = element_blank(),
          panel.spacing = unit(2, "lines"))
  
  ggsave(here::here("plots_revised", "suppl_int_bh.jpeg"), 
         units = "px", 
         width = 2540, 
         height = 2000, 
         dpi = 300)
  
  dev.off()
  
    
}

plot_bl = function(){
  
  int_bl = visreg::visreg(gamma_biomass()[[3]],
                                type = "conditional",
                                xvar = "Trophic.Level",
                                scale = "response", 
                                by = "MaxLength",
                                nn = 505,
                                gg = F)$fit
  
  int_bl = int_bl |> 
    dplyr::mutate(Length = ifelse(MaxLength == 10, "Small", 
                                  ifelse(MaxLength == 24, "Medium", "Large")))
  
  int_bl$Length <- factor(int_bl$Length, levels = c("Small", "Medium", "Large"))
  
  # With conf. int
  
  ggplot(int_bl, aes(x = Trophic.Level, y = visregFit)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr),
                alpha = 0.5,
                fill = "#03045e") +
    facet_wrap(vars(Length)) + 
    ylab("Effect size on biomass of low protection") +
    xlab("Trophic.Level") +
    theme_bw() +
    geom_hline(yintercept = 1, linetype = 2, linewidth = 0.5) +
    theme(strip.text = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 10),
          strip.background = element_blank(),
          panel.spacing = unit(2, "lines"))
  
  ggsave(here::here("plots_revised", "suppl_int_bl.jpeg"), 
         units = "px", 
         width = 2540, 
         height = 2000, 
         dpi = 300)
  
  dev.off()
  
}