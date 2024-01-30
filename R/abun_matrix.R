make_ssmatrix_melted_ab = function(n) {
  
  tab <- table(fish_data$valid_name_FishBase)
  boxplot(tab, outline = F)
  tab <- tab[tab>n]
  fish_filt  <- fish_data[fish_data$valid_name_FishBase %in% names(tab),]
  fish_filt$Species = as.factor(fish_filt$valid_name_FishBase)
  
  fish_filt1 = fish_filt |>
    group_by(Species) |>
    mutate(percentile_99 = ifelse(max(Num) > 499, quantile(Num, 0.99), quantile(Num, 1))) |>
    filter(Num < percentile_99) |>
    ungroup()
  
  ssmatrix_ab = t(fossil::create.matrix(as.data.frame(fish_filt1), 
                                        tax.name = "Species", 
                                        locality = "SurveyID",
                                        abund = T,
                                        abund.col = "Num"))
  
  melted_ab = ssmatrix_ab |>
    reshape2::melt() |>
    rename(SurveyID = "Var1", 
           Species = "Var2", 
           Abundance = "value") |>
    left_join(covariates_corrected |> 
                dplyr::select(SurveyID, 
                              Effectiveness, 
                              DiffYear,
                              Area_km2), by = "SurveyID") |>
    left_join(sites_info[,c(1:4,10)], by = "SurveyID") |>
    drop_na(Effectiveness)
  
  melted_ab$Effectiveness = fct_collapse(melted_ab$Effectiveness, 
                                         "Highly Protected" = c("High No take", 
                                                                "High No take multizoned",
                                                                "Medium No take", 
                                                                "Medium No take multizoned"),
                                         "Lightly Protected" = c("High Restricted take", 
                                                                 "High Restricted take multizoned", 
                                                                 "Medium Restricted take",
                                                                 "Medium Restricted take multizoned",
                                                                 "Low No take", 
                                                                 "Low No take multizoned",
                                                                 "Low Restricted take",
                                                                 "Low Restricted take multizoned"),
                                         "Unprotected" = c("Medium Fishing", "Low Fishing", "out"))
  
  melted_ab$Effectiveness = as.character(melted_ab$Effectiveness)
  melted_ab$Effectiveness[melted_ab$Effectiveness == "Highly Protected" & melted_ab$DiffYear >= 10 & melted_ab$Area_km2 >= 100] = "Fully Protected"
  melted_ab$Effectiveness = as.factor(melted_ab$Effectiveness)
  melted_ab$Effectiveness[melted_ab$Effectiveness != "Unprotected" & melted_ab$DiffYear < 1] = "Unprotected"
  
  return(melted_ab)
  
}

rid_rare_outside_ab = function() { # Get rid of species with 0-4 occurrences in the reference level
  
  # Some values were dropped, make sure we still have at least N occurrences
  
  melted_sum = melted_ab |>
    filter(Abundance > 0) |>
    group_by(Species) |>
    summarise(Abundance = n()) |>
    filter(Abundance >= 30)
  
  spec_vec = as.character(melted_sum$Species)
  
  melted_filt_ab = melted_ab[melted_ab$Species %in% spec_vec,]
  
  # Get rid of rares outside
  
  melted_fsum = melted_filt_ab |>
    filter(Effectiveness == "Unprotected" & Abundance > 0) |>
    group_by(Species) |>
    summarise(Abundance = n()) |>
    filter(Abundance > 4) 
  
  spec_vec = as.character(melted_fsum$Species)
  
  melted_filt_ab = melted_filt_ab[melted_filt_ab$Species %in% spec_vec,]
  melted_filt_ab$Species = as.factor(as.character(melted_filt_ab$Species))
  
  return(melted_filt_ab)
  
}

add_covariates_ab = function(){
  
  all_cov = hab_filt |>
    left_join(socio_filt, by = "SurveyID") |>
    left_join(env_filt, by = "SurveyID") |>
    left_join(fine_habitat[,-c(2:4)], by = "SurveyID")
  
  melted_cov = melted_filt_ab |> 
    left_join(all_cov, by = "SurveyID")
  
  return(melted_cov)
  
}