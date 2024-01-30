make_ssmatrix_melted = function(n) {
  
  tab <- table(fish_data$valid_name_FishBase)
  boxplot(tab, outline = F)
  tab <- tab[tab>n]
  fish_filt  <- fish_data[fish_data$valid_name_FishBase %in% names(tab),]
  fish_filt$Species = as.factor(fish_filt$valid_name_FishBase)
  
  ssmatrix = t(fossil::create.matrix(fish_filt, 
                                     tax.name = "Species", 
                                     locality = "SurveyID",
                                     abund = F))
  
  melted = ssmatrix |>
    reshape2::melt() |>
    rename(SurveyID = "Var1", 
           Species = "Var2", 
           Presence = "value") |>
    left_join(covariates_corrected |> 
                dplyr::select(SurveyID, 
                              Effectiveness, 
                              DiffYear,
                              Area_km2), by = "SurveyID") |>
    left_join(sites_info[,c(1:4,10)], by = "SurveyID") |>
    drop_na(Effectiveness)
  
  melted$Effectiveness = fct_collapse(melted$Effectiveness, 
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
  
  melted$Effectiveness = as.character(melted$Effectiveness)
  melted$Effectiveness[melted$Effectiveness == "Highly Protected" & melted$DiffYear >= 10 & melted$Area_km2 >= 100] = "Fully Protected"
  melted$Effectiveness = as.factor(melted$Effectiveness)
  melted$Effectiveness[melted$Effectiveness != "Unprotected" & melted$DiffYear < 1] = "Unprotected"
  
  return(melted)
  
}

rid_rare_outside = function() { # Get rid of species with 0-4 occurrences in the reference level
  
  # Some values were dropped, make sure we still have at least N occurrences
  
  melted_sum = melted |>
    group_by(Species) |>
    summarise(Presence = sum(Presence)) |>
    dplyr::filter(Presence >= 30)
  
  spec_vec = as.character(melted_sum$Species)
  
  melted_filt = melted[melted$Species %in% spec_vec,]
  
  # Get rid of rares outside
  
  melted_fsum = melted_filt |>
    filter(Effectiveness == "Unprotected") |>
    group_by(Species) |>
    summarise(Presence = sum(Presence)) |>
    filter(Presence > 4) 
  
  spec_vec = as.character(melted_fsum$Species)
  
  melted_filt = melted_filt[melted_filt$Species %in% spec_vec,]
  melted_filt$Species = as.factor(as.character(melted_filt$Species))
  
  return(melted_filt)
}

add_covariates = function(){
  
  all_cov = hab_filt |>
    left_join(socio_filt, by = "SurveyID") |>
    left_join(env_filt, by = "SurveyID") |>
    left_join(fine_habitat[,-c(2:4)], by = "SurveyID")
  
  melted_cov = melted_filt |> 
    left_join(all_cov, by = "SurveyID")
  
  return(melted_cov)
  
}