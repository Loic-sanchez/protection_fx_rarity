tidy_sites = function() {
  
  sites_correct$Effectiveness_cor = "out"  #Set to outside all sites that were unprotected when surveyed
  
  covariates_corrected = dplyr::left_join(covariates, 
                                          sites_correct[,c(1,21)], 
                                          by = "SurveyID")
  
  covariates_corrected$Effectiveness = replace(covariates_corrected$Effectiveness, 
                                               covariates_corrected$Effectiveness_cor == "out",
                                               "out") # Correct it on the original covariates data
  
  covar_jb_graham = covar_jb |> 
    left_join(covar_graham |> dplyr::select(SiteCode, area..km2.)) |>
    rename(Area_km2 = "area..km2.")
  
  covar_jb_graham$Area_km2 = as.numeric(covar_jb_graham$Area_km2)
  
  covar_jb_select = covar_jb_graham %>% select(SurveyID, DiffYear, Area_km2)
  covar_jb_sites = covar_jb_select[covar_jb_select$SurveyID %in% covariates_corrected$SurveyID,]
  
  covariates_corrected = covariates_corrected %>% 
    left_join(covar_jb_sites, by = "SurveyID")
  covariates_corrected = covariates_corrected[-which(duplicated(covariates_corrected$SurveyID)),]
  
  return(covariates_corrected)
  
}

all_cov_matrix = function()  {
  
  all_cov = hab_filt |>
    left_join(socio_filt, by = "SurveyID") |>
    left_join(env_filt, by = "SurveyID") |>
    left_join(fine_habitat[,-c(2:4)], by = "SurveyID")
  
  save(all_cov, file = here::here("data", "data", "all_cov.RData"))
}
