all_sp_glm_qual = lapply(1:length(unique(melted_cov$Species)), function(i) {
# all_sp_glm_qual = lapply(1:100, function(i) {
    
  subdf = filter(melted_cov, Species == unique(melted_cov$Species)[i])
  
  pres_sp1 = subdf %>%
    filter(Presence > 0)
  
  pres_sp1$Ecoregion = as.factor(as.character(pres_sp1$Ecoregion))
  
  eco_table = table(pres_sp1$Ecoregion)
  eco_tab = names(eco_table[eco_table > 4])
  
  abs_sp1 = subdf %>%
    filter(Presence < 1)
  
  abs_sp1 = abs_sp1[abs_sp1$Ecoregion %in% eco_tab,]
  pres_sp1 = pres_sp1[pres_sp1$Ecoregion %in% eco_tab,]
  
  presabs_sp1 = rbind(pres_sp1, abs_sp1)
  
  eco_weights = presabs_sp1 |>
    group_by(Ecoregion) |> 
    mutate(weights = (1/n())) |> 
    ungroup()
  
  pca_prep = all_cov[all_cov$SurveyID %in% unique(presabs_sp1$SurveyID),]
  pca_prep = pca_prep |> 
    left_join(eco_weights %>% dplyr::select(SurveyID, weights), by = "SurveyID")
  
  pca = ade4::dudi.pca(pca_prep[, -c(1, 97)],
                       center = T, scale = T, 
                       row.w = pca_prep$weights,
                       nf = 25, 
                       scannf = F)
  
  PC_coords = pca$li
  
  eigs = factoextra::get_eig(pca)
  cum_var = which(eigs[,3] > 70)[1]
  
  inertia_ax1 = eigs[1, 2]
  inertia_ax2 = eigs[2, 2]
  
  pca_prep_and_coords = cbind(pca_prep, PC_coords)
  presabs_sp1 = presabs_sp1[, 1:10] |> 
    left_join(pca_prep_and_coords[, c(1, 98:122)], by = "SurveyID") 
  
  istherefull = ifelse(sum(presabs_sp1$Presence[which(presabs_sp1$Effectiveness == "Fully Protected")]) > 4, T, F)
  istherehigh = ifelse(sum(presabs_sp1$Presence[which(presabs_sp1$Effectiveness == "Highly Protected")]) > 4, T, F)
  istherelight = ifelse(sum(presabs_sp1$Presence[which(presabs_sp1$Effectiveness == "Lightly Protected")]) > 4, T, F)
  isthereunp = ifelse(sum(presabs_sp1$Presence[which(presabs_sp1$Effectiveness == "Unprotected")]) > 4, T, F)
  
  presabs_sp1_nofull = presabs_sp1 %>%
    filter(Effectiveness != "Fully Protected")
  presabs_sp1_nohigh = presabs_sp1 %>%
    filter(Effectiveness != "Highly Protected")
  presabs_sp1_nolight = presabs_sp1 %>% 
    filter(Effectiveness != "Lightly Protected")
  presabs_sp1_nofullhigh = presabs_sp1 %>% 
    filter(Effectiveness != "Highly Protected" & Effectiveness != "Fully Protected")
  presabs_sp1_nofulllight = presabs_sp1 %>% 
    filter(Effectiveness != "Lightly Protected" & Effectiveness != "Fully Protected")
  presabs_sp1_nohighlight = presabs_sp1 %>% 
    filter(Effectiveness != "Lightly Protected" & Effectiveness != "Highly Protected")
  
  if(istherefull == F){presabs_sp1 = presabs_sp1_nofull}
  if(istherehigh == F){presabs_sp1 = presabs_sp1_nohigh}
  if(istherelight == F){presabs_sp1 = presabs_sp1_nolight}
  if(istherelight == F & istherefull == F){presabs_sp1 = presabs_sp1_nofulllight}
  if(istherehigh == F & istherefull == F){presabs_sp1 = presabs_sp1_nofullhigh}
  if(istherelight == F & istherehigh == F){presabs_sp1 = presabs_sp1_nohighlight}
  
  if(isthereunp == F){return(NULL)} 
  if(istherefull == F & istherehigh == F & istherelight == F){return(NULL)}
  
  names_var = c("Effectiveness", names(as.data.frame(PC_coords[, c(1:cum_var)])))
  formula_var = paste0(rev(names_var), collapse = "+")
  full_formula = paste("Presence ~ ", formula_var, sep = "")
  
  model = glm(formula = full_formula,
              family = "binomial",
              data = presabs_sp1)
  
  if(model$converged == F){return(NULL)}
  
  sum_isthere = sum(istherefull, istherehigh, istherelight)
  if(sum_isthere < 2){vif = sqrt((regclass::VIF(model)))} else{vif = regclass::VIF(model)[, 3]}
  
  if(max(vif) < sqrt(5)){no_multico = TRUE} else{no_multico = FALSE}
  # no_multico = as.data.frame(no_multico)
  
  # names_var = rev(names_var)
  
  while(max(vif) > sqrt(5)){ 
    
    u = which.max(vif)
    # u = ifelse(which.max(vif) == which(stringr::str_detect(names(vif), "Axis1")), 
    #            kit::topn(vif, 
    #                      2, 
    #                      decreasing = T)[2], 
    #            u) 
    names_var = names_var[-u]
    formula_var = paste0(names_var, collapse = "+")
    full_formula = paste("Presence ~ ", formula_var, sep = "")
    
    model = glm(formula = full_formula,
                family = "binomial",
                data = presabs_sp1)            
    
    # if(sum(str_detect(names_var, "Effectiveness")) == 0 | length(names_var) < 2){return(NULL)}
    if(sum_isthere < 2){vif = sqrt((regclass::VIF(model)))} else{vif = regclass::VIF(model)[, 3]}
    
    full_formula = full_formula
    
  }
  
  effectiveness_miss = stringr::str_detect(full_formula, "Axis1")
  
  
  
  cat(i, "\n")
  return(df_count = data.frame(no_multico, 
                               effectiveness_miss,
                               Species = unique(melted_cov$Species)[i]))
})

df_count_axes = data.table::rbindlist(all_sp_glm_qual)
df_count_axes_rev = data.table::rbindlist(all_sp_glm_qual)
sum(df_count_axes$effectiveness_miss)
sum(df_count_axes$axis1_count)
  
