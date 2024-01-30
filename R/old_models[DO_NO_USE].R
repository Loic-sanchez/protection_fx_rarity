occu_full = function(){
  
  df_occu_full = as.data.frame(data_occu[[1]])
  mod_occu_full = glm(log_full ~ (Trophic.Level * log_length) * Rarity + Importance,
                      data = df_occu_full)
  mod_occu_full_AIC = MASS::stepAIC(mod_occu_full, direction = "both")
  
  return(mod_occu_full_AIC)
  
}

occu_high = function(){
  
  df_occu_high = data_occu[[2]]
  
  mod_occu_high = glm(log_high ~ (Trophic.Level * log_length) * Rarity + Importance,
                      data = df_occu_high)
  mod_occu_high_AIC = MASS::stepAIC(mod_occu_high, direction = "both")
  
  return(mod_occu_high_AIC)
  
}

occu_light = function(){
  
  df_occu_light = data_occu[[3]]
  
  mod_occu_light = glm(log_light ~ (Trophic.Level * log_length) * Rarity + Importance,
                      data = df_occu_light)
  mod_occu_light_AIC = MASS::stepAIC(mod_occu_light, direction = "both")
  
  return(mod_occu_light_AIC)
  
}

abun_full = function(){
  
  df_abun_full = data_abun[[1]]
  
  mod_abun_full = glm(log_full ~ (Trophic.Level * log_length) * Rarity + Importance,
                      data = df_abun_full)
  mod_abun_full_AIC = MASS::stepAIC(mod_abun_full, direction = "both")
  
  return(mod_abun_full_AIC)
  
}

abun_high = function(){
  
  df_abun_high = data_abun[[2]]
  
  mod_abun_high = glm(log_high ~ (Trophic.Level * log_length) * Rarity + Importance,
                      data = df_abun_high)
  mod_abun_high_AIC = MASS::stepAIC(mod_abun_high, direction = "both")
  
  return(mod_abun_high_AIC)
  
}

abun_light = function(){
  
  df_abun_light = data_abun[[3]]
  
  mod_abun_light = glm(log_light ~ (Trophic.Level * log_length) * Rarity + Importance,
                       data = df_abun_light)
  
  mod_abun_light_AIC = MASS::stepAIC(mod_abun_light, direction = "both")
  
  return(mod_abun_light_AIC)
  
}

biom_full = function(){
  
  df_biom_full = data_biom[[1]]
  
  mod_biom_full = glm(log_full ~ (Trophic.Level * log_length) * Rarity + Importance,
                      data = df_biom_full)
  mod_biom_full_AIC = MASS::stepAIC(mod_biom_full, direction = "both")
  
  return(mod_biom_full_AIC)
  
}

biom_high = function(){
  
  df_biom_high = data_biom[[2]]
  
  mod_biom_high = glm(log_high ~ (Trophic.Level * log_length) * Rarity + Importance,
                      data = df_biom_high)
  mod_biom_high_AIC = MASS::stepAIC(mod_biom_high, direction = "both")
  
  return(mod_biom_high_AIC)
  
}

biom_light = function(){
  
  df_biom_light = data_biom[[3]]
  
  mod_biom_light = glm(log_light ~ (Trophic.Level * log_length) * Rarity + Importance,
                       data = df_biom_light)
  mod_biom_light_AIC = MASS::stepAIC(mod_biom_light, direction = "both")
  
  return(mod_biom_light_AIC)
  
}

fx_occu_full = function(){
  
  interaction_occu_full = visreg::visreg(occu_full, 
                                         type = "conditional",
                                         "Trophic.Level",
                                         by = "Rarity",
                                         trans = exp, 
                                         gg = F,
                                         nn = 2000)$fit
  
  interaction_occu_full$group = "Occurrence"
  interaction_occu_full$Protection = "Full"
  
  ggplot(interaction_occu_full, aes(x = Trophic.Level, y = visregFit)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.5) +
    facet_wrap(~Rarity) + 
    ylab("Effect size on occurrence of full protection") +
    xlab("Trophic level") +
    theme_minimal() +
    theme(strip.text.x = element_text(size = 10, face = "bold"))
  
}

fx_abun_high = function(){
  
  interaction_ab_high = visreg::visreg(abun_high, 
                                       type = "conditional",
                                       "Trophic.Level",
                                       by = "Rarity",
                                       trans = exp, 
                                       gg = F,
                                       nn = 2000)$fit
  
  interaction_ab_high$group = "Abundance"
  interaction_ab_high$Protection = "High"
  
  ggplot(interaction_ab_high, aes(x = Trophic.Level, y = visregFit)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.5) +
    facet_wrap(~Rarity) + 
    ylab("Effect size on abundance of high protection") +
    xlab("Trophic level") +
    theme_minimal() +
    theme(strip.text.x = element_text(size = 10, face = "bold"))
}

fx_biom_full = function(){
  
  interaction_bm_full = visreg::visreg(biom_full, 
                                       type = "conditional",
                                       "log_length",
                                       by = "Trophic.Level",
                                       trans = exp, 
                                       gg = F,
                                       nn = 2000)$fit
  
  supp.labs <- c("2" = "Low","3.3" = "Intermediate", "4" = "High")
  interaction_bm_full$Trophic.Level = as.factor(interaction_bm_full$Trophic.Level)
  
  ggplot(interaction_bm_full, aes(x = log_length, y = visregFit)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), 
                alpha = 0.5, 
                fill = "#0d0887") +
    ylab("Effect size on biomass of full protection") +
    xlab("Maximum length (log)") +
    facet_grid(~Trophic.Level, labeller = labeller(Trophic.Level = supp.labs)) +
    theme_minimal() +
    theme(strip.text.x = element_text(size = 10, face = "bold"),
          axis.title.y = element_text(vjust = +3),
          axis.title.x = element_text(vjust = -0.75))
  
}

fx_biom_full2 = function(){ 
  
interaction_bm_full2 = visreg::visreg(biom_full, 
                                      type = "conditional",
                                      "Trophic.Level",
                                      by = "Rarity",
                                      trans = exp, 
                                      gg = F,
                                      nn = 2000)$fit

interaction_bm_full2$group = "Biomass"
interaction_bm_full2$Protection = "Full"

ggplot(interaction_bm_full2, aes(x = Trophic.Level, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.5) +
  ylab("Effect size on biomass of high protection") +
  xlab("Trophic level") +
  facet_grid(~Rarity) +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 10, face = "bold"))

}

all_rarity_trophic = function(){
  
  interaction_occu_high = visreg::visreg(occu_high,
                                         type = "conditional",
                                         "Trophic.Level",
                                         by = "Rarity",
                                         trans = exp, 
                                         gg = F,
                                         nn = 2000)$fit
  interaction_occu_high$group = "Occurrence"
  
  interaction_ab_high = visreg::visreg(abun_high, 
                                       type = "conditional",
                                       "Trophic.Level",
                                       by = "Rarity",
                                       trans = exp, 
                                       gg = F,
                                       nn = 2000)$fit
  interaction_ab_high$group = "Abundance"
  
  
  interaction_bm_full2 = visreg::visreg(biom_full, 
                                        type = "conditional",
                                        "Trophic.Level",
                                        by = "Rarity",
                                        trans = exp, 
                                        gg = F,
                                        nn = 2000)$fit
  interaction_bm_full2$group = "Biomass"
  
    
  all_rarity_trophic = rbind(interaction_occu_high %>% dplyr::select(-c(log_high)),
                             interaction_ab_high %>% dplyr::select(-c(log_high)),
                             interaction_bm_full2 %>% dplyr::select(-c(log_length, log_full)))
  
  all_rarity_trophic = all_rarity_trophic %>% 
    mutate(group = fct_relevel(group, "Occurrence", "Abundance", "Biomass"))
  
  ggplot(all_rarity_trophic, aes(x = Trophic.Level, y = visregFit, group = group)) +
    geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr, fill = group), 
                alpha = 0.75,
                show.legend = F) +
    geom_line(linewidth = 1,
              show.legend = F) +
    ylab("Effect size of protection") +
    xlab("Trophic level") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_grid(group~Rarity, scales = "free") +
    theme_minimal() +
    theme(strip.text = element_text(size = 12, face = "bold"),
          panel.spacing.x = unit(2.5, "lines"),
          axis.title.y = element_text(vjust = +2.65, 
                                      size = 15),
          axis.title.x = element_text(vjust = -0.75,
                                      size = 15),
          axis.text = element_text(size = 12)) +
    scale_fill_manual(values = c("#00b4d8", "#00b4d8", "#03045e")) +
    labs(fill = "Type")
}

threeway_interaction = function(){
  
  interaction_bm_high_rare = visreg::visreg(biom_high, 
                                            type = "conditional",
                                            "log_length",
                                            by = "Trophic.Level",
                                            cond = list(Rarity = "Rare"),
                                            trans = exp, 
                                            gg = F,
                                            nn = 2000)$fit
  interaction_bm_high_rare$Rarity = "Rare"
  
  interaction_bm_high_common = visreg::visreg(biom_high, 
                                              type = "conditional",
                                              "log_length",
                                              by = "Trophic.Level",
                                              cond = list(Rarity = "Common"),
                                              trans = exp, 
                                              gg = F,
                                              nn = 2000)$fit
  interaction_bm_high_common$Rarity = "Common"
  
  both_bm_high = rbind(interaction_bm_high_common, interaction_bm_high_rare)
  
  supp.labs <- c("2" = "Low trophic level","3.3" = "Intermediate trophic level", "4" = "High trophic level")
  
  ggplot(both_bm_high, aes(x = 10^log_length, y = visregFit)) +
    geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), 
                alpha = 0.75,
                fill = "#00b4d8") +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    ylab("Effect size on biomass of medium protection") +
    xlab("Maximum length (log)") +
    scale_x_log10() +
    facet_grid(Rarity ~ Trophic.Level, scales = 'free_y', labeller = labeller(Trophic.Level = supp.labs)) +
    theme_minimal() +
    theme(strip.text = element_text(size = 10, face = "bold"),
          panel.spacing = unit(0.5, "lines"),
          axis.title.y = element_text(vjust = +3),
          axis.title.x = element_text(vjust = -0.75))

}