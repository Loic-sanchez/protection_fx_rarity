add_rarity = function(){ 

melted_rares = melted_cov %>%
  filter(Presence > 0) %>% 
  group_by(Species) %>% 
  summarise(n_occu = n()) %>% 
  mutate(decile_rank = ntile(n_occu, 10)) %>% 
  # mutate(Rarity = ifelse(decile_rank > 5, "Common", "Rare"))
  mutate(Rarity = ifelse(n_occu > 100, "Common", "Rare"))

ggplot(melted_rares, aes(x = n_occu, colour = Rarity)) + 
  geom_histogram(binwidth = 5) + 
  xlim(0, 500)

return(melted_rares)

}

# targets = c("Acanthuridae", "Caesionidae", "Carangidae", 
#             "Ephippidae", "Haemulidae", "Kyphosidae", 
#             "Labridae", "Lethrinidae", "Lutjanidae", 
#             "Mullidae", "Nemipteridae", "Scaridae", 
#             "Scombridae", "Serranidae", "Siganidae", 
#             "Sparidae", "Sphyraenidae")
# 
# others_20cm = c("Balistidae", "Holocentridae", "Pomacanthidae",
#                 "Priacanthidae") 
# 
# traits$target[traits$Family %in% targets] = "target"
# traits$target[traits$Family %in% others_20cm & traits$MaxLength > 20] = "target"
# traits$target = replace(traits$target, is.na(traits$target), "not target")

traits_fish = rfishbase::species(unique(traits$CURRENT_SPECIES_NAME))
traits_fish_sel = traits_fish |> dplyr::select(Species, Importance, Vulnerability)

traits_and_imp = traits |> left_join(traits_fish_sel, by = "Species") |> 
  filter(Species != "Hypoplectrus aberrans")

traits_fish_sel$Importance[traits_fish_sel$Importance == "commercial" |
                             traits_fish_sel$Importance == "highly commercial" |
                             traits_fish_sel$Importance == "minor commercial"] = "Commercial"
traits_fish_sel$Importance[traits_fish_sel$Importance == "of no interest" |
                             traits_fish_sel$Importance == "subsistence fisheries" |
                             traits_fish_sel$Importance == "of potential interest"] = "Not_commercial"
traits_fish_sel$Importance = replace(traits_fish_sel$Importance, is.na(traits_fish_sel$Importance), "Not_commercial")

# traits_fish_sel$Importance[traits_fish_sel$Importance == "commercial"] = 3
# traits_fish_sel$Importance[traits_fish_sel$Importance == "highly commercial"] = 4
# traits_fish_sel$Importance[traits_fish_sel$Importance == "minor commercial"] = 2
# traits_fish_sel$Importance[traits_fish_sel$Importance == "of no interest" |
#                              traits_fish_sel$Importance == "of potential interest"] = 0
# traits_fish_sel$Importance[traits_fish_sel$Importance == "subsistence fisheries"] = 1
# traits_fish_sel$Importance = replace(traits_fish_sel$Importance, is.na(traits_fish_sel$Importance), 0)
# traits_fish_sel$Importance = as.numeric(traits_fish_sel$Importance)

prep_data_occu = function() {
  
  traits$Species = traits$CURRENT_SPECIES_NAME
  
  df_occu_traits = OS_occu %>% 
    mutate(log_full = log(RR_full),
           log_high = log(RR_high),
           log_light = log(RR_light)) %>% 
    left_join(traits, by = "Species") %>% 
    left_join(traits_fish_sel, by = "Species") |> 
    left_join(melted_rares %>% dplyr::select(Species, Rarity, decile_rank), by = "Species")  %>% 
    mutate(log_length = log10(MaxLength)) 
  
  df_occu_full = df_occu_traits %>% 
    drop_na(log_full, Trophic.Level, log_length, Importance)
  df_occu_high = df_occu_traits %>% 
    drop_na(log_high, Trophic.Level, log_length, Importance)
  df_occu_light = df_occu_traits %>% 
    drop_na(log_light, Trophic.Level, log_length, Importance)
  
  list_occu = list(df_occu_full, df_occu_high, df_occu_light)
}

prep_data_abun = function(){
  
  traits$Species = traits$CURRENT_SPECIES_NAME
  
  df_ab_traits = df_ab %>%
    mutate(log_full = log(IRR_full),
           log_high = log(IRR_high),
           log_light = log(IRR_light)) %>% 
    left_join(traits, by = "Species") %>% 
    left_join(traits_fish_sel, by = "Species") |> 
    left_join(melted_rares %>% dplyr::select(Species, Rarity, decile_rank), by = "Species")  %>% 
    mutate(log_length = log10(MaxLength))
  
  df_ab_full = df_ab_traits %>% 
    drop_na(log_full, Trophic.Level, log_length, Importance)
  df_ab_high = df_ab_traits %>% 
    drop_na(log_high, Trophic.Level, log_length, Importance)
  df_ab_light = df_ab_traits %>% 
    drop_na(log_light, Trophic.Level, log_length, Importance) 
  
  list_abun = list(df_ab_full, df_ab_high, df_ab_light)
  
}

prep_data_bm = function() {
  
  traits$Species = traits$CURRENT_SPECIES_NAME
  
  df_bm_product = ES_bm %>% 
    left_join(ES_occu %>% 
                dplyr::select(Species, 
                              estimate_outside, 
                              estimate_full, 
                              estimate_high,
                              estimate_light, 
                              AUC), 
              by = "Species") %>% 
    rowwise() %>% 
    mutate(product_full = estimate_full.x*estimate_full.y,
           product_high = estimate_high.x*estimate_high.y,
           product_light = estimate_light.x*estimate_light.y,
           product_outside = estimate_outside.x*estimate_outside.y) %>% 
    mutate(IRR_product_full = product_full/product_outside,
           IRR_product_high = product_high/product_outside,
           IRR_product_light = product_light/product_outside) %>% 
    drop_na(estimate_outside.y)
  
  df_bm_traits = df_bm_product %>%
    mutate(log_full = log(IRR_product_full),
           log_high = log(IRR_product_high),
           log_light = log(IRR_product_light)) %>% 
    left_join(traits, by = "Species") %>% 
    left_join(traits_fish_sel, by = "Species") |> 
    left_join(melted_rares %>% dplyr::select(Species, Rarity, decile_rank), by = "Species")  %>% 
    mutate(log_length = log10(MaxLength))
  
  df_bm_full = df_bm_traits %>% 
    drop_na(log_full, Trophic.Level, log_length, Importance)  
  df_bm_high = df_bm_traits %>% 
    drop_na(log_high, Trophic.Level, log_length, Importance)
  df_bm_light = df_bm_traits %>% 
    drop_na(log_light, Trophic.Level, log_length, Importance)
  
  list_bm = list(df_bm_full, df_bm_high, df_bm_light)
  
}

boxplot(df_occu_full$log_full[df_occu_full$Rarity == "Rare" & df_occu_full$Trophic.Level < 3.3] ~ df_occu_full$Importance[df_occu_full$Rarity == "Rare" & df_occu_full$Trophic.Level < 3.3],
        xlab = "", ylab = "", main = "Occurrence - High")
boxplot(df_ab_full$log_full[df_ab_full$Rarity == "Rare" & df_ab_full$Trophic.Level < 3.3] ~ df_ab_full$Importance[df_ab_full$Rarity == "Rare" & df_ab_full$Trophic.Level < 3.3],
        xlab = "", ylab = "", main = "Abundance - High")
boxplot(df_bm_full$log_full[df_bm_full$Rarity == "Rare" & df_bm_full$Trophic.Level < 3.3] ~ df_bm_full$Importance[df_bm_full$Rarity == "Rare" & df_bm_full$Trophic.Level < 3.3],
        xlab = "", ylab = "", main = "Biomass - High")

boxplot(df_occu_high$log_high[df_occu_high$Rarity == "Rare" & df_occu_high$Trophic.Level < 3.3] ~ df_occu_high$Importance[df_occu_high$Rarity == "Rare" & df_occu_high$Trophic.Level < 3.3],
        xlab = "", ylab = "Log-transformed Effect Size", main = "Occurrence - Medium")
boxplot(df_ab_high$log_high[df_ab_high$Rarity == "Rare" & df_ab_high$Trophic.Level < 3.3] ~ df_ab_high$Importance[df_ab_high$Rarity == "Rare" & df_ab_high$Trophic.Level < 3.3],
        xlab = "", ylab = "", main = "Abundance - Medium")
boxplot(df_bm_high$log_high[df_bm_high$Rarity == "Rare" & df_bm_high$Trophic.Level < 3.3] ~ df_bm_high$Importance[df_bm_high$Rarity == "Rare" & df_bm_high$Trophic.Level < 3.3],
        xlab = "", ylab = "", main = "Biomass - Medium")

boxplot(df_occu_light$log_light[df_occu_light$Rarity == "Rare" & df_occu_light$Trophic.Level < 3.3] ~ df_occu_light$Importance[df_occu_light$Rarity == "Rare" & df_occu_light$Trophic.Level < 3.3],
        xlab = "", ylab = "", main = "Occurrence - Low")
boxplot(df_ab_light$log_light[df_ab_light$Rarity == "Rare" & df_ab_light$Trophic.Level < 3.3] ~ df_ab_light$Importance[df_ab_light$Rarity == "Rare" & df_ab_light$Trophic.Level < 3.3],
        xlab = "Importance of rare species", ylab = "", main = "Abundance - Low")
boxplot(df_bm_light$log_light[df_bm_light$Rarity == "Rare" & df_bm_light$Trophic.Level < 3.3] ~ df_bm_light$Importance[df_bm_light$Rarity == "Rare" & df_bm_light$Trophic.Level < 3.3],
        xlab = "", ylab = "", main = "Biomass - Low")

par(mfrow = c(3,3))

