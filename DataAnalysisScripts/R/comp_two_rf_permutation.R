source(here::here("DataAnalysisScripts","R", "pseudo_log_breaks.R"))

comp_two_rf_permutation <- function(all, 
                                    vars_simple, 
                                    vars, 
                                    gas_name, 
                                    layer_names,
                                    log_vars,
                                    reps = 100){
  
  vars <- vars[!vars %in% c("LakeID", "Latitude", "Longitude")]
  colors <- c("#40476D", "#1098F7", "#0C7C59", "gray70")
  names(colors) <- vars_simple
  
  log_select <- log_vars[log_vars %in% c(vars, vars_simple)]
  
  ### SET UP ###
  focal <- all %>% 
    filter(if_all(all_of(log_select), \(x) x>0)
    ) %>%
    ungroup() %>%
    mutate(across(all_of(log_select), log),
           mutate(across(where(is.character), as.factor))) %>%
    select(all_of(c("value", vars_simple, vars, "LakeID", "name", "Layer")))
  
  for_rf <- focal %>%
    na.omit()
  
  ns <- for_rf %>%
    group_by(name, Layer) %>%
    summarize(n = length(unique(LakeID))) %>%
    rename(gas = name)
  
  ### RUN RANDOM FOREST ###
  rf_out <- lapply(1:reps, function(x) {
    bot_ch4 <- permutation_rf_comp(x, for_rf, vars_simple, vars, 
                                   gas_name = "CH[4]", layer_names = "bot")
    bot_co2 <- permutation_rf_comp(x, for_rf, vars_simple, vars, 
                                   gas_name = "CO[2]", layer_names = "bot")
    surf_ch4 <- permutation_rf_comp(x, for_rf, vars_simple, vars, 
                                   gas_name = "CH[4]", layer_names = "surf")
    surf_co2 <- permutation_rf_comp(x, for_rf, vars_simple, vars, 
                                   gas_name = "CO[2]", layer_names = "surf")
    
    eval <- list(bot_ch4[[1]],
                 bot_co2[[1]],
                 surf_ch4[[1]],
                 surf_co2[[1]]) %>%
      bind_rows()
    
    ImpData <- list(bot_ch4[[2]],
                    bot_co2[[2]],
                    surf_ch4[[2]],
                    surf_co2[[2]]) %>%
      bind_rows()
    
    list(eval, ImpData)
  })

  #Unpack
  eval <- lapply(rf_out, "[[", 1) %>%
    bind_rows()
  
  eval_join <- eval %>%
    group_by(gas, Layer) %>%
    summarize(r2 = mean(r2_simple)) 
  
  ImpData <- bind_rows(sapply(rf_out, "[", 2))
  
  ### PLOT ###
  # Visualize variable importance ----------------------------------------------
  var_names_imp = sub("\n", "~", var_names)
  imp_levels_used <- rev(recode(vars_simple, !!!var_names_imp))
  names(colors) <- var_names_imp[vars_simple]
  
  imp <- ImpData %>%
    left_join(ns) %>%
    left_join(eval_join) %>%
    mutate(Var.Names = recode(Var.Names, !!!var_names_imp),
           Var.Names = factor(Var.Names, levels = imp_levels_used),
           Layer = case_match(Layer,
                              "surf"~"Surface",
                              "bot"~"Bottom")) %>%
    group_by(Var.Names, gas, Layer) %>%
    summarize(sd = sd(`%IncMSE`, na.rm = T),
              `%IncMSE` = mean(`%IncMSE`, na.rm = T),
              title = paste0("'", unique(Layer), "'~", unique(gas), 
                             "~`(`*italic(n)*` =`~", unique(n), 
                             "*'; '*mean~R^2*' = '~", round(unique(r2),2), 
                             "*')'"),
              .groups = "drop") %>%
    mutate(fct_order = case_when(gas == "CH[4]" & Layer == "Surface"~1,
                             gas == "CO[2]" & Layer == "Surface"~2,
                             gas == "CH[4]" & Layer == "Bottom"~3,
                             gas == "CO[2]" & Layer == "Bottom"~4)) %>%
    mutate(title = factor(title),
           title = fct_reorder(title, fct_order)) %>%
    ggplot(aes(x = `%IncMSE`, y = Var.Names, fill = Var.Names))+
    geom_col()+
    geom_errorbar(aes(xmin = `%IncMSE`-sd, xmax = `%IncMSE`+sd), 
                  width = 0.1)+
    scale_fill_manual(values = colors)+
    xlab("Variable importance (% increase in MSE)")+
    egg::theme_article()+
    theme(axis.title.y = element_blank(),
          plot.title.position = "plot",
          plot.title = element_text(hjust = 0))+
    scale_x_continuous(n.breaks = 3)+
    scale_y_discrete(labels = parse(text = imp_levels_used))+
    facet_wrap(~title, labeller = "label_parsed")+
    theme(legend.position = "none",
          panel.grid.major = element_line(color = "grey90", size = 0.5),
          panel.grid.minor = element_line(color = "grey95", size = 0.25))
  
  eval_sum <- eval %>%
    pivot_longer(c(rmse_simple, rmse_complex), names_prefix = "rmse_",
                 values_to = "rmse") %>%
    group_by(gas, Layer, name) %>%
    summarize(mean = mean(rmse),
              sd = sd(rmse))
  
  min <- min(eval_sum$mean - eval_sum$sd)
  max <- max(eval_sum$mean + eval_sum$sd)
  
  rmse <- eval_sum %>%
    mutate(Layer = case_match(Layer,
                              "surf"~"Surface~",
                              "bot"~"Bottom~"),
           label = paste0(Layer, gas)) %>%
    pivot_wider(names_from = name, values_from = c(mean, sd)) %>%
    ggplot(aes(x = mean_simple, y = mean_complex))+
    geom_point()+
    geom_errorbar(aes(ymin = mean_complex - sd_complex, 
                      ymax = mean_complex + sd_complex),
                  width = 0.1)+
    geom_errorbar(aes(xmin = mean_simple - sd_simple, 
                      xmax = mean_simple + sd_simple),
                  width = 0.1)+
    geom_abline(color= "grey50")+
    xlab("RMSE without in-lake\ndriver data (µM)")+
    ylab("RMSE with in-lake\ndriver data (µM)")+
    ggrepel::geom_label_repel(aes(label = label), parse = T, 
                              size = 3, label.padding = 0.08, 
                              force_pull = 0.2,
                              direction = "both")+
    xlim(c(min,max))+
    ylim(c(min,max))+
    theme_bw()+
    theme(aspect.ratio = 1)
  
  eval_paired <- eval %>%
    mutate(pct_decrease = (rmse_complex - rmse_simple)/rmse_simple*100) %>%
    group_by(gas, Layer) %>%
    summarize(mean = mean(pct_decrease),
              sd = sd(pct_decrease))
  
  bar <- eval_paired %>%
    mutate(Layer = case_match(Layer,
                              "surf"~"Surface",
                              "bot"~"Bottom"),
           Layer = factor(Layer, levels = c("Surface", "Bottom"))) %>%
    ggplot(aes(y = mean, x = Layer))+
    ylab("Percentage decrease in RMSE\nusing in-lake driver data")+
    egg::theme_article()+
    geom_hline(color= "grey60", yintercept = 0)+
    geom_col(width = 0.8)+
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1)+
    facet_wrap(~gas, labeller = label_parsed)+
    theme(axis.title.x = element_blank(),
          panel.grid.major = element_line(color = "grey90", size = 0.5),
          panel.grid.minor = element_line(color = "grey95", size = 0.25),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  
  s<-ggpubr::ggarrange(imp, 
                       ggpubr::ggarrange(rmse, bar,
                                         labels = c("b", "c"),
                                         label.x = c(0, -0.05)), 
                       nrow = 2, labels = c("a", NA))
  return(s)
}

permutation_rf_comp <- function(rep_n, for_rf, vars_simple, 
                                vars, gas_name, layer_names){
  for_rf_lim <- for_rf %>%
    filter(name == gas_name,
           Layer %in% layer_names) %>%
    group_by(LakeID) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    select(-LakeID, -Layer, -name)
  
  for_rf_simple <- for_rf_lim %>%
    select(all_of(c("value", vars_simple)))
  
  for_rf_complex <- for_rf_lim %>%
    select(all_of(c("value", vars)))
  
  rf_simple <- randomForest(value ~ ., data = for_rf_simple, 
                          proximity=TRUE, importance = T) 
  
  rf_complex <- randomForest(value ~ ., data = for_rf_complex, 
                      proximity=TRUE, importance = T) 
  
  ImpData_simple <- as.data.frame(randomForest::importance(rf_simple)) %>% 
    mutate(rep = rep_n,
           gas = gas_name,
           Layer = layer_names
           )
  ImpData_simple$Var.Names <- row.names(ImpData_simple)
  
  preds <- data.frame(rf_simple_pred = predict(rf_simple),
                      rf_complex_pred = predict(rf_complex),
                      obs = for_rf_lim$value)
  
  r2_simple <- rf_simple$rsq[length(rf_simple$rsq)]
  
  out <- preds %>%
    summarize(rmse_simple = Metrics::rmse(obs, rf_simple_pred),
              rmse_complex = Metrics::rmse(obs, rf_complex_pred),
              mae_simple = Metrics::mae(obs, rf_simple_pred),
              mae_complex = Metrics::mae(obs, rf_complex_pred),
              median = median(abs(obs)),
              mean = mean(abs(obs)),
              r2_simple = r2_simple,
              n = n(),
              gas = gas_name,
              Layer = layer_names) %>%
    mutate(across(rmse_simple:mean, exp))
  
  return(list(out, ImpData_simple))
}

