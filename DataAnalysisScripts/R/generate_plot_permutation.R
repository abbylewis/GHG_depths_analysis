source(here::here("DataAnalysisScripts", "R", "pseudo_log_breaks.R"))
source(here::here("DataAnalysisScripts", "R", "permutation_rf.R"))

#' Plot drivers of gas concentrations or saturation
#'
#' @param all Dataframe of focal data and drivers
#' @param vars Driver variables
#' @param var_names Variable names for figure
#' @param gas_name Which gas to analyze ('CH[4]' or 'CO[2]')
#' @param log_vars Which variables should be log transformed
#' @param pSat Are we plotting gas saturation? (changes axis label)
#' @param reps How many replications (permutation approach)
#'
#' @returns Drivers figure
#'
generate_plot_permutation <- function(all, vars, var_names, gas_name, log_vars,
                                      pSat = F, reps = 100) {
  vars <- vars[!vars %in% c("LakeID", "Latitude", "Longitude")]
  colors <- c("#69140E", "#40476D", "#1098F7", "#0C7C59", "#E65F5C", "gray70", "grey70", "grey70")
  names(colors) <- vars

  var_names_part <- c(
    "DO_mgL" = "DO~(mg~L^-1)",
    "SurfaceArea_km2" = "SA~(km^2)",
    "MaximumDepth_m" = "Max~depth~(m)",
    "MeanDepth_m" = "Mean~depth~(m)",
    "TP_ugL_mean" = "TP~(µg~L^-1)",
    "Temp_C" = "Temp~(ºC)",
    "Osgood" = "Osgood",
    "Absolute_latitude" = "Abs~lat~(DD)",
    "buoyancy_frequency" = "Buoy~freq~(s^-2)",
    "daylength" = "Day~length~(h)"
  )

  levels_used <- recode(vars, !!!var_names_part)
  log_select <- log_vars[log_vars %in% vars]

  ### SET UP ###
  for_rf_surf <- all %>%
    filter(
      name == gas_name,
      Layer == "surf",
      if_all(all_of(log_select), \(x) x > 0)
    ) %>%
    ungroup() %>%
    mutate(
      across(all_of(log_select), log),
      mutate(across(where(is.character), as.factor))
    ) %>%
    select(all_of(c("value", vars, "LakeID"))) %>%
    na.omit()

  for_rf_bot <- all %>%
    filter(
      name == gas_name,
      Layer == "bot",
      if_all(all_of(log_select), \(x) x > 0)
    ) %>%
    ungroup() %>%
    mutate(
      across(all_of(log_select), log),
      mutate(across(where(is.character), as.factor))
    ) %>%
    select(all_of(c("value", vars, "LakeID"))) %>%
    na.omit()

  num_vars <- vars[sapply(for_rf_surf[vars], is.numeric)]
  cat_vars <- vars[sapply(for_rf_surf[vars], is.factor)]

  ### RUN RANDOM FOREST ###
  rf_surf <- lapply(1:reps, function(x) {
    permutation_rf(x, for_rf_surf, num_vars, cat_vars)
  })
  rf_bot <- lapply(1:reps, function(x) {
    permutation_rf(x, for_rf_bot, num_vars, cat_vars)
  })

  # Unpack
  surf_r2 <- mean(sapply(rf_surf, "[[", 1))
  partials_df_surf <- bind_rows(sapply(rf_surf, "[", 2))
  ImpData_surf <- bind_rows(sapply(rf_surf, "[", 3)) %>%
    mutate(Layer = "Surface")

  # Unpack
  bot_r2 <- mean(sapply(rf_bot, "[[", 1))
  partials_df_bot <- bind_rows(sapply(rf_bot, "[", 2))
  ImpData_bot <- bind_rows(sapply(rf_bot, "[", 3)) %>%
    mutate(Layer = "Bottom")

  ### PLOT ###
  # Visualize variable importance ----------------------------------------------
  ImpData <- bind_rows(ImpData_surf, ImpData_bot)

  var_names_imp <- sub("\n", " ", var_names)
  imp_levels_used <- rev(recode(vars, !!!var_names_imp))
  names(colors) <- var_names_imp[vars]

  surf_name <- paste0(
    "'Surface ('*italic(n)*' = ", length(unique(for_rf_surf$LakeID)),
    "; mean '*R^2*' = ", round(surf_r2, 2), ")'"
  )
  bot_name <- paste0(
    "'Bottom ('*italic(n)*' = ", length(unique(for_rf_bot$LakeID)),
    "; mean '*R^2*' = ", round(bot_r2, 2), ")'"
  )

  importance <- ImpData %>%
    mutate(
      Var.Names = recode(Var.Names, !!!var_names_imp),
      Var.Names = factor(Var.Names, levels = imp_levels_used),
      Layer = case_match(
        Layer,
        "Surface" ~ surf_name,
        "Bottom" ~ bot_name
      ),
      Layer = factor(Layer, levels = c(surf_name, bot_name))
    ) %>%
    group_by(Var.Names, Layer) %>%
    summarize(
      sd = sd(`%IncMSE`, na.rm = T),
      `%IncMSE` = mean(`%IncMSE`, na.rm = T)
    ) %>%
    ggplot(aes(x = `%IncMSE`, y = Var.Names, fill = Var.Names)) +
    geom_col() +
    geom_errorbar(
      aes(
        xmin = `%IncMSE` - sd,
        xmax = `%IncMSE` + sd
      ),
      width = 0.1
    ) +
    xlab("Variable importance (% increase in MSE)") +
    scale_fill_manual(values = colors) +
    egg::theme_article() +
    theme(
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      panel.grid.minor = element_line(color = "grey95", size = 0.25),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    scale_x_continuous(n.breaks = 4) +
    scale_y_discrete(labels = parse(text = imp_levels_used)) +
    facet_wrap(~Layer, labeller = label_parsed)

  # Set up rug plots
  for_rug_surf <- for_rf_surf %>%
    select(all_of(num_vars)) %>%
    pivot_longer(everything(), names_to = "var", values_to = "x") %>%
    mutate(
      var = recode(var, !!!var_names_part),
      var = factor(var, levels = levels_used)
    )
  for_rug_bot <- for_rf_bot %>%
    select(all_of(num_vars)) %>%
    pivot_longer(everything(), names_to = "var", values_to = "x") %>%
    mutate(
      var = recode(var, !!!var_names_part),
      var = factor(var, levels = levels_used)
    )

  for_rug <- for_rug_surf %>%
    mutate(Layer = "Surface") %>%
    bind_rows(for_rug_bot %>% mutate(Layer = "Bottom")) %>%
    mutate(Layer = factor(Layer, levels = c("Surface", "Bottom"))) %>%
    mutate(x = ifelse(var %in% recode(log_vars, !!!var_names_part),
      exp(x),
      x
    ))

  # Plot partials
  partials_sorted <- partials_df_surf %>%
    mutate(Layer = "Surface") %>%
    bind_rows(partials_df_bot %>% mutate(Layer = "Bottom")) %>%
    mutate(
      var = recode(var, !!!var_names_part),
      var = factor(var, levels = levels_used),
      y = exp(y),
      Layer = factor(Layer, levels = c("Surface", "Bottom")),
      x = ifelse(var %in% recode(log_vars, !!!var_names_part),
        exp(x),
        x
      )
    )

  names(colors) <- var_names_part[vars]
  type <- ifelse(pSat, "% saturation", "concentration (µM)")

  part <- partials_sorted %>%
    group_by(Layer, x, var) %>%
    summarize(
      sd = sd(y, na.rm = T),
      y = mean(y, na.rm = T)
    ) %>%
    ggplot(aes(x = x)) +
    {
      if (pSat) geom_hline(yintercept = 100, lty = "11")
    } +
    geom_ribbon(
      aes(
        ymax = y + sd,
        ymin = y - sd,
        fill = var
      ),
      alpha = 0.5
    ) +
    geom_line(aes(y = y, color = var)) +
    geom_rug(aes(color = var),
      sides = "b", outside = T,
      data = for_rug,
      length = unit(0.05, "inch"), alpha = 0.2
    ) +
    coord_cartesian(clip = "off") +
    ggh4x::facet_grid2(Layer ~ var,
      scales = "free",
      switch = "both",
      labeller = label_parsed
    ) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(x = NULL) +
    egg::theme_article() +
    theme(
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      strip.placement = "outside", # format to look like title
      strip.background.x = element_blank(),
      strip.text.x = element_text(
        margin = unit(c(0, 0, 0, 0), units = "cm"),
        size = 8
      ),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1)
    ) +
    ylab(substitute(
      "Predicted" ~ x ~ type,
      list(
        x = parse(text = gas_name)[[1]],
        type = type
      )
    )) +
    facetted_pos_scales(
      x = list(
        var == "Max~depth~(m)" ~ scale_x_continuous(
          trans = scales::pseudo_log_trans(0.01, 10),
          labels = scales::label_comma(drop0trailing = T),
          breaks = c(0.1, 1, 10, 100, 1000)
        ),
        var == "Mean~depth~(m)" ~ scale_x_continuous(
          trans = scales::pseudo_log_trans(0.01, 10),
          labels = scales::label_comma(drop0trailing = T),
          breaks = c(0.1, 1, 10, 100, 1000)
        ),
        var == "Osgood" ~ scale_x_continuous(
          trans = scales::pseudo_log_trans(0.001, 10),
          labels = scales::label_comma(drop0trailing = T),
          breaks = c(0.1, 1, 10, 100, 1000)
        ),
        var == "TP~(µg~L^-1)" ~ scale_x_continuous(
          trans = scales::pseudo_log_trans(1, 10),
          labels = scales::label_comma(drop0trailing = T),
          breaks = c(1, 10, 100, 1000)
        ),
        var == "Buoy~freq~(1/s)" ~ scale_x_continuous(
          trans = scales::pseudo_log_trans(0.0001, 10),
          labels = scales::label_comma(drop0trailing = T),
          breaks = c(0.001, 0.01, 0.1)
        ),
        var == "SA~(km^2)" ~ scale_x_continuous(
          trans = scales::pseudo_log_trans(0.0001, 10),
          labels = scales::label_comma(drop0trailing = T),
          breaks = c(0.001, 0.1, 10, 1000)
        ),
        !var %in% recode(log_vars, !!!var_names_part) ~ scale_x_continuous(labels = scales::label_comma())
      )
    )

  # Combine
  plot_gas <- ggpubr::ggarrange(
    importance,
    part,
    ncol = 1, heights = c(1, 1.7),
    labels = c("a", "b"),
    label.x = 0.03
  )
  return(plot_gas)
}
