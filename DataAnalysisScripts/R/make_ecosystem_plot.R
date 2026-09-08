make_ecosystem_plot <- function(comb, gas, text) {
  
  ggplot_data <- comb %>%
    filter(name == gas)
  
  ggplot(ggplot_data, aes(y = Ecosystem, x = value)) +
    # Without quantiles for shading
    geom_violin(aes(color = Ecosystem, fill = Ecosystem), alpha = 0.1) +
    # 25th and 75th percentiles
    geom_violin(
      aes(color = Ecosystem),
      fill = NA,
      quantiles = c(0.25, 0.75),
      quantile.linetype = "11"
    ) +
    # Median
    geom_violin(
      aes(color = Ecosystem),
      fill = NA,
      quantiles = 0.5,
      quantile.linetype = "solid"
    ) +
    geom_text(
      data = text %>% filter(name == gas),
      aes(x = xpos, label = label, color = Ecosystem,
          hjust = hjust),
      vjust = -0.5,
      size = 2.5
    ) +
    scale_x_log10() +
    xlab(bquote(.(parse(text = gas)[[1]]) ~ "concentration (µM)")) +
    scale_color_manual(values = c(
      rep("grey70", length(unique(ggplot_data$Ecosystem)) - 2),
      "turquoise4", "turquoise4"
    )) +
    scale_fill_manual(values = c(
      rep("grey70", length(unique(ggplot_data$Ecosystem)) - 2),
      "turquoise4", "turquoise4"
    )) +
    egg::theme_article() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "none",
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor.x = element_line(color = "grey95", linewidth = 0.25)
    )
}

make_ecosystem_plot_wide <- function(comb, text){
  letters <- comb %>%
    group_by(name) %>%
    group_modify(~{
      
      # Dunn's test
      dunn <- rstatix::dunn_test(
        .x,
        value ~ Ecosystem
      )
      
      # Named vector of adjusted p-values
      pvals <- dunn$p.adj
      names(pvals) <- paste(dunn$group1, dunn$group2, sep = "-")
      
      # Compact letter display
      cld <- multcompView::multcompLetters(pvals)
      
      tibble(
        Ecosystem = names(cld$Letters),
        letter = cld$Letters,
      )
    }) %>%
    left_join(
      comb %>%
        filter(value >= 0) %>%
        group_by(name, Ecosystem) %>%
        summarize(
          ymax = max(value, na.rm = TRUE),
          .groups = "drop"
        ),
      by = c("name", "Ecosystem")
    ) %>%
    mutate(Ecosystem = factor(Ecosystem, levels = factor_levels))
  
  comb %>% 
    filter(value >= 0) %>% 
    mutate(value = if_else(value == 0, 
                           min(value[value > 0], na.rm = T), 
                           value, missing = NA)) %>% 
    ggplot(aes(x = Ecosystem, y = value)) + 
    # Without quantiles for shading 
    geom_violin(aes(color = Ecosystem, fill = Ecosystem), alpha = 0.1) + 
    # With 25 and 75% 
    geom_violin(aes(color = Ecosystem), fill = NA, quantiles = c(0.25, 0.75), 
                quantile.linetype = "11") + 
    # With median 
    geom_violin(aes(color = Ecosystem), fill = NA, quantiles = c(0.5), 
                quantile.linetype = "solid") + 
    geom_text(aes(y = ypos, label = label), 
              data = text, color = "grey20", fontface = "italic",
              hjust = 0.5, vjust = 0.5, size = 2, lineheight = 0.8) +
    geom_text(aes(y = ypos*7, label = label_median), 
              data = text, color = "grey20", fontface = "bold",
              hjust = 0.5, vjust = 0.5, size = 2, lineheight = 0.7) +
    geom_text(
      data = letters,
      vjust = -0.5,
      size = 3,
      aes(
        y = ymax, 
        label = letter,
        color = Ecosystem
      )) + 
    scale_y_log10(labels = scales::label_comma(drop0trailing = T,
                                               accuracy = 0.001),
                  minor_breaks = minor_log_breaks)+ 
    ylab("Concentration (µM)")+ 
    scale_color_manual(values = c(rep("grey70", n_distinct(comb$Ecosystem) - 2), 
                                  "turquoise4", "turquoise4") ) + 
    scale_fill_manual(values = c(rep("grey70", n_distinct(comb$Ecosystem) - 2), 
                                 "turquoise4", "turquoise4") ) + 
    egg::theme_article()+ 
    facet_wrap(~name, scales = "free", space = "free_x",
               labeller = "label_parsed")+
    theme(axis.title.x = element_blank(), 
          legend.position = "none", 
          axis.text.x = ggtext::element_markdown(
            angle = 45,
            hjust = 1,
            vjust = 1
          ),
          panel.grid.major.y = element_line(color = "grey90", size = 0.5), 
          panel.grid.minor.y = element_line(color = "grey95", size = 0.25))
}

minor_log_breaks <- function(lims) {
  exponents <- seq(
    floor(log10(min(lims))),
    ceiling(log10(max(lims)))
  )
  
  breaks <- outer(1, 10^exponents)
  
  breaks <- as.vector(breaks)
  breaks[breaks > min(lims) & breaks < max(lims)]
}
