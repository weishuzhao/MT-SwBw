###
#' @Date: 2022-09-01 14:55:02
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-09-01 15:38:44
#' @FilePath: /2021_09-MT10kSW/workflow/utils/RLib.local/R/abundance_group_bar.r
#' @Description:
###

#' @title compare significance of abundance between 2*2 groups
#'
#' @param relative_abundance.class data.frame of at least 5 columns.
#' @param TOTAL_NUM_ADJUST total number of value 'x'
#' @param df_x <data.frame> describe how each bar grouped by bar in figure (x and the right var in grid_formula)
#' @param df_fill <data.frame> describe how to split 'group' to compare (fill and the left var in grid_formula)
#' @param x rowname of each bar in figure
#' @param y <float> abundance to stat
#' @param fill color of bar in each x in figure
#' @param grid_formula formula used in 'facet_grid', such as `formula("Location ~ Phylum")`
#' @param p.value.char.filter decide what significant label will show in figuer
#'
#' @return figure
report__abundance_bar_group <- function(
  relative_abundance.class,
  TOTAL_NUM_ADJUST = NA, df_x = NA, df_fill = NA,
  x = "Taxa_label", y = "Relative_abundance", fill = "Group",
  grid_formula = formula("Location ~ Phylum"),
  p.value.char.filter = . %>% .[!grepl("[-?]", .$p.value.char),]
) {
  if (is.na(TOTAL_NUM_ADJUST)) {
    TOTAL_NUM_ADJUST =
      relative_abundance.class$Taxa_label %>% unique %>% length %>% {. * 2}
  }
  if (is.na(df_x)) {
    df_x =
      relative_abundance.class %>%
      .[c(all.vars(grid_formula)[2], x)] %>%
      unique
  }
  if (is.na(df_fill)) {
    df_fill =
      relative_abundance.class %>%
      .[c(fill, all.vars(grid_formula)[1])] %>%
      unique
  }

  relative_abundance.class.signif =
    relative_abundance.class %>%
    location.group.signif(x, y, all.vars(grid_formula)[1], fill, TOTAL_NUM_ADJUST) %>%
    merge(df_x) %>%
    p.value.char.filter

  relative_abundance.class =
    relative_abundance.class %>%
    {merge(merge(df_x, df_fill), ., all.x = TRUE)} %>%
    {.[is.na(.[y]), "Hide"] = "TRUE"; .[is.na(.)] = 1; .}

  #### add paint                                                            ####
  font_size_1 = 16
  font_size_2 = 14
  axis.ticks.length = 0.1

  p =
    ggplot(data = relative_abundance.class) +
    geom_boxplot(mapping = aes_string(x = x, y = y, fill = fill,
                                      linetype = "Hide")) +
    scale_linetype_manual(values = c("1" = 1, "TRUE" = 0)) +
    guides(linetype = "none") +
    scale_fill_manual(values = sample_meta_col) +

    facet_grid(grid_formula,
               scales = "free_x", space = "free_x") +

    geom_label(data =
                 relative_abundance.class.signif %>%
                 .[.$p.value.char != "", ],
               mapping = aes_string(x = x, y = "Max. + max(Max.) * 0.05",
                                    label = "p.value.char"),
               label.padding = unit(0.05, "lines"), label.size = 0,
               fill = "#7f7f7f3f") +
    labs(x = "")

  p1 =
    p +
    scale_y_log10(labels =
                    ~format(.x, scientific = TRUE) %>%
                    str_replace("^0e\\+0", "0e+0") %>%
                    str_replace("e\\+0", "%*%10^") %>%
                    parse(text = .)) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = 'white', size = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    ) +
    theme(
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = font_size_2, colour = "black", face = "bold"),
      axis.title = element_text(size = font_size_1, face = "bold", colour = "black"),
      axis.ticks.length = unit(axis.ticks.length, 'cm'),

      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "plain")
    ) +
    theme(
      legend.title = element_text(size = font_size_2, face = "bold"),
      legend.text = element_text(size = font_size_2, face = "bold")
      #legend.position = "bottom"
    ) +
    theme(text = element_text(size = font_size_1,
                              hjust = 0.5,
                              lineheight = 0.5)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.key = element_blank()) +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.x = element_blank())

  p1
}
