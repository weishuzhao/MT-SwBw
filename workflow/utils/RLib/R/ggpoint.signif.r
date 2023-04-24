###
#* @Date: 2022-04-25 00:51:52
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2022-12-21 14:13:31
#' @FilePath: /metaSC/R/RLib/R/ggpoint.signif.r
#* @Description:
###
# FUNCTION FOR SIGNIFICANT TEST >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
suppressMessages(library(ggExtra))
suppressMessages(library(multcomp))
suppressMessages(library(tidyr))
suppressMessages(library(patchwork))


#' @title test difference between `group` on `value`
#'
#' @param data: long format, include 2 cols (value, group)
#' @param value: one of colname in data, numeric value to compare
#' @param group: one of colname in data, group the value
#' @return: data.frame(
#'
#'            index = group,
#'
#'            columns = c(
#'
#'              char,  # marker
#'
#'              locat,  # max locat of the value
#'
#'              locat.mark,  # place marker shall locat
#'
#'              group
#'
#'            )
#'
#'          )
#'
#' @export
group.signif.mark <- function(
  data, value, group, signif.alpha = 0.05,
  max.locat.scale = 0.1, glht.mcp.test = "Tukey"
) {
  test.b = cld(
    glht(aov(formula(paste0(value, " ~ .glht_mcp_factor")),
             cbind(data,
                   .glht_mcp_factor = as.factor(data[,group]))),
         linfct = mcp(.glht_mcp_factor = glht.mcp.test)),
    alpah = signif.alpha)
  stat =
    data %>%
    {split(.[, value], .[, group])} %>%
    lapply(. %>% {c("min" = min(.), "mean" = mean(.), "max" = max(.))}) %>%
    bind_rows(.id = "group")
  signif.mark = data.frame(char = test.b$mcletters$Letters, stat)
  signif.mark["locat.mark"] =
    signif.mark$max %>%
    {. + . * max.locat.scale}
  signif.mark[group] = rownames(signif.mark)
  signif.mark
}

#' @title add signif figure on a given figure
#'
#' @param p a ggplot geom_point plot
#'          g must define `mapping` in `ggplot` but not `geom_point`.
#'
#' @param igroup group and color information
#'               such as `list("env" = sample_meta_col)`
#' @param x.geom,y.geom geom plot type
#' @param static.mark if NULL, show nothing
#'
#'                    if "*", show paired wilcox test between any pairs
#'                            in the group
#'
#'                    if numeric, show character difference given by
#'                                turkey test with adjust P-value lower
#'                                than that level.
#' @return a patchwork picture
#' @export
ggpoint.signif <- function(
  p, igroup,
  x.geom = geom_boxplot, y.geom = geom_boxplot,
  static.mark = "*",
  p.width = 3, p.height = 3,
  p_theme = NULL,
  scale_x = scale_x_continuous,
  scale_y = scale_y_continuous
) {
  # extract data from plot OBJECT >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
  pbuild = ggplot_build(p)

  pdata = pbuild$plot$data
  xname = quo_name(pbuild$plot$mapping$x)
  yname = quo_name(pbuild$plot$mapping$y)
  # extract data from plot OBJECT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

  # get groups AND comparisons >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
  fillname = names(igroup)
  comparisons = combn(names(igroup[[fillname]]), 2, FUN = list)
  # get groups AND comparisons <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

  # draw plot at side >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
  px = ggplot(data = pdata,
              mapping = aes_string(
                y = fillname, fill = fillname, #color = fillname,
                x = xname)) +
    x.geom() +
    scale_fill_manual(values = igroup[[fillname]])
  py = ggplot(data = pdata,
              mapping = aes_string(
                x = fillname, fill = fillname, #color = fillname,
                y = yname)) +
    y.geom() +
    scale_fill_manual(values = igroup[[fillname]])
  # add static marks such as '*' or 'a'/'b' ================================== #
  if (is.null(static.mark)) {
  } else if (static.mark == "*") {
    px = px +
      ggsignif::geom_signif(comparisons = comparisons, map_signif_level = T,
                            step_increase = 0.1)
    py = py +
      ggsignif::geom_signif(comparisons = comparisons, map_signif_level = T,
                            step_increase = 0.1)
  } else if (is.double(static.mark) && 0 < static.mark && static.mark < 1) {
    px = px +
      geom_text(data = group.signif.mark(pdata, xname, fillname,
                                         signif.alpha = static.mark),
                aes_string(x = "locat.mark", y = fillname, label = "char"))
    py = py +
      geom_text(data = group.signif.mark(pdata, yname, fillname,
                                         signif.alpha = static.mark),
                aes_string(y = "locat.mark", x = fillname, label = "char"))
  } else {
    warning("unknown static marker: ", static.mark)
  }
  # draw plot at side <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

  # reset limits >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
  xlimits = c(pbuild$layout$panel_scales_x[[1]]$get_limits(),
              ggplot_build(px)$layout$panel_scales_x[[1]]$get_limits())
  xlimits = c(min(xlimits), max(xlimits))
  ylimits = c(pbuild$layout$panel_scales_y[[1]]$get_limits(),
              ggplot_build(py)$layout$panel_scales_y[[1]]$get_limits())
  ylimits = c(min(ylimits), max(ylimits))
  # reset limits <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

  subadjust = function(x) {
    x + labs(x = "", y = "") +
      guides(fill = "none") +
      theme_void() +
      theme(axis.text = element_blank())
  }

  # use patchwork to organize picture >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
  {p + scale_x(limits = xlimits) + scale_y(limits = ylimits)} %>%
    {if (is.null(p_theme)) . else p_theme(.)} +
    subadjust(px) + scale_x(limits = xlimits) +
    subadjust(py) + scale_y(limits = ylimits) +
    plot_layout(
      design = c(patchwork::area(t = 2, l = 1, b = 1 + p.height, r = p.width),
                 patchwork::area(t = 1, l = 1, b = 1, r = p.width),
                 patchwork::area(t = 2, l = 1 + p.width,
                                 b = 1 + p.height, r = 1 + p.width)),
      guides = 'collect')
  # use patchwork to organize picture <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
}
