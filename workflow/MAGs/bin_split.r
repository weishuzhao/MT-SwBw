###
#* @Date: 2022-03-30 00:30:15
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-05-21 12:55:58
#' @FilePath: /2022_05-ZFMG-release/workflow/MAGs/bin_split.r
#* @Description:
#    Those both in MT and ME | only in MT
#    at class level (Fig.3)
# TODO: Siegel test
###
source("Scripts/utils/RLib.r")

sample_meta_col = c(
  "ME" = "#d7191c", "MT" = "#2c7bb6"
)

### ######################################################################## ###
#### LOADING MAGs INFORMATIONS                                              ####
### ######################################################################## ###
Stdb =
  merge(read.csv("workflow/MAGs/Stdb.csv"),
        read.csv("workflow/MAGs/Wtdb.csv")[, c("classification",
                                               "cluster")],
        by.x = "secondary_cluster", by.y = "cluster") %>%
  {.$location = ifelse(grepl("^1", .$sample), "ME", "MT"); .} %>%
  {.$size = .$Genome.size..bp.; .$Genome.size..bp. <- NULL; .}


MAG_info =
  read.csv("workflow/MAGs/growthsnake.csv") %>%
  {.$genome = .$MAG; .$MAG <- NULL; .} %>%
  merge(Stdb, by = "genome")
names(MAG_info)


# FUNCTION FOR SIGNIFICANT TEST >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
library(ggExtra)
library(multcomp)
library(tidyr)
library(patchwork)


# function: test difference between `group` on `value`
# data: long format, include 2 cols (value, group)
# value: numeric value
# group: to group the value
# return: data.frame(
#   index = group,
#   columns = c(
#     char,  # marker
#     locat,  # max locat of the value
#     locat.mark,  # place marker shall locat
#     group
#   )
# )
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
  signif.mark = data.frame(char = test.b$mcletters$Letters,
                           locat = unlist(lapply(split(data[,value],
                                                       data[,group]),
                                                 max)))
  signif.mark["locat.mark"] =
    signif.mark$locat %>%
    {. + . * max.locat.scale}
  signif.mark[group] = rownames(signif.mark)
  signif.mark
}


ggpoint_signif <- function(
  p, igroup,
  x.geom = geom_boxplot, y.geom = geom_boxplot,
  static.mark = "*",
  p.width = 3, p.height = 3,
  p_theme = NULL
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
      geom_signif(comparisons = comparisons, map_signif_level = T,
                  step_increase = 0.1)
    py = py +
      geom_signif(comparisons = comparisons, map_signif_level = T,
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
  {p +
    scale_x_continuous(limits = xlimits) +
    scale_y_continuous(limits = ylimits)} %>%
    {if (is.null(p_theme)) . else p_theme(.)} +
    subadjust(px) + scale_x_continuous(limits = xlimits) +
    subadjust(py) + scale_y_continuous(limits = ylimits) +
    plot_layout(
      design = c(patchwork::area(t = 2, l = 1, b = 1 + p.height, r = p.width),
                 patchwork::area(t = 1, l = 1, b = 1, r = p.width),
                 patchwork::area(t = 2, l = 1 + p.width,
                                 b = 1 + p.height, r = 1 + p.width)),
      guides = 'collect')
  # use patchwork to organize picture <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

}
# <<-<<-<<                                                            <<-<<-<< #


# plot >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
p_theme <- function(p) {
  font_size_1 = 18
  font_size_2 = 14
  axis.ticks.length = 0.3

  p +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(color = 'black', fill = 'transparent'),
      axis.title = element_text(
        size = font_size_1, face = "bold", colour = "black"),
      axis.text = element_text(colour = "black",
                               face = "bold", size = font_size_2),
      legend.title = element_text(size = font_size_2, face = "bold"),
      legend.text = element_text(size = font_size_2, face = "bold"),
      axis.ticks.length = unit(axis.ticks.length, 'cm')
      #legend.position = "bottom"
    ) +
    theme(text = element_text(family = "Arial",
                              size = 16,
                              hjust = 0.5,
                              lineheight = 0.5)) +
    theme(plot.title = element_text(hjust = 0.5))
}


### ######################################################################## ###
#### report information of MAGs                                             ####
### ######################################################################## ###
MAG_info %>% {split(.$size / 1e6, .$location)} %>% {lapply(., mean)} %>% bind_rows(.id = "location")
MAG_info %>% {split(.$GC, .$location)} %>% {lapply(., mean)} %>% bind_rows(.id = "location")
MAG_info %>% {split(.$OGT, .$location)} %>% {lapply(., mean)} %>% bind_rows(.id = "location")
MAG_info %>% {split(.$MGT, .$location)} %>% {lapply(., mean)} %>% bind_rows(.id = "location")
p = ggplot(MAG_info,
           aes_string(x = "size / 1e6",
                      y = "GC",
                      color = "location",
                      fill = "location")) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = sample_meta_col) +
  scale_fill_manual(values = sample_meta_col) +
  labs(x = "Genome Size (Mbp)",
       y = "G+C contant (%)")

psignif = ggpoint_signif(
  p, igroup = list("location" = sample_meta_col),
  x.geom = geom_boxplot, y.geom = geom_boxplot,
  static.mark = "*", p.width = 5, p.height = 5,
  p_theme = p_theme)
ggsave(filename = "figs/3B.svg",
       plot = psignif, dpi = 300, width = 7.5, height = 6)

# ============================================================================ #
p = ggplot(MAG_info,
           aes_string(x = "OGT",
                      y = "MGT",
                      color = "location",
                      fill = "location")) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = sample_meta_col) +
  scale_fill_manual(values = sample_meta_col) +
  labs(x = "Optimal growth temperature (℃)",
       y = "Minimum generation time (h)")

psignif = ggpoint_signif(
  p, igroup = list("location" = sample_meta_col),
  x.geom = geom_boxplot, y.geom = geom_boxplot,
  static.mark = "*", p.width = 5, p.height = 5,
  p_theme = p_theme)
ggsave(filename = "figs/3C.svg",
       plot = psignif, dpi = 300, width = 7.5, height = 6)


colnames(MAG_info)
### ######################################################################## ###
#### report information of MAGs grouped by Venn.class                       ####
### ######################################################################## ###
# group MAGs according to Venn at class >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
collapse.div <- function(div.raw, collapse.name) {
  div.collapse = aggregate(
    div.raw, by = list(collapse.name = collapse.name), FUN = sum)
  rownames(div.collapse) = div.collapse$collapse.name
  div.collapse$collapse.name <- NULL
  return(div.collapse)
}

venn.class =
  reshape2::acast(
    MAG_info, formula = formula("classification ~ location"),
    fun.aggregate = length) %>%
  collapse.div(div.raw = ., taxon.split(rownames(x = .), 1, "class")) %>%
  {. > 0} %>%
  as.data.frame(.) %>%
  {function(x)
    ifelse((x["ME"]) & (!x["MT"]), "ME.only", "") %>%
      ifelse((!x["ME"]) & (x["MT"]), "MT.only", .) %>%
      ifelse((x["ME"]) & (x["MT"]), "both", .)
  }(.)
MAG_info$venn.class =
  MAG_info$classification %>%
  taxon.split(., 1, "class") %>%
  venn.class[.] %>%
  {ifelse(. == "both", str_glue("{MAG_info$location}.both"), .)} %>%
  factor(x = ., levels = c("ME.only", "ME.both", "MT.both", "MT.only"))

venn.class.color = c(
  "ME.only" = "#d7191c",
  "ME.both" = "#f28d8e",
  "MT.both" = "#98c5e6",
  "MT.only" = "#2c7bb6"
)
# <<-<<-<<                                                            <<-<<-<< #

##### PLOT                                                                 #####
p = ggplot(MAG_info,
           aes_string(x = "size",
                      y = "GC",
                      color = "venn.class",
                      fill = "venn.class",
                      shape = "venn.class")) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = venn.class.color) +
  scale_fill_manual(values = venn.class.color) +
  scale_shape_manual(values = c("ME.only" = 24,
                                "ME.both" = 21,
                                "MT.both" = 21,
                                "MT.only" = 25)) +
  labs(x = "Genome Size",
       y = "GC")

psignif = ggpoint_signif(
  p, igroup = list("venn.class" = venn.class.color),
  x.geom = geom_boxplot, y.geom = geom_boxplot,
  static.mark = 0.05, p.width = 5, p.height = 5,
  p_theme = function(p) p_theme(
    p +
      scale_x_continuous(labels =
                           ~format(.x, scientific = TRUE) %>%
                           str_replace("^0e\\+0", "0e+0") %>%
                           str_replace("e\\+0", "%*%10^") %>%
                           parse(text = .))
  ))
ggsave(filename = "figs/S8A.svg",
       plot = psignif, dpi = 300, width = 7.5, height = 6)

# ============================================================================ #
p = ggplot(MAG_info,
           aes_string(x = "OGT",
                      y = "MGT",
                      color = "venn.class",
                      fill = "venn.class",
                      shape = "venn.class")) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = venn.class.color) +
  scale_fill_manual(values = venn.class.color) +
  scale_shape_manual(values = c("ME.only" = 24,
                                "ME.both" = 21,
                                "MT.both" = 21,
                                "MT.only" = 25)) +
  labs(x = 'optimal growth temperature',
       y = 'minimum generation time (h)')

psignif = ggpoint_signif(
  p, igroup = list("venn.class" = venn.class.color),
  x.geom = geom_boxplot, y.geom = geom_boxplot,
  static.mark = 0.05, p.width = 5, p.height = 5,
  p_theme = p_theme)
ggsave(filename = "figs/S8B.svg",
       plot = psignif, dpi = 300, width = 7.5, height = 6)
# <<-<<-<<                                                            <<-<<-<< #
