###
#' @Date: 2022-06-28 15:37:09
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2023-04-09 17:24:58
#' @FilePath: /2022_09-M_mem/workflow/utils/libs/metaSC/R/RLib/R/total_percent_bar.r
#' @Description:
###

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))


#' @title get total abundance of div.otu
#'
#' @param div.raw a long table of diversity, like:
#'
#'                colnames: <sample.name>, <value.name>, <*keey.cols>
#'                values  : <any>        , <numeric>   , <*any>
#'                example : Sample       , Reads       , Group, Layer, ...
#'
#'                warning: these clumns name will be overwritten:
#'                    ..Sample, ..All, Total
#'
#' @return a smaller row of total diversity, like:
#'
#'                colnames: (<index>), <sample.name>, <value.name>, <*keey.cols>
#'                values  : (<any>  ), <any>        , <numeric>   , <*any>
#'                example : (Sample ), Sample       , Total       , Group, Layer
#'
get_total_ce <-
  function(div.raw, value.name,
           sample.name = "Sample",
           keep.cols = c("Group", "Layer")) {
    keep.div.otu <-
      div.raw %>%
      {
        .$..Sample <- .[, sample.name]
        .
      } %>%
      .[c("..Sample", sample.name, keep.cols)] %>%
      unique()

    div.raw %>%
      group_by(..Sample = .[, sample.name]) %>%
      dplyr::summarise(Total = sum(get(value.name))) %>%
      merge(keep.div.otu) %>%
      data.frame(row.names = .$..Sample) %>%
      {
        .$..Sample <- NULL
        .
      }
  }


get_total_plot <- function(
    div.raw, value.name,
    sample.name = "Sample", fill = "Group",
    labs.x = "", labs.y = "",
    font_size_1 = 15, font_size_2 = 13, font_size_3 = 10,
    axis.ticks.length = 0.1) {
  total_ce <- get_total_ce(div.raw, value.name)

  p_total_mapped_rate <-
    ggplot(
      data = total_ce,
      mapping = aes_string(x = sample.name, y = "Total", fill = fill)
    ) +
    geom_bar(
      stat = "identity",
      col = "black", size = 0.3
    ) +
    scale_x_discrete(limits = rownames(total_ce), labels = NULL) +
    labs(title = "", x = labs.x, y = labs.y) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.text = element_text(
        size = font_size_2, colour = "black", face = "bold"
      ),
      axis.title = element_text(
        size = font_size_1, face = "bold", colour = "black"
      ),
      axis.ticks.length = unit(0, "cm"),
      legend.title = element_text(size = font_size_2, face = "bold"),
      legend.key = element_blank()
    )

  return(p_total_mapped_rate)
}


get_relative_ce <- function(
    div.raw, value.name, sample.name = "Sample",
    total_ce = NULL) {
  if (is.null(total_ce)) {
    total_ce <- get_total_ce(div.raw, value.name, sample.name = sample.name)
  }
  div.raw %>%
    as.data.frame() %>%
    {
      .[, value.name] / total_ce[as.character(.[, sample.name]), "Total"] * 100
    }
}

#' @title draw a picture of relative abundance
#'
#' @param div.raw a long table of diversity, like:
#'
#'                colnames: <sample.name>, <value.name>, <*keey.cols>
#'                values  : <any>        , <numeric>   , <*any>
#'                example : Sample       , Reads       , Group, Layer, ...
#'
#' @param value.name abundance column name
#' @param fill.name name to group by, and to be colored
#' @param sample.name name of sample to define row names
#' @return plot of relative abundance in each environment
get_percent_plot <- function(
    div.raw, value.name, fill.name, sample.name = "Sample",
    labs.x = "", labs.y = "",
    TOP_N_TAXON_PER_LAYER = NULL, MIN_REPORT_TAXON_PERCENT = 0,
    font_size_1 = 15, font_size_2 = 13, font_size_3 = 10,
    axis.ticks.length = 0.1,
    geom_bar.color = "#454545") {
  total_ce <- get_total_ce(div.raw, value.name, sample.name)

  ce <-
    div.raw %>%
    as.data.frame() %>%
    mutate(
      Abundance = get_relative_ce(., value.name, sample.name, total_ce)
    ) %>%
    {
      .$name <- .[, fill.name]
      .
    } %>%
    {
      if (!is.null(TOP_N_TAXON_PER_LAYER)) {
        top_taxon <-
          reshape2::acast(
            ., formula(paste0("name ~ ", sample.name)),
            value.var = "Abundance", fill = 0, fun.aggregate = sum
          ) %>%
          otu.div.top(TOP_N_TAXON_PER_LAYER, min.value = 0) %>%
          sort()
      } else if (MIN_REPORT_TAXON_PERCENT > 0) {
        top_taxon <-
          split(.[c("Abundance", "name")], .[sample.name]) %>%
          lapply(
            . %>%
              group_by(name = get("name")) %>%
              summarise(Abundance = sum(get("Abundance")))
          ) %>%
          bind_rows(.id = sample.name) %>%
          {
            .$name[.$Abundance >= MIN_REPORT_TAXON_PERCENT]
          } %>%
          unique() %>%
          as.character() %>%
          sort()
      } else {
        top_taxon <- NULL
      }
      if (!is.null(top_taxon)) {
        .$name <- .$name %>%
          as.character() %>%
          {
            ifelse(. %in% top_taxon, ., "others")
          } %>%
          factor(levels = c(top_taxon, "others"))
      }
      .
    } %>%
    {
      .
    }

  p_relative_abundance <-
    ggplot(
      data = ce,
      mapping = aes_string(
        x = sample.name, y = "Abundance",
        fill = "name"
      )
    ) +
    geom_bar(
      position = position_stack(reverse = TRUE),
      stat = "identity",
      color = geom_bar.color, size = 0.3
    ) +
    paletteer::scale_fill_paletteer_d("ggsci::default_igv") +
    scale_x_discrete(limits = rownames(total_ce)) +
    # guides(fill = "none") +
    guides(fill = guide_legend(title = fill.name, ncol = 1, reverse = TRUE)) +
    # theme(axis.text.x = element_text(colour = axis.x.color)) +
    labs(title = "", x = labs.x, y = labs.y)

  p1 <-
    p_relative_abundance +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
    ) +
    theme(
      axis.line = element_line(colour = "black"),
      axis.text = element_text(
        size = font_size_2, colour = "black", face = "bold"
      ),
      axis.title = element_text(
        size = font_size_1, face = "bold", colour = "black"
      ),
      axis.ticks.length = unit(axis.ticks.length, "cm"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    theme(
      legend.title = element_text(size = font_size_2, face = "bold"),
      legend.text = element_text(size = font_size_3),
      legend.key = element_blank()
    ) +
    theme(text = element_text(
      family = "Arial",
      size = font_size_1,
      hjust = 0.5,
      lineheight = 0.5
    )) +
    theme(plot.title = element_text(hjust = 0.5))

  return(p1)
}
