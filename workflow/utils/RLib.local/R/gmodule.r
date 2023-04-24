###
#' @Date: 2022-07-02 22:51:37
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-08-27 16:34:24
#' @FilePath: /2021_09-MT10kSW/workflow/utils/RLib.local/R/gmodule.r
#' @Description:
###


load_gmodule <- function(gmodule_path, symnum = TRUE) {
  gmodule=
    gmodule_path %>%
    {
      df = read.csv(., row.names = 1)
      colnames(df) =
        read.csv(., header = FALSE)[1, -1] %>%
        gsub("^.+\\.\\.(.+\\.\\..+).bam", "\\1", .)
      df
    } %>%
    as.matrix
  if (symnum) {
    return(symnum.gmodule(gmodule))
  } else {
    return(gmodule)
  }
}
symnum.gmodule <- function(gmodule) {
  gmodule %>%
    symnum(cutpoints = c(0, 0.2, 0.75, 1), symbols = c(0, 0.5, 1), legend = FALSE)
}


gmodule_heatmap_color = c(rep("#4DBBD5", 4), rep("#FFDC91", 11), rep("#E64B35", 5))
gmodule_heatmap_legend_breaks = c(0, 0.2, 0.5, 0.75, 1)


get_module.name <- function() {
  module.name = read.csv("data/module_name.tsv",
                         sep = "\t", header = 1) %>%
    {row.names(.) = .$entry; .}
}


get_taxon_group <- function(genome.relative_abundance, row.name = "Genome", col.name = "Group") {
  genome.relative_abundance %>%
    reshape2::dcast(formula(paste(row.name, "~", col.name)),
                    value.var = col.name, fill = "", fun.aggregate = unique)
}
