###
#' @Date: 2022-03-03 15:56:58
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-08-05 22:17:52
#' @FilePath: /2021_09-MT10kSW/workflow/MAGs_tpm/relative_abundance_heatmap.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


Wtdb_abd = argv[1]
#Wtdb_abd = stringr::str_glue("Stdb.relative_abundance.tsv") %>% file_path$file_path$results() %>% as.character
fig_out = argv[2]
existance_or_abundance = argv[3]
taxon.level.spec = argv[4]


### ######################################################################## ###
#### extract information of MAGs                                            ####
### ######################################################################## ###
genome_taxonomy = load__genome_taxonomy(load__Stdb(), load__Wtdb())
genome.relative_abundance = get_relative_abundance(Wtdb_abd, genome_taxonomy)


taxon.common.envs <- function(genome.relative_abundance, envs) {
  genome.relative_abundance %>%
    #{message.print(head(.), dim(.))}
    .[c("name", "Group")] %>%
    unique %>%
    .[.[, "Group"] %in% envs, ] %>%
    .$name %>%
    .[duplicated(.)] %>%
    unique
}


### ######################################################################## ###
#### DRAW pheatmap and output                                               ####
### ######################################################################## ###
##### MAKE OTU TANLE OF taxonomy ~ sample                                  #####
div.otu =
  genome.relative_abundance %>%
  {.$name = taxon.split(.$Taxonomy, 1, taxon.level.spec); .} %>%
  {.[.$name %in% c(taxon.common.envs(., c("Ss", "Sw", "Bs", "Bw"))),]} %>%
  {.$Abundance = get_relative_ce(., "Relative_abundance"); .} %>%
  reshape2::acast(formula("name ~ Sample"),
                  value.var = "Abundance", fill = 0, fun.aggregate = sum)

##### MAKE annotation for cols and rows                                    #####
annotation_col =
  sample_meta %>%
  {
    .$Layer =
      .$Layer %>%
      #gsub("water", "", .) %>%
      gsub("^(.+)\\.(.+)$", "\\2", .) %>%
      gsub("^(.+)-(.+)$", "\\1", .) %>%
      gsub("^0$", "surf", .) %>%
      gsub("^\\d$", "mid", .) %>%
      gsub("^\\d\\d+$", "deep", .) %>%
      factor(levels = c("surf", "mid", "deep"))
    .
  } %>%
  {
    data.frame(.[c("Depth", "Group", "Layer")],
               row.names = str_glue("{.$Group}_{.$Site}..{.$Layer}"))
  #} %>% {.[.$Layer != -1, ]} %>% reshape2::acast("Site ~ Layer") %>% t %>% pheatmap::pheatmap(cluster_rows = FALSE)
  }
annotation_row =
  div.otu %>%
  reshape2::melt(varnames = c("name", "Sample"), value.name = "Abundance") %>%
  .[.$Abundance > 0,] %>%
  {.$Group = annotation_col[as.character(.$Sample), "Group"]; .} %>%
  get_taxon_group("name", "Group")
annotation_colors = append(
  sample_meta_col %>%
    Map(function(x, i) {c(i, "#FFFFFF") %>% {names(.) = c(x, ""); .}}, names(.), .),
  list(
    Group = sample_meta_col,
    Layer =
      paletteer::paletteer_c("ggthemes::Brown", 3) %>%
      c %>%
      {names(.) = c("surf", "mid", "deep"); .}
    )
)


##### PLOT pheatmap OF abundance markered BY site                          #####
plot_height = 3 + round(0.134 * dim(div.otu)[1])
plot_width = 18

svg(filename = fig_out, width = plot_width, height = plot_height)
if (existance_or_abundance == "existance") {
  p =
    div.otu %>%
    {
      d = .

      p = pheatmap::pheatmap(
        (. > 0) + 0,
        cellwidth = 10, cellheight = 10,
        annotation_col = annotation_col,
        annotation_row = annotation_row,
        annotation_colors = annotation_colors
      )
      p
    }
} else if (existance_or_abundance == "Abundance") {
  p =
    div.otu %>%
    {
      legend_labels = c(0, 0.2, 1, 2, 5, 10, 20, 30, 50)
      d = .

      p = pheatmap::pheatmap(
        {.[. > 50] = 50; .} %>% {log10(. + 1)} %>% {.[. == 0] = NA; .},
        cellwidth = 10, cellheight = 10,
        cluster_rows = hclust(dist(d, method = "euclidean"), "complete"),
        cluster_cols = hclust(dist(t(d), method = "euclidean"), "complete"),
        legend_breaks = legend_labels %>% {log10(. + 1)},
        legend_labels = c(legend_labels[-length(legend_labels)], ">50"),

        annotation_col = annotation_col,
        annotation_row = annotation_row,
        annotation_colors = annotation_colors
      )
      p
    }
}
dev.off()
