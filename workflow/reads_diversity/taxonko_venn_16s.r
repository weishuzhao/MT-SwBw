###
#' @Date: 2022-03-03 15:56:58
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2023-09-08 18:27:31
#' @FilePath: /2021_09-MT10kSW/workflow/reads_diversity/taxonko_venn_16s.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, fig_out_path, keyword_args                         #####
fig_mag_venn <- argv[2]

##### LOAD data AND transform TO basic format                              #####
sample_meta_cross <-
  read.csv("results/reads_diversity/metadata.tsv", sep = "\t") %>%
  mutate(
    X = get("Sample"),
    Layer = ifelse(
      get("Type") == "16s", get("Sample"), gsub("^[^_]+_", "", get("Sample"))
    ),
    Sample = paste0(get("Group"), "_", get("Layer"))
  ) %>%
  merge(unique(sample_meta[c("Location", "Group")]))

otu_count <-
  "results/reads_diversity/level-7.csv" %>%
  {
    df <- read.csv(.)
    colnames(df) <-
      c("X", read.csv(., header = FALSE)[1, -1])
    df
  } %>%
  dplyr::select(
    -c("Layers", "Depth", "Latitude", "Longitude", "Location", "Group")
  ) %>%
  pivot_longer(!c("X"), names_to = "Taxonomy", values_to = "ReadsCount") %>%
  merge(sample_meta_cross) %>%
  filter(get("ReadsCount") > 0)

otu_rltabd <-
  otu_count %>%
  group_by(Sample = get("Sample")) %>%
  mutate(Abundance = get("ReadsCount") / sum(get("ReadsCount")) * 100)

### ######################################################################## ###
#### Define function AND Calculate data                                     ####
### ######################################################################## ###
##### Define function                                                      #####
#' DRAW Venn plot ACCORDING TO div.otu
#'
#' taxon that unannotated at given level
#'  will be regarded different across envirnment
#' e.g. if the taxon is "Bacteria;Proteobaceria;" after split,
#'      it will result in "Bacteria;Proteobaceria;env1"
#'       as it is obtained in env1
venn_grid_taxon <- function(taxon_level_spec) {
  otu_rltabd %>% # nolint: object_usage_linter.
    mutate(name = taxon.split(get("Taxonomy"), 1, taxon_level_spec)) %>% # nolint
    {
      split(.$name, .$Group) # nolint
    } %>%
    lapply(unique) %>%
    {
      VennDiagram::venn.diagram(
        x = ., # nolint
        filename = NULL, imagetype = "png",
        fill = sample_meta_col[names(.)], alpha = 0.75, # nolint
        lwd = 3,
        label.col = "black",
        cex = 2, fontfamily = "Arial", fontface = "bold",
        main = taxon_level_spec, main.cex = 2, main.fontfamily = "Arial",
        main.pos = c(0.5, 0), main.just = c(0.5, 1),
        # cat.col = {names(.) %>% },
        ext.line.lty = "dotted", ext.dist = -0.1,
        disable.logging = TRUE
      )
    } %>%
    {
      ggpubr::as_ggplot(.) # nolint
    } %>%
    {
      . + theme(plot.margin = unit(rep(0.3, 4), "in")) # nolint
    }
}

##### Calculate data                                                       #####

### ######################################################################## ###
#### Plot figures and OUTPUT                                                ####
### ######################################################################## ###
##### Plot figures                                                         #####
p_all <- NULL
for (taxon.level.spec in taxon.levels) {
  if (is.null(p_all)) {
    p_all <- venn_grid_taxon(taxon.level.spec)
  } else {
    p_all <- p_all + venn_grid_taxon(taxon.level.spec)
  }
}

##### OUTPUT                                                               #####
p <-
  p_all +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")
ggsave(
  filename = fig_mag_venn,
  plot = p,
  width = 10, height = 10, dpi = 300
)
