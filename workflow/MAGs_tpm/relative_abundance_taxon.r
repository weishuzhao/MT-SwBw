###
#' @Date: 2022-06-25 10:52:06
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2024-02-19 14:27:50
#' @FilePath: /2021_09-MT10kSW/workflow/MAGs_tpm/relative_abundance_taxon.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, fig_out_path, keyword_args                         #####
Wtdb_abd <- argv[1]
# Wtdb_abd = stringr::str_glue("Wtdb.relative_abundance.tsv") %>% file_path$file_path$results() %>% as.character
fig_out <- argv[2]


##### GLOBAL CONST vars                                                    #####
font_size_1 <- 15
font_size_2 <- 13
font_size_3 <- 10
axis_ticks_length <- 0.1


##### LOAD data AND transform TO basic format                              #####
genome_taxonomy <- load__genome_taxonomy(load__Stdb(), load__Wtdb())
genome_rltabd <- get_relative_abundance(Wtdb_abd, genome_taxonomy)

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

### ######################################################################## ###
#### Define function AND Calculate data                                     ####
### ######################################################################## ###
##### Define function                                                      #####

##### Calculate data                                                       #####
genome_rltabd %>%
  mutate(FakeCount = ceiling(get("Relative.abundance") * 1000)) %>%
  pivot_wider(
    id_cols = "Taxonomy", names_from = "Sample",
    values_from = "FakeCount", values_fill = 0
  ) %>%
  column_to_rownames("Taxonomy") %>%
  t() %>%
  vegan::estimateR() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  dplyr::arrange(get("Sample"))

otu_count %>%
  pivot_wider(
    id_cols = "Taxonomy", names_from = "Sample",
    values_from = "ReadsCount", values_fill = 0
  ) %>%
  column_to_rownames("Taxonomy") %>%
  t() %>%
  vegan::estimateR() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  dplyr::arrange(get("Sample")) %>%
  mutate(
    Area = ifelse(grepl("^(S[^_])_.+$", get("Sample")), "slope", "bottom")
  ) %>%
  wilcox.test(formula("S.chao1 ~ Area"), data = ., paired = FALSE)

### ######################################################################## ###
#### Plot figures and OUTPUT                                                ####
### ######################################################################## ###
##### Plot figures                                                         #####
p1 <- get_total_plot(genome_rltabd, "Relative_abundance",
  labs.x = "total mapping rate"
)

taxon_color <- load__taxon_color()

p2 <-
  genome_rltabd %>%
  mutate(
    Taxonomy = sapply(taxon.split(get("Taxonomy"), 1, 7), get_taxon_color)
  ) %>%
  get_percent_plot("Relative_abundance", "Taxonomy",
    labs.x = "Layer", labs.y = "relative abundance"
  ) +
  scale_fill_manual(
    values = c(taxon_color$LEGEND_COLORS) %>%
      `names<-`(taxon_color$LEGEND_LABELS) %>%
      .[order(names(.))]
  )


##### OUTPUT                                                               #####
pout <- p1 + p2 + guides(fill = "none") +
  plot_layout(heights = c(1, 4), nrow = 2, guides = "collect") +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 18))
ggsave(filename = fig_out, plot = pout, width = 14, height = 10, dpi = 300)
