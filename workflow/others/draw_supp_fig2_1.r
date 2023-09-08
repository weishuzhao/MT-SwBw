###
#' @Date: 2022-07-20 13:43:25
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2023-09-08 17:40:26
#' @FilePath: /2021_09-MT10kSW/workflow/others/draw_supp_fig2_1.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, fig_out_path, keyword_args                         #####
wtdb_abd <- argv[1]
# wtdb_abd = stringr::str_glue("Wtdb.relative_abundance.tsv") %>% file_path$file_path$results() %>% as.character

##### GLOBAL CONST vars                                                    #####
font_size_1 <- 13
font_size_2 <- 10
font_size_3 <- 10
axis_ticks_length <- 0.1

##### LOAD data AND transform TO basic format                              #####
genome_taxonomy <- load__genome_taxonomy(load__Stdb(), load__Wtdb())
genome_rltabd <- get_relative_abundance(wtdb_abd, genome_taxonomy)
sample_meta_cross <-
  read.csv("results/reads_diversity/metadata.tsv", sep = "\t") %>%
  mutate(
    X = get("Sample"),
    Layer = ifelse(
      get("Type") == "16s", get("Sample"), gsub("^[^_]+_", "", get("Sample"))
    ),
    Sample = paste0(get("Group"), "_", get("Layer"))
  ) %>%
  merge(unique(sample_meta[c("Location", "Group")])) %>%
  mutate(Group = factor(get("Group"), names(sample_meta_col))) %>%
  .[order(.$Group, .$Sample), ]

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

##### Calculate data                                                       #####
p_nmds_s <-
  list(
    "16S" = otu_rltabd,
    "metagenome" = mutate(genome_rltabd, Abundance = get("Relative_abundance"))
  ) %>%
  lapply(function(div_otu_long_raw) {
    div_otu_long <-
      div_otu_long_raw[c("Sample", "Group", "Taxonomy", "Abundance")]
    taxon.levels[-1] %>%
      as.list() %>%
      `names<-`(unlist(.)) %>%
      lapply(function(taxon_level_spec) {
        div_otu <-
          div_otu_long %>%
          mutate(
            name = taxon.split(get("Taxonomy"), 1, taxon_level_spec)
          ) %>%
          reshape2::acast(
            formula = formula("name ~ Sample"),
            value.var = "Abundance",
            fun.aggregate = sum, fill = 0
          )
        list("jaccard" = "jaccard", "bray" = "bray") %>%
          lapply(function(dist) {
            p <-
              plot.beta.div(div_otu,
                pname = "relative abundance",
                method = "nmds", dist = dist, area = "polygon"
              ) +
              scale_color_manual(values = sample_meta_col) +
              scale_fill_manual(values = sample_meta_col)
            p$labels$title <-
              p$labels$title %>%
              gsub("^[^\n]+\n(.+) (p\\(Pr\\(>F\\)\\)[^\n]+)\n.+$", "\\1\n\\2", .)
            p
          })
      })
  })


### ######################################################################## ###
#### Plot figures and OUTPUT                                                ####
### ######################################################################## ###
p <- NULL %>%
  {
    p <- .
    for (taxon_level_spec in taxon.levels[-1]) {
      for (dist in c("jaccard", "bray")) {
        for (data_type in c("16S", "metagenome")) {
          if (is.null(p)) {
            p <- p_nmds_s[[data_type]][[taxon_level_spec]][[dist]] +
              labs(x = "", y = taxon_level_spec)
          } else if (
            (dist == "jaccard") &&
              (data_type == "16S") &&
              (taxon_level_spec != taxon.levels[length(taxon.levels)])) {
            p <- p + p_nmds_s[[data_type]][[taxon_level_spec]][[dist]] +
              labs(x = "", y = taxon_level_spec)
          } else if (
            (dist == "jaccard") &&
              (data_type == "16S") &&
              (taxon_level_spec == taxon.levels[length(taxon.levels)])) {
            p <- p + p_nmds_s[[data_type]][[taxon_level_spec]][[dist]] +
              labs(x = paste(dist, data_type), y = taxon_level_spec)
          } else if ((taxon_level_spec == taxon.levels[length(taxon.levels)])) {
            p <- p + p_nmds_s[[data_type]][[taxon_level_spec]][[dist]] +
              labs(x = paste(dist, data_type), y = "")
          } else {
            p <- p + p_nmds_s[[data_type]][[taxon_level_spec]][[dist]] +
              labs(x = "", y = "")
          }
        }
      }
    }
    p +
      plot_layout(ncol = 4, guides = "collect") +
      plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")
  }

ggsave("results/figs/figs2_relative_nmds_all.svg", p, width = 14, height = 18)

### ######################################################################## ###
#### Plot figures and OUTPUT                                                ####
### ######################################################################## ###
##### Plot figures                                                         #####
cross_rltabd <-
  read.csv("results/reads_diversity/abundance.csv") %>%
  merge(sample_meta_cross[c("X", "Sample")]) %>%
  mutate(X = NULL) %>%
  column_to_rownames("Sample")

p2 <-
  #' otu_rltabd %>%
  #' mutate(Taxa_label = taxon.split(get("Taxonomy"), 1, 3)) %>%
  cross_rltabd %>%
  rownames_to_column("Sample") %>%
  pivot_longer(
    !c("Sample"),
    names_to = "Taxa_label", values_to = "Abundance"
  ) %>%
  merge(sample_meta_cross) %>%
  filter(get("Type") == "16S") %>%
  as.data.frame() %>%
  get_percent_plot(
    "Abundance",
    fill.name = "Taxa_label",
    sample.name = "Layer",
    labs.x = "16S sample", labs.y = "relative abundance",
    TOP_N_TAXON_PER_LAYER = 10,
    font_size_1 = font_size_1, font_size_2 = font_size_2,
    font_size_3 = font_size_3, axis.ticks.length = axis_ticks_length
  )
# p2$data[c("Group", "Site", "Layer", "Sample", "Abundance", "name", "Genome", "Taxonomy")] %>% write.csv("results/figs/fig2_relative_nmds_class.C.csv", row.names = FALSE, quote = FALSE)
p2_x <- p2 +
  scale_x_discrete(
    limits = filter(sample_meta_cross, get("Type") == "16S")$Layer
  ) +
  theme(
    axis.text.x = element_text(
      color = sample_meta_col[
        filter(sample_meta_cross, get("Type") == "16S")$Group
      ]
    )
  )

##### OUTPUT                                                               #####
ggsave(
  filename = "results/figs/figs2_class_16s.svg",
  plot = p2_x, width = 6, height = 7
)
