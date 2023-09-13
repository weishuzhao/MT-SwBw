###
#' @Date: 2023-09-12 23:12:57
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2023-09-12 23:20:19
#' @FilePath: /2021_09-MT10kSW/workflow/others/draw_supp_fig_1_2.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, fig_out_path, keyword_args                         #####
wtdb_abd <- argv[1]
#' wtdb_abd = stringr::str_glue("Wtdb.relative_abundance.tsv") %>% file_path$file_path$results() %>% as.character # nolint
otu_count_file <- argv[2]
#' otu_count_file <- "results/reads_diversity/level-7.csv"
fig_out <- argv[3]

##### GLOBAL CONST vars                                                    #####
share_col <-
  c("Water only" = "#66E5FF", "Share" = "#91B2BE", "Sediment only" = "#DBD4C5")

# depressed
group1_col <- c(
  "Bw|Sw",
  "Bw, Sw",
  "Bw, Ss", "Bs, Bw",
  "Bs, Sw", "Ss, Sw",
  "Bs, Bw, Sw", "Bw, Ss, Sw",
  "Bs, Bw, Ss, Sw",
  "Bs, Bw, Ss", "Bs, Ss, Sw",
  "Bs, Ss",
  "Bs|Ss"
) %>%
  {
    col <- c(paletteer::paletteer_d("ggsci::default_igv")[seq_along(.)])
    names(col) <- .
    col
  }

font_size_1 <- 13
font_size_2 <- 10
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
  merge(unique(sample_meta[c("Location", "Group")])) %>%
  mutate(Group = factor(get("Group"), names(sample_meta_col))) %>%
  .[order(.$Group, .$Sample), ] %>%
  mutate(Group = as.character(get("Group")))

otu_count <-
  otu_count_file %>%
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
assign_share <-
  . %>%
  gsub("^(Bw|Sw)$", "Water only", .) %>%
  gsub("^(Bs|Ss)$", "Sediment only", .) %>%
  unique() %>%
  sort() %>%
  paste(collapse = ", ") %>%
  gsub("^(.+, .+)$", "Share", .) %>%
  factor(names(share_col) %>% rev())

##### Calculate data                                                       #####
otu_rltabd <-
  otu_count %>%
  group_by(Sample = get("Sample")) %>%
  mutate(Abundance = get("ReadsCount") / sum(get("ReadsCount")) * 100)

### ######################################################################## ###
#### Plot figures and OUTPUT                                                ####
### ######################################################################## ###
##### Plot figures                                                         #####

p1 <-
  genome_rltabd %>%
  filter(get("Relative_abundance") > 0.01) %>%
  group_by(Genome = get("Genome")) %>%
  mutate(Group1 = assign_share(get("Group"))) %>%
  group_by(
    Group1 = get("Group1"),
    Site = get("Site"), Layer = get("Layer"),
    Sample = get("Sample"), Group = get("Group")
  ) %>%
  summarise(Abundance = sum(get("Relative_abundance"))) %>%
  as.data.frame() %>%
  get_percent_plot(
    "Abundance",
    fill.name = "Group1", sample.name = "Layer",
    labs.y = "MAG prevalence"
  ) +
  scale_x_discrete(
    limits =
      sample_meta[c("Site", "Layers")] %>%
      apply(1, . %>% paste(collapse = ".."))
  ) +
  theme(axis.text.x = element_text(color = sample_meta_col[sample_meta$Group]))
p2 <-
  genome_rltabd %>%
  filter(get("Relative_abundance") > 0.01) %>%
  group_by(Genus = taxon.split(get("Taxonomy"), 1, 6)) %>%
  mutate(Group1 = assign_share(get("Group"))) %>%
  group_by(
    Group1 = get("Group1"),
    Site = get("Site"), Layer = get("Layer"),
    Sample = get("Sample"), Group = get("Group")
  ) %>%
  summarise(Abundance = sum(get("Relative_abundance"))) %>%
  as.data.frame() %>%
  get_percent_plot(
    "Abundance",
    fill.name = "Group1", sample.name = "Layer",
    labs.y = "MAG genus prevalence"
  ) +
  scale_x_discrete(
    limits =
      sample_meta[c("Site", "Layers")] %>%
      apply(1, . %>% paste(collapse = ".."))
  ) +
  theme(axis.text.x = element_text(color = sample_meta_col[sample_meta$Group]))

p3 <- otu_rltabd %>%
  filter(get("Abundance") > 0.01) %>%
  group_by(Genus = taxon.split(get("Taxonomy"), 1, 6)) %>%
  mutate(Group1 = assign_share(get("Group"))) %>%
  group_by(
    Group1 = get("Group1"),
    Layer = get("Layer"),
    Sample = get("Sample"), Group = get("Group")
  ) %>%
  summarise(Abundance = sum(get("Abundance"))) %>%
  as.data.frame() %>%
  get_percent_plot(
    "Abundance",
    fill.name = "Group1", sample.name = "Layer",
    labs.y = "16S genus prevalence"
  ) +
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


p1 + p2 + p3 +
  plot_layout(design = "A\nB\nC", guides = "collect") &
  scale_fill_manual(values = share_col)

##### OUTPUT                                                               #####
ggsave(filename = fig_out, plot = p2_x, width = 6, height = 10)
