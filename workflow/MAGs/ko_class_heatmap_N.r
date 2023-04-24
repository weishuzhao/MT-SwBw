###
#' @Date: 2022-05-18  17:46:00
#' @Editor: Wang Jing
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-09-22 17:04:35
#' @FilePath: /2021_09-MT10kSW/workflow/MAGs/ko_class_heatmap_N.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, keyword_args, fig_out_path                         #####
Wtdb_abd = argv[1]
#Wtdb_abd = stringr::str_glue("Stdb.relative_abundance.tsv") %>% file_path$file_path$results() %>% as.character
fig_out = argv[2]

##### LOAD data AND transform TO basic format                              #####
Stdb = load__Stdb()
Wtdb = load__Wtdb()
genome_taxonomy = load__genome_taxonomy(Stdb, Wtdb)
genome.relative_abundance = get_relative_abundance(Wtdb_abd, genome_taxonomy)

key_genes = load__key_genes()
key_genes_N =
  key_genes %>% .[.$Pathway == "N",] %>%
  {.$Pathway = .$Arrow %>% factor(levels = unique(.)); .}

pathway_col = set__pathway_col(key_genes_N$Pathway)

genomeko = load__genomeko()
genomeko.key =
  genomeko[key_genes_N$KO %>% as.character, ] %>%
  .[apply(., 1, sum) > 0, apply(., 2, sum) > 0]


### ######################################################################## ###
#### Define function AND Calculate data                                     ####
### ######################################################################## ###
##### Define function                                                      #####

##### Calculate data                                                       #####
annotation_class =
  genome.relative_abundance %>%
  {.$name = taxon.split(.$Taxonomy, 1, 3); .} %>%
  get_taxon_group("name") %>%
  {rownames(.) = .$name; .$name <- NULL; .}

classlocko_pct =
  Wtdb %>%
  {.$Genome = .$Genome %>% gsub(".fa$", "", .); .} %>%
  .[.$Genome %in% colnames(genomeko.key), ] %>%
  {split(genomeko.key[, .$Genome] %>% t %>% data.frame, .$Taxonomy %>% taxon.split(1, 3))} %>%
  lapply(. %>% apply(2, . %>% {sum(. > 0) / length(.)})) %>%
  bind_rows(.id = "Class") %>% column_to_rownames("Class")


### ######################################################################## ###
#### Plot figures and OUTPUT                                                ####
### ######################################################################## ###
##### Plot figures                                                         #####
p =
  classlocko_pct %>%
  {
    annotation_col =
      key_genes_N %>%
      {
        data.frame(
          Pathway = .$Pathway, Label = .$Label, row.names = .$KO
        )
      }
    d = .
    pheatmap::pheatmap(
      d %>% {.[. == 0] = NA; .},
      cellwidth = 10, cellheight = 10,
      cluster_cols = FALSE, cluster_rows = FALSE,
      na_col = "white",
      color = RColorBrewer::brewer.pal(n = 9, name = "OrRd")[c(3, 5, 8)],
      breaks = c(0, 0.2, 0.5, 1),
      legend_breaks = c(0, 0.2, 0.5, 1),
      gaps_col = annotation_col["Pathway"] %>% table %>% .[. > 0] %>% cumsum,
      annotation_row = annotation_class,
      annotation_col = annotation_col[c("Pathway", "Label")] %>% {.["Label"] = NULL; .},
      labels_col = annotation_col[, "Label"] %>% as.character,
      annotation_colors =
        append(
          sample_meta_col %>%
            Map(function(x, i) {c(i, "#FFFFFF") %>% {names(.) = c(x, ""); .}}, names(.), .),
          list(
            "Group" = sample_meta_col,
            "Pathway" = pathway_col[annotation_col$Pathway %>% unique]
          )
        ),
      silent = TRUE
    ) %>%
      ggplotify::as.ggplot(.)
  }

##### OUTPUT                                                               #####
ggsave(filename = fig_out, plot = p, width = 16, height = 9)
