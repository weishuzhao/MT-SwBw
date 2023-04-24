###
#' @Date: 2022-05-18  17:46:00
#' @Editor: Wang Jing
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-08-05 22:17:50
#' @FilePath: /2021_09-MT10kSW/workflow/MAGs_tpm/ko_tpm_relationship_N.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)

library(igraph)


Wtdb_abd = argv[1]
#Wtdb_abd = stringr::str_glue("Stdb.relative_abundance.tsv") %>% file_path$file_path$results() %>% as.character


### ######################################################################## ###
#### loading data                                                           ####
### ######################################################################## ###
Stdb = load__Stdb()
Wtdb = load__Wtdb()
genome_taxonomy = load__genome_taxonomy(Stdb, Wtdb)
genome.relative_abundance = get_relative_abundance(Wtdb_abd, genome_taxonomy)

key_genes = load__key_genes()


### ######################################################################## ###
#### modify data for analysis                                               ####
### ######################################################################## ###
Stdb.spec =
  get_focus_species(genome.relative_abundance, Stdb, genome_taxonomy, 1) %>%
  # set order according to group distribution
  {.[.[unique(.$group)] %>% apply(1, . %>% paste(collapse = "")) %>% order, ]}

genomeko = load__genomeko()
genomeko.key =
  genomeko[key_genes$KO %>% as.character, ] %>%
  .[apply(., 1, sum) > 0, apply(., 2, sum) > 0]

### ######################################################################## ###
#### function for analysis                                                  ####
### ######################################################################## ###
genome_jaccard_network <- function(
  genomeko, ko.spec = NA, genome.spec = NA, threshold = 0.3
) {
  if (any(is.na(ko.spec))) ko.spec = rownames(genomeko)
  if (any(is.na(genome.spec))) genome.spec = colnames(genomeko)

  genomeko.key =
    genomeko %>%
    .[ko.spec, genome.spec] %>%
    .[apply(., 1, sum) > 0, apply(., 2, sum) > 0]

  key_genes = key_genes %>%
    {rownames(.) = .$KO; .} %>% .[rownames(genomeko.key), ]

  m.d = vegdist(genomeko.key, method = "jaccard", binary = TRUE)

  g =
    igraph::graph_from_adjacency_matrix(
      1 - as.matrix(m.d), weighted = T, mode = "undirected", diag = FALSE
    ) %>%
    delete.edges(E(.) %>% .[.$weight < threshold]) %>%
    delete.vertices(V(.)[degree(.) == 0])

  plot.igraph(
    g,
    layout = layout_with_fr,
    edge.width = E(g)$weight,
    vertex.size = 7,
    vertex.label = key_genes$label[match(V(g)$name, key_genes$KO)],
    vertex.label.cex = 0.7,
    vertex.color =
      paletteer::paletteer_dynamic("cartography::pastel.pal", 17)[
        key_genes$pathway[match(V(g)$name, key_genes$KO)]]
  )
  g
}

op <- par(mfrow = c(2, 2))
genome_jaccard_network(
  genomeko,
  key_genes %>% .[.$new == "False", "KO"] %>% as.character(.),
  Stdb.spec %>% .[.[c("Sw")] %>% paste(collapse = "") != "", "genome"]
)
genome_jaccard_network(
  genomeko,
  key_genes %>% .[.$new == "False", "KO"] %>% as.character(.),
  Stdb.spec %>% .[.[c("Ss")] %>% paste(collapse = "") != "", "genome"]
)
genome_jaccard_network(
  genomeko,
  key_genes %>% .[.$new == "False", "KO"] %>% as.character(.),
  Stdb.spec %>% .[.[c("Bw")] %>% paste(collapse = "") != "", "genome"]
)
genome_jaccard_network(
  genomeko,
  key_genes %>% .[.$new == "False", "KO"] %>% as.character(.),
  Stdb.spec %>% .[.[c("Bs")] %>% paste(collapse = "") != "", "genome"]
)
par(op)



genome_jaccard_network(
  genomeko,
  genome.spec = grepl("f__Nitrosopumilaceae", Stdb.spec$Taxonomy)
)


### ######################################################################## ###
#### depressed                                                              ####
### ######################################################################## ###
ko_jaccard <- function(genomeko, loc) {
  index = genomeko.key %>% apply(., 1, sum) > 0
  m.d = vegdist(genomeko[index, ], method = "jaccard", binary = TRUE)

  dend <- m.d %>%
    hclust(method = "complete") %>%
    as.dendrogram() %>%
    dendextend::set("labels", paste(ko_label$Pathway[index], ko_label$Lable[index], sep = "-"))

  dend
}


ko_label =
  read.csv("workflow/MAGs_tpm/key_gene.list", sep = "\t") %>%
  mutate(Pathway = factor(Pathway, levels = c(
    "WL", "3HP-4HB", "rTCA", "CBB",
    "Alkane", "Aromatic", "Complex sugar",
    "D-AA", "L-sugar",
    "N", "S",
    "As", "Se",
    "ROS", "Glycine betaine", "TMAO", "DMSP"
  ))) %>%
  arrange(Pathway, Lable)

genomeko =
  "workflow/MAGs/genomeko1.csv" %>%
  {
    df = read.csv(., row.names = 1)
    colnames(df) = read.csv(., header = FALSE)[1, -1]
    df
  } %>%
  .[ko_label$Class, Wtdb$genome] %>%
  t %>%
  data.frame(
    .,
    genome = Wtdb$genome,
    location = grepl("^1", Wtdb$genome) %>% ifelse("ME", "MT")
  ) %>%
  as_tibble()


ko_ME <- ko_jaccard(genomeko, "ME", 0.3, "figs/5B1.pdf")
ko_MT <- ko_jaccard(genomeko, "MT", 0.3, "figs/5B2.pdf")

dends <- dendlist(ko_MT, ko_ME)
tanglegram(
  dends,
  main_left = "ME", main_right = "MT",
  common_subtrees_color_branches = TRUE,
  lwd = 2, edge.lwd = 1.5,
  columns_width = c(5,1,5),
  margin_inner = 9)


library(dendextend)
