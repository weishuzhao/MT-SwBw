###
#' @Date: 2022-05-18  17:46:00
#' @Editor: Wang Jing
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-09-01 23:47:15
#' @FilePath: /2021_09-MT10kSW/workflow/MAGs/ko_relationship_N.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)

suppressMessages(library(igraph))


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, keyword_args, fig_out_path                         #####
Wtdb_abd = argv[1]
#Wtdb_abd = stringr::str_glue("Stdb.relative_abundance.tsv") %>% file_path$file_path$results() %>% as.character
fig_out = argv[2]
domain.spec = argv[3]
#domain.spec = "Archaea"  # domain.spec = "Bacteria"  # domain.spec = "All"

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
genome_jaccard_network <- function(
  genomeko,
  ko.spec = NA, genome.spec = NA,
  threshold = 0.3, key_genes = key_genes
) {
  if (any(is.na(ko.spec))) ko.spec = rownames(genomeko)
  else ko.spec = ko.spec %>% .[. %in% rownames(genomeko)]
  if (any(is.na(genome.spec))) genome.spec = colnames(genomeko)
  else genome.spec = genome.spec %>% .[. %in% colnames(genomeko)]

  genomeko.key =
    genomeko %>%
    .[ko.spec, genome.spec] %>%
    .[apply(., 1, sum) > 0, apply(., 2, sum) > 0]

  key_genes.spec = key_genes %>%
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
    vertex.size = 18,
    vertex.label = key_genes.spec$Label[match(V(g)$name, key_genes.spec$KO)],
    vertex.label.cex = 0.8,
    vertex.label.font = 2,
    vertex.label.family = "Arial",
    vertex.color =
      pathway_col[key_genes.spec$Pathway[match(V(g)$name, key_genes.spec$KO)]],
    vertex.frame.color = "gray",
    asp = 0
  )
}

##### Calculate data                                                       #####
Stdb.domain.spec =
  genome.relative_abundance %>%
  {
    if (domain.spec == "All") .
    else .[.$Taxonomy %>% taxon.split(1, 1) %>% {. == domain.spec}, ]
  } %>%
  split(.$Cluster) %>%
  lapply(. %>% {.["Genome"] = unique(.[.$Rep == TRUE, "Genome"]); .}) %>%
  bind_rows() %>%
  get_taxon_group %>%
  {.[.[unique(.$Group)] %>% apply(1, . %>% paste(collapse = "")) %>% order, ]}

### ######################################################################## ###
#### Plot figures and OUTPUT                                                ####
### ######################################################################## ###
svg(filename = fig_out, width = 10, height = 7)
par(mfrow = c(2, 2)) %>% {
  for (i in names(sample_meta_col)) {
    genome_jaccard_network(
      genomeko.key,
      genome.spec =
        Stdb.domain.spec %>%
        .[.[i] %>% apply(1, . %>% {any(. != "")}), ] %>%
        .[, "Genome"],
      key_genes = key_genes_N
    )
    title(i)
  }
  par(.)
}
dev.off()
