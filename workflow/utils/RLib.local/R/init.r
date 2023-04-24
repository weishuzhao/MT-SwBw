###
#' @Date: 2022-06-09 22:55:53
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-09-22 16:59:39
#' @FilePath: /2021_09-MT10kSW/workflow/utils/RLib.local/R/init.r
#' @Description:
###

suppressWarnings(library(tidyverse))
suppressWarnings(library(patchwork))
suppressWarnings(library(vegan))


if(grepl("RLib.local$", getwd())) {

  #### LOAD general functions                                               ####
  . <- plyr::.
  message.print = function(...) {capture.output(...) %>% paste0(collapse = "\n") %>% {base::message(.)}}


  #### LOAD general vars                                                    ####
  argv = commandArgs(trailingOnly = TRUE)

  file_path =
    setwd("../../../") %>%
    {
      reticulate::use_condaenv("python39", Sys.getenv("CONDA_EXE", "~/software/miniconda3/bin/conda"))
      file_path = reticulate::import_from_path("workflow.utils.file_path")
      setwd(.)
      file_path
    }

  ##### sample meta AND color                                              #####
  paletteer::paletteer_d("RColorBrewer::Paired")[c(2, 1, 6, 5)]
  sample_meta_col =
    c("Sw" = "#78c86c", "Ss" = "#fc8715", "Bw" = "#2b83ba", "Bs" = "#d7191c")

  sample_meta =
    file_path$sample_meta %>%
    split(.$Group) %>%
    {lapply(names(sample_meta_col), function(x) .[[x]])} %>%
    bind_rows

  site_group =
    sample_meta %>% dplyr::select(c("Site", "Group")) %>% unique


  #### LOADER function                                                      ####
  layer_2_group <- function(layer) {
    site = gsub("^(.+)\\.\\..+$", "\\1", layer)
    group = site_group %>% {rownames(.) = .$Site; .} %>% .[site, "Group"]
  }

  get_focus_species <- function(
      genome.relative_abundance, Stdb, genome_taxonomy, taxon.level.spec
  ) {
    annotation_species =
      get_taxon_group(genome.relative_abundance, "Genome", "Group") %>%
      {.$Genome = .$Genome %>% gsub("^(.+)$", "\\1.fa", .); .}

    genome.relative_abundance %>%
      {.$name = taxon.split(.$Taxonomy, 1, taxon.level.spec); .} %>%
      {.[c("Group", "name")]} %>%
      unique %>%
      {.$Location = .$Group %>% grepl("s$", .) %>% ifelse("sediment", "water"); .} %>%
      {.[c("Location", "name")]} %>%
      unique %>% .$name %>% .[duplicated(.)] %>% sort %>% unique %>%
      {
        focus_species = .
        genome_taxonomy %>%
          .[taxon.split(.$Taxonomy, 1, taxon.level.spec) %>% {. %in% focus_species}, ]
      } %>%
      {.$name.spec = taxon.split(.$Taxonomy, 1, taxon.level.spec) %>% as.character; .} %>%
      {.$name = taxon.split(.$Taxonomy, 1, 7) %>% as.character; .} %>%
      {.$Genome = .$Genome %>% gsub("^(.+)$", "\\1.fa", .); .} %>%
      merge(Stdb) %>%
      {.$quality = assign_quality(.$Completeness, .$Contamination); .} %>%
      {merge(., annotation_species, by = "Genome")}
  }

  load__gene_ko_tpm  <- function() {
    file_path$file_path$results("gene_ko_tpm.csv") %>% as.character %>%
      read.csv(col.names = c("KO", "Genome", "Layer", "TPM")) %>%
      {.$Layer = gsub("^.+\\.\\.(.+\\.\\..+).bam", "\\1", .$Layer); .}
  }


  load__genomeko <- function(genome_taxonomy = NA) {
    gene_map =
      file_path$file_path$results("gene_map.csv") %>% as.character %>%
      read.csv %>% {colnames(.) = colnames(.) %>% str_to_title; .} %>%
      .[c("Genome", "Gene")]
    all_gene_annots =
      file_path$file_path$results("all_gene_annots.csv") %>% as.character %>%
      read.csv %>% {colnames(.) = c("All", "KO"); .} %>%
      .[c("KO", "All")]
    genomeko.long = merge(gene_map, all_gene_annots, by.x = "Gene", by.y = "All")
    if (is.na(genome_taxonomy)) {
      genomeko.long %>%
        {.$Genome = .$Genome %>% gsub(".fa$", "", .); .} %>%
        reshape2::acast(formula("KO ~ Genome"))
    } else {
      genomeko.long %>%
        {.$Genome = .$Genome %>% gsub(".fa$", "", .); .} %>%
        merge(genome_taxonomy[c("Cluster", "Genome")]) %>%
        {.$Genome <- NULL; .} %>%
        merge(
          genome_taxonomy[c("Cluster", "Genome", "Rep")] %>%
          .[.$Rep, ]
        )
    } %>%
    reshape2::acast(formula("KO ~ Genome"))
  }

  load__key_genes <- function() {
    "data/key_gene.tsv" %>%
      read.csv(sep = "\t") %>%
      {.$KO = .$KO %>% factor(levels = unique(.)); .} %>%
      {.$Label = .$Label %>% factor(levels = unique(.)); .} %>%
      {.$Pathway = .$Pathway %>% factor(levels = unique(.)); .}
  }

  set__pathway_col <- function(pathways, n = 24) {
    if (is.null(levels(pathways))) {
      pathways = unique(pathways)
    } else {
      pathways = levels(pathways)
    }

    cols = c(paletteer::paletteer_d("ggsci::default_igv")[1:length(pathways)])
    names(cols) = pathways
    cols
  }

  load__taxon_color <- function() {
    "data/taxon_color.tsv" %>%
      read.csv(sep = "\t") %>%
      {.$Taxa_label = .$LEGEND_LABELS %>% taxon.split(1, 2) %>% as.character; .} %>%
      split(.$Taxa_label) %>%
      lapply(. %>% {
        .$Taxa_label =
          if (length(.$Taxa_label) > 1) ifelse(
            .$LEGEND_LABELS == .$Taxa_label,
            paste0(.$LEGEND_LABELS, ";others"),
            .$LEGEND_LABELS
          ) else .$LEGEND_LABELS
        .
      }) %>%
      bind_rows %>%
      {rownames(.) = .$Taxa_label; .} %>%
      .[order(.$LEGEND_LABELS, decreasing = TRUE), c("LEGEND_COLORS", "LEGEND_LABELS", "Taxa_label")]
  }

  get_taxon_color <- function(taxonomy) {
    apply(
      taxon_color, 1,
      . %>% {ifelse(grepl(pattern = .[2], x = taxonomy), .[3], "")}
    ) %>% {.[. != ""][1]} %>% {ifelse(is.na(.), "others", .)}
  }

} else {

  make.package <- function(dir) {
    suppressWarnings(roxygen2::roxygenize(dir))
    suppressWarnings(install.packages(dir, quiet = TRUE, repos = NULL))
  }

  make.package("../../RLib")
  library(RLib)
  make.package("../../RLib.local")
  library(RLib.local)

}
