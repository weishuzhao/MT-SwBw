###
#' @Date: 2022-07-31 20:51:00
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-08-05 22:17:50
#' @FilePath: /2021_09-MT10kSW/workflow/MAGs_tpm/ko_tpm_MAG_N-layer_table.r
#' @Description: Genomes grouped by prevalence
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


Wtdb_abd = argv[1]
#Wtdb_abd = stringr::str_glue("Stdb.relative_abundance.tsv") %>% file_path$file_path$results() %>% as.character
tab_out = argv[2]
#tab_out = ""


### ######################################################################## ###
#### loading data                                                           ####
### ######################################################################## ###
Stdb = load__Stdb()
genome_taxonomy = load__genome_taxonomy(load__Stdb(), load__Wtdb())
genome.relative_abundance = get_relative_abundance(Wtdb_abd, genome_taxonomy)

gene_ko_tpm = load__gene_ko_tpm()

key_genes = load__key_genes()

### ######################################################################## ###
#### modify data for analysis                                               ####
### ######################################################################## ###
##### from workflow/MAGs/gmodule_heatmap_venn.r                            #####
Wtdb.existance =
  genome.relative_abundance %>%
  {data.frame("genome" = .$genome %>% paste0(".fa"),
              "Taxonomy" = .$Taxonomy %>% taxon.split(1, 7),
              "Layer" = .$Layer, "sample" = .$sample,
              "exist" = TRUE)}

key_gene_genome_layer =
  key_genes %>%
  .[.$pathway == "N", ] %>%
  {.$ko = .$KO; .$KO <- NULL; .} %>%
  merge(gene_ko_tpm) %>%
  merge(Wtdb.existance[c("Layer", "sample")] %>% unique) %>%
  merge(Wtdb.existance[c("genome", "Layer", "exist")], all.x = TRUE) %>%
  {.$tpm.adj = .$tpm * (1 - 2 * is.na(.$exist)); .}


### ######################################################################## ###
#### OUTPUT                                                                 ####
### ######################################################################## ###
key_gene_genome_layer %>%
  merge(Wtdb.existance[c("genome", "Taxonomy")] %>% unique) %>%
  {.$max = 0; .} %>%
  reshape2::dcast(formula("arrow + ko + label + genome + Taxonomy + max ~ sample"),
                  value.var = "tpm.adj") %>%
  {.$max = apply(.[-(1:5)], 1, . %>% .[!is.na(.)] %>% abs %>% max); .} %>%
  split(.$arrow) %>%
  lapply(function(df) {df$arrow <- NULL; df}) %>%
  lapply(function(df) df[-(1:5)] %>% abs %>% {.[is.na(.)] = 0; .} %>% apply(2, mean) %>% bind_rows(df)) %>%
  lapply(function(df) {df[1, "genome"] = "mean (abs)"; df}) %>%
  {
    cat(.[[1]] %>% colnames, sep = "\t", file = tab_out)
    cat("\n", file = tab_out, append = TRUE)
    for (group in names(.)) {
      cat(group, file = tab_out, append = TRUE)
      cat("\n", file = tab_out, append = TRUE)
      .[[group]] %>%
        write.table(file = tab_out, sep = "\t", quote = FALSE, na = "",
                    row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }
