###
#' @Date: 2022-07-31 20:51:00
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-08-05 14:46:52
#' @FilePath: /2021_09-MT10kSW/workflow/MAGs_tpm/ko_tpm_MAG_N_table.r
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
  {data.frame("Genome" = .$Genome %>% paste0(".fa"),
              "Taxonomy" = .$Taxonomy %>% taxon.split(1, 7),
              "Layer" = .$Layer, "Group" = .$Group,
              "Exist" = TRUE)}

key_gene_genome_layer =
  key_genes %>%
  .[.$Pathway == "N", ] %>%
  merge(gene_ko_tpm) %>%
  merge(Wtdb.existance[c("Layer", "Group")] %>% unique) %>%
  merge(Wtdb.existance[c("Genome", "Layer", "Exist")], all.x = TRUE) %>%
  {.$TPM.adj = .$TPM * (1 - 2 * is.na(.$Exist)); .}


### ######################################################################## ###
#### OUTPUT                                                                 ####
### ######################################################################## ###
key_gene_genome_layer %>%
  merge(Wtdb.existance[c("Genome", "Taxonomy")] %>% unique) %>%
  {.$sum = 0; .} %>%
  reshape2::dcast(formula("Genome + Taxonomy + Group + Layer + sum ~ KO"),
                  value.var = "TPM.adj") %>%
  {.$sum = apply(.[-(1:5)], 1, . %>% .[!is.na(.)] %>% abs %>% sum); .} %>%
  split(.$Group) %>%
  #lapply(function(df) {df$Group <- NULL; df}) %>%
  lapply(function(df) {df[order(df$Layer, df$Taxonomy, df$sum),]}) %>%
  lapply(function(df) df[-(1:5)] %>%
           apply(2, . %>% .[!is.na(.)] %>% abs %>% max) %>%
           bind_rows(df)) %>%
  lapply(function(df) {df[1, "Genome"] = "max"; df}) %>%
  {
    cat(.[[1]] %>% colnames, sep = "\t", file = tab_out)
    cat("\n", file = tab_out, append = TRUE)
    rownames(key_genes) = key_genes$KO
    key_genes[.[[1]] %>% colnames, c("Label", "Arrow")] %>% t %>%
      write.table(file = tab_out, sep = "\t", quote = FALSE, na = "",
                  row.names = FALSE, col.names = FALSE, append = TRUE)
    for (group in names(.)) {
      cat(group, file = tab_out, append = TRUE)
      cat("\n", file = tab_out, append = TRUE)
      .[[group]] %>%
        write.table(file = tab_out, sep = "\t", quote = FALSE, na = "",
                    row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }
