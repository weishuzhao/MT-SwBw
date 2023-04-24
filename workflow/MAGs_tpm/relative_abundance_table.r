###
#' @Date: 2022-07-31 20:51:00
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-08-05 14:44:38
#' @FilePath: /2021_09-MT10kSW/workflow/MAGs_tpm/relative_abundance_table.r
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

### ######################################################################## ###
#### modify data for analysis                                               ####
### ######################################################################## ###
##### from workflow/MAGs/gmodule_heatmap_venn.r                            #####
Stdb.spec =
  get_focus_species(genome.relative_abundance, Stdb, genome_taxonomy, 1) %>%
  # set order according to group distribution
  {.[.[unique(.$Group)] %>% apply(1, . %>% paste(collapse = "")) %>% order, ]}


### ######################################################################## ###
#### OUTPUT                                                                 ####
### ######################################################################## ###
Stdb.spec %>%
  #{split(.[c("Genome", "Taxonomy")], .[unique(.$Group)])} %>%
  split(.[unique(.$Group)] %>% apply(1, . %>% paste(collapse = "."))) %>%
  {
    div.otu =
      genome.relative_abundance %>%
      reshape2::dcast(formula("Genome ~ Sample"),
                      value.var = "Relative_abundance",
                      fun.aggregate = sum) %>%
      {.$Genome = paste0(.$Genome, ".fa"); .}
    lapply(., . %>% .[c("Genome", "Taxonomy")] %>% merge(div.otu))
  } %>%
  lapply(function(df) df[c(-1, -2)] %>% apply(2, mean) %>% bind_rows(df)) %>%
  lapply(function(df) {df[1, "Genome"] = "mean"; df}) %>%
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
