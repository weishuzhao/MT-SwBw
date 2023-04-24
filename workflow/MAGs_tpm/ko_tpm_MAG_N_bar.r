###
#' @Date: 2022-05-04 09:42:54
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-08-05 14:59:01
#' @FilePath: /2021_09-MT10kSW/workflow/MAGs_tpm/ko_tpm_MAG_N_bar.r
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
genome_taxonomy = load__genome_taxonomy(load__Stdb(), load__Wtdb())
genome.relative_abundance = get_relative_abundance(Wtdb_abd, genome_taxonomy)

gene_ko_tpm = load__gene_ko_tpm()

key_genes = load__key_genes()


### ######################################################################## ###
#### modify data for analysis                                               ####
### ######################################################################## ###
Stdb.spec =
  get_focus_species(genome.relative_abundance, Stdb, genome_taxonomy, 1) %>%
  # set order according to group distribution
  {.[.[unique(.$Group)] %>% apply(1, . %>% paste(collapse = "")) %>% order, ]}


tpm.key =
  gene_ko_tpm %>%
  merge(key_genes, by = "KO") %>%
  {.$Site = .$Layer %>% gsub("^(.+)\\.\\.(.+)$", "\\1", .); .} %>%
  {.$Pathway = factor(.$Pathway, levels = levels(key_genes$Pathway)); .} %>%
  merge(site_group) %>%
  merge(sample_meta[c("Location", "Group")] %>% unique)


### ######################################################################## ###
#### specifit modify to plot                                                ####
### ######################################################################## ###

tpm.key.N.sum.mean =
  tpm.key %>%
  {.[.$Pathway == "N", ]} %>%
  {.$Group = factor(.$Group, levels = names(sample_meta_col)); .} %>%
  split(paste(.$KO, .$Layer)) %>%
  lapply(function(df) {
    data.frame(
      df[1, c("Group", "Site", "KO", "Layer", "Label", "Arrow", "Pathway", "Location")],
      TPM = sum(df$TPM)
    )
  }) %>%
  bind_rows() %>%
  split(paste(.$Group, .$KO)) %>%
  lapply(function(df) {
    data.frame(
      df[1, c("Group", "KO", "Label", "Arrow", "Pathway", "Location")],
      TPM = mean(df$TPM)
    )
  }) %>%
  bind_rows() %>%
  .[order(.$Group, .$Arrow, -.$TPM), ]

p =
  tpm.key.N.sum.mean %>%
  split(.$Arrow) %>%
  lapply(. %>% data.frame(Width = sum(.$TPM))) %>%
  bind_rows() %>%
  {.$Width = .$Width / max(.$Width) * 0.95; .} %>%
  {
    ggplot(data = .,
           mapping = aes_string(x = "Arrow", y = "TPM",
                                color = "Group", fill = "Group")) +
      geom_bar(mapping = aes_string(width = "Width"),
               stat = "identity", alpha = 0.3,
               position = get(paste0("position_", "fill"))(reverse = FALSE)
               ) +  # stack
      geom_text(mapping = aes_string(label = "Label"),
                color = "black",
                position = get(paste0("position_", "fill"))(vjust = 0.5)) +
      scale_color_manual(values = sample_meta_col) +
      scale_fill_manual(values = sample_meta_col) +
      labs(y = "gene relative abundance") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }

ggsave(filename = fig_out, plot = p, width = 12, height = 8, limitsize = FALSE)
