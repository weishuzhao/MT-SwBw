###
#' @Date: 2022-06-25 10:52:06
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-07-20 10:51:46
#' @FilePath: /2021_09-MT10kSW/workflow/MAGs_tpm/relative_abundance_taxon.1.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


Wtdb_abd = argv[1]
#Wtdb_abd = stringr::str_glue("Stdb.relative_abundance.tsv") %>% file_path$file_path$results() %>% as.character
fig_out = argv[2]


#### GLOBAL vars                                                            ####
font_size_1 = 15
font_size_2 = 13
font_size_3 = 10
axis.ticks.length = 0.1


#### LOAD table                                                             ####
genome_taxonomy = load__genome_taxonomy(load__Stdb(), load__Wtdb())
genome.relative_abundance = get_relative_abundance(Wtdb_abd, genome_taxonomy)


taxon_color = load__taxon_color()
get_taxon_color = function(taxonomy) {
  apply(
    taxon_color, 1, . %>% {ifelse(grepl(
      pattern = .[2],
      x = taxonomy
    ), .[2], "")}
  ) %>% {.[. != ""][1]} %>% {ifelse(is.na(.), "others", .)}
}


Stdb.spec =
  get_focus_species(genome.relative_abundance, load__Stdb(), genome_taxonomy, 1)

relative_abundance.max =
  genome.relative_abundance %>%
  split(.$cluster) %>% lapply(. %>% .[which.max(.$relative_abundance), ]) %>% bind_rows %>%
  .[c("cluster", "relative_abundance")] %>%
  merge(Stdb.spec, by = "cluster") %>%
  {.$taxon_color = sapply(.$name, get_taxon_color); .} %>%
  {.$group_set = .[unique(.$group)] %>% apply(1, . %>% paste(collapse = ".")); .}


p =
  ggplot(data = relative_abundance.max %>% {.$y = .[, "relative_abundance"]; .},  # totalAvgDepth
       mapping = aes_string(x = "group_set", y = "rank(-y)")) +
  geom_point(mapping = aes_string(size = "y", fill = "taxon_color"),
             shape = 21) +
  geom_text(mapping = aes_string(label = "cluster")) +
  scale_fill_manual(
    values =
      c(taxon_color$LEGEND_COLORS) %>%
      {names(.) = taxon_color$LEGEND_LABELS; .} %>%
      .[order(names(.))]
  ) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, face = "plain")
  )

ggsave("results/cache/genome_cross_abundance.svg", p, width = 10, height = 60, limitsize = FALSE)
