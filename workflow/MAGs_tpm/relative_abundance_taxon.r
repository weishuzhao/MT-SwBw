###
#' @Date: 2022-06-25 10:52:06
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-08-05 22:52:54
#' @FilePath: /2021_09-MT10kSW/workflow/MAGs_tpm/relative_abundance_taxon.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, fig_out_path, keyword_args                         #####
Wtdb_abd = argv[1]
#Wtdb_abd = stringr::str_glue("Stdb.relative_abundance.tsv") %>% file_path$file_path$results() %>% as.character
fig_out = argv[2]


##### GLOBAL CONST vars                                                    #####
font_size_1 = 15
font_size_2 = 13
font_size_3 = 10
axis.ticks.length = 0.1


##### LOAD data AND transform TO basic format                              #####
genome_taxonomy = load__genome_taxonomy(load__Stdb(), load__Wtdb())
genome.relative_abundance = get_relative_abundance(Wtdb_abd, genome_taxonomy)


### ######################################################################## ###
#### Define function AND Calculate data                                     ####
### ######################################################################## ###
##### Define function                                                      #####

##### Calculate data                                                       #####

### ######################################################################## ###
#### Plot figures and OUTPUT                                                ####
### ######################################################################## ###
##### Plot figures                                                         #####
p1 = get_total_plot(genome.relative_abundance, "Relative_abundance",
                    labs.x = "total mapping rate")

taxon_color = load__taxon_color()

p2 =
  genome.relative_abundance %>%
  {.$Taxonomy = sapply(taxon.split(.$Taxonomy, 1, 7), get_taxon_color); .} %>%
  get_percent_plot("Relative_abundance", "Taxonomy",
                   labs.x = "Layer", labs.y = "relative abundance") +
  scale_fill_manual(
    values =
      c(taxon_color$LEGEND_COLORS) %>%
      {names(.) = taxon_color$LEGEND_LABELS; .} %>%
      .[order(names(.))]
  )


##### OUTPUT                                                               #####
pout =
  p1 + p2 + guides(fill = "none") +
  plot_layout(heights = c(1, 4), nrow = 2, guides = "collect") +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 18))
ggsave(filename = fig_out, plot = pout, width = 14, height = 10, dpi = 300)
