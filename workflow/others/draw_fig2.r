###
#' @Date: 2022-07-20 13:43:25
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-11-16 16:45:35
#' @FilePath: /2021_09-MT10kSW/workflow/others/draw_fig2.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, fig_out_path, keyword_args                         #####
Wtdb_abd = argv[1]
#Wtdb_abd = stringr::str_glue("Wtdb.relative_abundance.tsv") %>% file_path$file_path$results() %>% as.character
fig_out = argv[2]

##### GLOBAL CONST vars                                                    #####
font_size_1 = 13
font_size_2 = 10
font_size_3 = 10
axis.ticks.length = 0.1

##### LOAD data AND transform TO basic format                              #####
genome_taxonomy = load__genome_taxonomy(load__Stdb(), load__Wtdb())
genome.relative_abundance = get_relative_abundance(Wtdb_abd, genome_taxonomy)

taxon_color = load__taxon_color()

### ######################################################################## ###
#### Define function AND Calculate data                                     ####
### ######################################################################## ###
##### Define function                                                      #####

##### Calculate data                                                       #####

### ######################################################################## ###
#### Plot figures and OUTPUT                                                ####
### ######################################################################## ###
##### Plot figures                                                         #####
set.seed(589)
p1s =
  {
    div.otu =
      genome.relative_abundance %>%
      {.$name = .$Taxonomy %>% taxon.split(1, 3); .} %>%
      reshape2::acast(formula = formula("name ~ Sample"),
                      value.var = "Relative_abundance",
                      fun.aggregate = sum, fill = 0)
    list("jaccard" = "jaccard", "bray" = "bray") %>%
      lapply(function(dist) {
        p =
          plot.beta.div(div.otu, pname = "relative abundance",
                        method = "nmds", dist = dist,
                        draw_labels = FALSE, area = "p") +
          scale_color_manual(values = sample_meta_col) +
          scale_fill_manual(values = sample_meta_col) +

          theme(
            axis.line = element_line(colour = "black"),
            axis.text = element_text(size = font_size_2, colour = "black", face = "bold"),
            axis.title = element_text(size = font_size_1, face = "bold", colour = "black"),
            axis.ticks.length = unit(axis.ticks.length, 'cm'),

            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
          ) +
          theme(
            legend.title = element_text(size = font_size_2, face = "bold"),
            legend.text = element_text(size = font_size_3),
            legend.key = element_blank()  # element_rect(fill = "gray")
            #legend.position = "bottom"
          )
        p$labels$title =
          p$labels$title %>%
          gsub("^.+plot of ", "", .) %>%
          gsub("\nrelative abundance.*$", "", .)
        p
      })
  }
# p1s$jaccard + p1s$bray

p2 =
  genome.relative_abundance %>%
  #{.$Taxa_label = taxon.split(.$Taxonomy, 1, 3); .} %>%
  {.$Taxa_label = sapply(taxon.split(.$Taxonomy, 1, 7), get_taxon_color); .} %>%
  get_percent_plot("Relative_abundance", "Taxa_label",
                   labs.x = "sample", labs.y = "relative abundance",
                   font_size_1 = font_size_1, font_size_2 = font_size_2, font_size_3 = font_size_3, axis.ticks.length = axis.ticks.length)
# p2$data[c("Group", "Site", "Layer", "Sample", "Abundance", "name", "Genome", "Taxonomy")] %>% write.csv("results/figs/fig2_relative_nmds_class.C.csv", row.names = FALSE, quote = FALSE)
p2.x = p2 +
  scale_fill_manual(
    values =
      c(taxon_color$LEGEND_COLORS) %>%
      {names(.) = taxon_color$Taxa_label; .} %>%
      .[order(names(.))]
  ) +
  scale_x_discrete(
    limits =
      sample_meta[c("Site", "Layers")] %>%
      apply(1, . %>% paste(collapse = ".."))
  ) +
  theme(axis.text.x = element_text(color = sample_meta_col[sample_meta$Group]))


##### OUTPUT                                                               #####
pout =
  p1s$jaccard + p1s$bray + p2.x +
  plot_layout(
    design = "AB
              CC",
    guides = 'collect'
  ) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 18))
ggsave(filename = fig_out, plot = pout, width = 13, height = 10)
