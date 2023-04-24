###
#' @Date: 2022-05-04 09:42:54
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-09-01 22:37:59
#' @FilePath: /2021_09-MT10kSW/workflow/others/draw_supp_fig2.r
#' @Description:
#'  类群taxa相对丰度在Ss、Sw、Bs、Bw中的boxplot比较图，标sig
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, keyword_args, fig_out_path                         #####
Wtdb_abd = argv[1]
#Wtdb_abd = stringr::str_glue("Stdb.relative_abundance.tsv") %>% file_path$file_path$results() %>% as.character
fig_out = argv[2]

##### GLOBAL CONST vars                                                    #####

##### LOAD data AND transform TO basic format                              #####
Stdb = load__Stdb()
Wtdb = load__Wtdb()
genome_taxonomy = load__genome_taxonomy(Stdb, Wtdb)
genome.relative_abundance = get_relative_abundance(Wtdb_abd, genome_taxonomy)

taxon_color = load__taxon_color()

key_genes = load__key_genes()
key_genes_N =
  key_genes %>% .[.$Pathway == "N",] %>%
  {.$Pathway = .$Arrow %>% factor(levels = unique(.)); .}

genomeko = load__genomeko()
key_genes_N.genome =
  genomeko %>%
  reshape2::melt(varnames = c("KO", "Genome")) %>%
  .[.$value > 0,] %>%
  merge(key_genes_N) %>%
  {.$Genome = .$Genome %>% gsub(".fa$", "", .); .} %>%
  merge(genome_taxonomy[c("Cluster", "Genome")]) %>%
  {.$Genome <- NULL; .} %>%
  merge(genome_taxonomy[c("Cluster", "Genome", "Rep")] %>% .[.$Rep, ])

### ######################################################################## ###
#### Define function AND Calculate data                                     ####
### ######################################################################## ###
##### Define function                                                      #####
report__tpm.key_bar_group <- function(
  relative_abundance.class,
  TOTAL_TAXA_NUM = NA, df_x = NA, df_fill = NA,
  x = "Taxa_label", y = "Relative_abundance", fill = "Group",
  grid_formula = formula("Location ~ Phylum"),
  p.value.char.filter = . %>% .[!grepl("[-?]", .$p.value.char),]
) {
  if (is.na(TOTAL_TAXA_NUM)) {
    TOTAL_TAXA_NUM =
      relative_abundance.class$Taxa_label %>% unique %>% length %>% {. * 2}
  }
  if (is.na(df_x)) {
    df_x =
      relative_abundance.class %>%
      .[c(all.vars(grid_formula)[2], x)] %>%
      unique
  }
  if (is.na(df_fill)) {
    df_fill =
      relative_abundance.class %>%
      .[c(fill, all.vars(grid_formula)[1])] %>%
      unique
  }

  relative_abundance.class.signif =
    relative_abundance.class %>%
    location.group.signif(x, y, all.vars(grid_formula)[1], fill, TOTAL_TAXA_NUM) %>%
    merge(df_x) %>%
    p.value.char.filter

  relative_abundance.class =
    relative_abundance.class %>%
    {merge(merge(df_x, df_fill), ., all.x = TRUE)} %>%
    {.[is.na(.[y]), "Hide"] = "TRUE"; .[is.na(.)] = 1; .}

  #### add paint                                                            ####
  font_size_1 = 16
  font_size_2 = 14
  axis.ticks.length = 0.1

  p =
    ggplot(data = relative_abundance.class) +
    geom_boxplot(mapping = aes_string(x = x, y = y, fill = fill,
                                      linetype = "Hide")) +
    scale_linetype_manual(values = c("1" = 1, "TRUE" = 0)) +
    guides(linetype = "none") +
    scale_fill_manual(values = sample_meta_col) +

    facet_grid(grid_formula,
               scales = "free_x", space = "free_x") +

    geom_label(data =
                 relative_abundance.class.signif %>%
                 .[.$p.value.char != "", ],
               mapping = aes_string(x = x, y = "Max. + max(Max.) * 0.05",
                                    label = "p.value.char"),
               label.padding = unit(0.05, "lines"), label.size = 0,
               fill = "#7f7f7f3f") +
    labs(x = "")

  p1 =
    p +
    scale_y_log10(labels =
                    ~format(.x, scientific = TRUE) %>%
                    str_replace("^0e\\+0", "0e+0") %>%
                    str_replace("e\\+0", "%*%10^") %>%
                    parse(text = .)) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = 'white', size = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    ) +
    theme(
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = font_size_2, colour = "black", face = "bold"),
      axis.title = element_text(size = font_size_1, face = "bold", colour = "black"),
      axis.ticks.length = unit(axis.ticks.length, 'cm'),

      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "plain")
    ) +
    theme(
      legend.title = element_text(size = font_size_2, face = "bold"),
      legend.text = element_text(size = font_size_2, face = "bold")
      #legend.position = "bottom"
    ) +
    theme(text = element_text(family = "Arial",
                              size = font_size_1,
                              hjust = 0.5,
                              lineheight = 0.5)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.key = element_blank()) +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.x = element_blank())

  p1
}


##### Calculate data                                                       #####
relative_abundance.class =
  genome.relative_abundance %>%
  {.$Taxa_label = sapply(taxon.split(.$Taxonomy, 1, 7), get_taxon_color); .} %>%
  group_by(Taxa_label, Layer) %>% summarise(Relative_abundance = sum(Relative_abundance)) %>%
  {.$Phylum = .$Taxa_label %>% taxon.split(1, 2); .} %>%
  {.$Site = .$Layer %>% gsub("^(.+)\\.\\.(.+)$", "\\1", .); .} %>%
  merge(site_group) %>%
  merge(sample_meta[c("Location", "Group")] %>% unique) %>%
  {.$Location = factor(.$Location, levels = c("Slope", "Bottom")); .}

TOTAL_TAXA_NUM =
  relative_abundance.class$Taxa_label %>% unique %>% length %>% {. * 2}


### ######################################################################## ###
#### Plot figures and OUTPUT                                                ####
### ######################################################################## ###
p1 =
  relative_abundance.class %>%
  as.data.frame() %>%
  {report__tpm.key_bar_group(.)}

ggsave(filename = fig_out, plot = p1, width = 18, height = 9, dpi = 300)
