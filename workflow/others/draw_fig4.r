###
#' @Date: 2022-07-20 13:43:25
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-08-28 20:00:29
#' @FilePath: /2021_09-MT10kSW/workflow/others/draw_fig4.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, keyword_args, fig_out_path                         #####
fig_out = argv[1]

##### GLOBAL CONST vars                                                    #####
fun.aggregate = "sum"

##### LOAD data AND transform TO basic format                              #####
gene_ko_tpm = load__gene_ko_tpm()

key_genes = load__key_genes()
key_genes_N =
  key_genes %>% .[.$Pathway == "N",] %>%
  {.$Pathway = .$Arrow %>% factor(levels = unique(.)); .}
TOTAL_GENE_NUM = key_genes_N$KO %>% unique %>% length %>% {. * 2}

pathway_col = set__pathway_col(key_genes_N$Pathway)


set.seed(404)
site_group.sample =
  site_group %>%
  split(.$Group) %>%
  lapply(. %>% sample_n(3)) %>%
  bind_rows

### ######################################################################## ###
#### Define function AND Calculate data                                     ####
### ######################################################################## ###
##### Define function                                                      #####
report__tpm.key_heatmap_site <- function(tpm_key, fun.aggregate) {
  tpm_key %>%
    {.$Sample = str_glue("{.$Group}_{.$Layer}"); .} %>%
    reshape2::acast("Sample ~ KO", value.var = "TPM",
                    fun.aggregate = get(fun.aggregate)) %>%
    {
      sorted.ko = key_genes_N$KO %>% levels
      .[, sorted.ko[sorted.ko %in% colnames(.)]]
    } %>%
    {.[is.na(.)] = 0; .} %>%
    {
      annotation_col =
        key_genes_N %>%
        {
          data.frame(
            Pathway = .$Pathway, Label = .$Label, row.names = .$KO
          )
        }
      pheatmap::pheatmap(
        {d = .; d[d == 0] <- NA; d},
        color = colorRampPalette(c("seashell", "khaki", "orange", "red", "red3",
                                   "red4", "purple4", "purple3", "purple2", "purple1"))(100),
        #main = fun.aggregate,
        cellwidth = 10, cellheight = 10,
        cluster_cols = FALSE,
        cluster_rows = hclust(dist(., method = "euclidean"), "complete"),
        gaps_col = annotation_col["Pathway"] %>% table %>% .[. > 0] %>% cumsum,
        annotation_row =
          sample_meta %>%
          {data.frame(.[c("Group")], row.names = str_glue("{.$Group}_{.$Site}..{.$Layer}"))},
        annotation_col = annotation_col[c("Pathway", "Label")] %>% {.["Label"] = NULL; .},
        labels_col = annotation_col[, "Label"] %>% as.character,
        annotation_colors =
          list(
            "Group" = sample_meta_col,
            "Pathway" = pathway_col[annotation_col$Pathway %>% unique]
            #"diff" = c("others" = "#4DBBD5", "all.signif" = "#FFDC91", "mixed" = "#E64B35")
          ),
        silent = TRUE
      ) %>%
        ggplotify::as.ggplot(.)
    }
}


report__tpm.key_bar_group <- function(tpm_key) {
  tpm_key.signif =
    tpm_key %>%
    location.group.signif("KO", "TPM", "Location", "Group", TOTAL_GENE_NUM) %>%
    merge(key_genes_N, by = "KO") %>%
    .[.$p.value.char != "-",]

  tpm_key =
    tpm_key %>%
    {
      merge(
        merge(.[c("KO", "Label", "Pathway")] %>% unique,
              .[c("Group", "Location")] %>% unique),
        ., all.x = TRUE
      )
    } %>%
    {.[is.na(.$TPM), "Hide"] = "TRUE"; .[is.na(.)] = 1; .}

  #### add paint                                                            ####
  font_size_1 = 14
  font_size_2 = 11
  font_size_3 = 2.5
  axis.ticks.length = 0.1

  p =
    ggplot(data = tpm_key) +
    #geom_signif(mapping = aes_string(x = "Location", y = "TPM"),
    #            map_signif_level = TRUE,
    #            comparisons = list(c("ME", "MT"))) +
    geom_boxplot(mapping = aes_string(x = "Label", y = "TPM",
                                      fill = "Group",
                                      linetype = "Hide")) +
    scale_linetype_manual(values = c("1" = 1, "TRUE" = 0)) +
    guides(linetype = "none") +
    scale_fill_manual(values = sample_meta_col) +

    facet_grid(formula("Location ~ Pathway"),
               scales = "free_x", space = "free_x") +

    geom_label(data =
                 tpm_key.signif %>%
                 .[.$p.value.char != "", ],
               mapping = aes_string(x = "Label", y = "Max. + max(Max.) * 0.05",
                                    label = "p.value.char"),
               label.padding = unit(0.05, "lines"), label.size = 0,
               fill = "#7f7f7f3f", size = font_size_3) +
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
tpm.key =
  gene_ko_tpm %>%
  group_by(KO, Layer) %>%
  summarise(TPM = sum(TPM)) %>%
  merge(key_genes_N, by = "KO") %>%
  {.$Pathway = factor(.$Pathway, levels = levels(key_genes_N$Pathway)); .} %>%
  {.$Site = .$Layer %>% gsub("^(.+)\\.\\.(.+)$", "\\1", .); .} %>%
  merge(site_group) %>%
  merge(sample_meta[c("Location", "Group")] %>% unique) %>%
  {.$Location = factor(.$Location, levels = c("Slope", "Bottom")); .}

### ######################################################################## ###
#### Plot figures and OUTPUT                                                ####
### ######################################################################## ###
##### Plot figures                                                         #####
annotation_col =
  key_genes_N %>%
  {
    data.frame(
      Pathway = .$Pathway, Label = .$Label, row.names = .$KO
    )
  }

p1 =
  tpm.key %>%
  report__tpm.key_heatmap_site(fun.aggregate)
p2 =
  tpm.key %>%
#  .[.$Site %in% site_group.sample$Site, ] %>%
  {report__tpm.key_bar_group(.) + guides(fill = "none")}


##### OUTPUT                                                               #####
pout =
  p1 + p2 + plot_layout(design = c(
    patchwork::area(1, 1, 5, 8),
    patchwork::area(6, 1, 9, 7)
  )) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 18))
ggsave(filename = fig_out, plot = pout, width = 14, height = 12)
