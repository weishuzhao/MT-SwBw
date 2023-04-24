###
#' @Date: 2022-03-23 00:39:47
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-08-09 09:07:10
#' @FilePath: /2021_09-MT10kSW/workflow/others/draw_supp_fig1.r
#' @FromFilePath: /2021_09-MT10kSW/workflow/reads_diversity/rarefy_phyloFlash.r
#' @Description:
#'  Fig. S1. Annotation percentage of mOTU of metagenomes of seawater and sediment samples used in this study.
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, keyword_args, fig_out_path                         #####
fig_out = argv[1]

##### GLOBAL CONST vars                                                    #####

##### LOAD data AND transform TO basic format                              #####
div.phyloflash.raw =
  file_path$file_path$otus("phyloFlash_raw") %>% as.character %>%
  read.csv() %>%
  {
    .$Site = gsub(
      "^.+phyloFlash\\.\\.(.+)\\.\\.(.+)\\.csv$", "\\1..\\2", .$File
    )
    .
  } %>%
  reshape2::acast(formula("OTU ~ Site"), value.var = "Reads", fill = 0)


### ######################################################################## ###
#### Define function AND Calculate data                                     ####
### ######################################################################## ###
##### Define function                                                      #####

##### Calculate data                                                       #####
div.lastannot = sapply(
  as.character(taxon.levels),
  function(x)
    apply(div.phyloflash.raw[!taxon.split(rownames(div.phyloflash.raw), x) %>%
                               grepl(pattern = "^\\([^;]+\\)$", .),],  # > 0,#,
          2, sum))

div.lastannot.min = div.lastannot %>%
  (function(x) x / x[, 1] * 100)(.) %>%
  {
    Layer = rownames(.)
    data.frame(Min = apply(., 2, min),
               Layer = apply(., 2, . %>% {Layer[which.min(.)]}),
               Taxon = factor(colnames(.), levels = taxon.levels))
  } %>%
  {
    .$Taxon.last =
      .$Taxon %>% as.character %>% {c(NA, .)[1:length(.)]} %>%
      factor(levels = taxon.levels)
    .
  }

div.raw =
  div.lastannot %>%
  apply(., 1, function(x) x - c(x[-1], 0)) %>%
  reshape2::melt(varnames = c("Annotated level", "Layer"), value.name = "Reads",
                 as.is = TRUE) %>%
  as.data.frame() %>%
  {
    .[, "Annotated level"] = factor(.[, "Annotated level"],
                                  levels = rev(taxon.levels))
    .
  } %>%
  {.$Group = layer_2_group(.$Layer); .} %>%
  {.$Sample = paste0(.$Group, "_", .$Layer); .}

paletteer::paletteer_dynamic("cartography::pastel.pal", 20) %>%
  {
    pathways = taxon.levels
    cols = c(.[1:length(pathways)])
    names(cols) = pathways
    cols
  }

### ######################################################################## ###
#### Plot figures and OUTPUT                                                ####
### ######################################################################## ###
p = get_percent_plot(div.raw, "Reads", "Annotated level",
                      labs.x = "", labs.y = "relative abundance") +
  scale_x_discrete(
    limits =
      sample_meta[c("Site", "Layers")] %>%
      apply(1, . %>% paste(collapse = ".."))
  ) +
  theme(axis.text.x = element_text(color = sample_meta_col[sample_meta$Group]))

p2 =
  p +
  geom_hline(
    data = div.lastannot.min[2:6,],
    mapping = aes_string(yintercept = "Min",
                         linetype = "Taxon",
                         color = "Taxon"),
    size = 0.7
  ) +
  scale_color_manual(
    values =
      paletteer::paletteer_d("ggsci::default_igv")[2:7] %>% rev %>% c %>%
      {names(.) = taxon.levels[-1]; .}
  ) +
  guides(color = "none") +
  geom_point(
    data = div.lastannot.min[2:6,],
    mapping = aes_string(x = "Layer", y = "Min", fill = NULL),
    size = 5, shape = 10
  )


ggsave(filename = fig_out, plot = p2, width = 8, height = 5, dpi = 300)
