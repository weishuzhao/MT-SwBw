###
#' @Date: 2022-03-23 00:39:47
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-08-05 14:16:34
#' @FilePath: /2021_09-MT10kSW/workflow/reads_diversity/rarefy_phyloFlash.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, keyword_args, fig_out_path                         #####
fig_out = argv[1]
tab_out = argv[2]
#tab_out = as.character(file_path$file_path$otus("phyloFlash_percent"))

##### GLOBAL CONST vars                                                    #####
font_size_1 = 15
font_size_2 = 13
font_size_3 = 10
axis.ticks.length = 0.1

TOP_N_TAXON_PER_LAYER = 7

##### LOAD data AND transform TO basic format                              #####
div.phyloflash.raw =
  file_path$file_path$otus("phyloFlash_raw") %>% as.character %>%
  read.csv() %>%
  {.$Site = gsub("^.+phyloFlash\\.\\.(.+)\\.\\.(.+)\\.csv$", "\\1..\\2", .$File); .} %>%
  reshape2::acast(formula("OTU ~ Site"), value.var = "Reads", fill = 0)


### ######################################################################## ###
#### Define function AND Calculate data                                     ####
### ######################################################################## ###
##### Define function                                                      #####
div.lastannot = sapply(
  as.character(taxon.levels),
  function(x)
    apply(div.phyloflash.raw[!taxon.split(rownames(div.phyloflash.raw), x) %>%
                               grepl(pattern = "^\\([^;]+\\)$", .),],  # > 0,#,
          2, sum))

##### Calculate data                                                       #####
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

p1 = get_total_plot(div.raw, "Reads", labs.x = "total annot reads") +
  scale_x_discrete(
    limits =
      sample_meta[c("Site", "Layers")] %>%
      apply(1, . %>% paste(collapse = "..")),
    labels = NULL
  )
p2 = get_percent_plot(div.raw, "Reads", "Annotated level",
                      labs.x = "Layer", labs.y = "relative abundance") +
  scale_x_discrete(
    limits =
      sample_meta[c("Site", "Layers")] %>%
      apply(1, . %>% paste(collapse = ".."))
  ) +
  theme(axis.text.x = element_text(color = sample_meta_col[sample_meta$Group]))



pout =
  p1 + p2 + plot_layout(heights = c(1, 4), nrow = 2, guides = "collect")

ggsave(filename = fig_out, plot = pout, width = 8, height = 6, dpi = 300)


### ######################################################################## ###
#### generate div.phyloflash.rarefy FOR all following study                 ####
### ######################################################################## ###
div.rarefy.save <- function(div.phyloflash.raw,
                            min_rarefy_threshold = 50000,
                            file = NULL) {
  set.seed(499)
  div.phyloflash.rarefy =
    div.phyloflash.raw %>%
    {t(.)} %>%
    {.[apply(., 1, sum) > min_rarefy_threshold, ]} %>%
    {vegan::rrarefy(., apply(., 1, sum) %>% min)} %>%
    {t(.)}
  if (!is.null(file)) {
    write.csv(div.phyloflash.rarefy, file = file)
  }
  return(invisible(div.phyloflash.rarefy))
}

div.percent.save <- function(div.phyloflash.raw, file = NULL) {
  div.phyloflash.percent =
    div.phyloflash.raw %>%
    {t(.) / apply(., 2, sum) * 100}  %>%
    {t(.)}
  if (!is.null(file)) {
    write.csv(div.phyloflash.percent, file = file)
  }
  return(invisible(div.phyloflash.percent))
}
div.percent.save(div.phyloflash.raw, tab_out)
