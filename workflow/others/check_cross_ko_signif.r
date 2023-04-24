###
#' @Date: 2022-07-19 11:21:27
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-07-19 19:44:01
#' @FilePath: /2021_09-MT10kSW/workflow/others/check_cross_ko_signif.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


key_genes = load__key_genes()

gene_ko_tpm = load__gene_ko_tpm()


#TOTAL_GENE_NUM = gene_ko_tpm$KO %>% unique %>% length %>% {. * 2}
TOTAL_GENE_NUM = key_genes$KO %>% length %>% {. * 2}


site_group %>%
  {.[c("location", "property")] = .$group %>% str_split_fixed("", 2); .} %>%
  {.$location = ifelse(.$location == "B", "Bottom", "Slope"); .} %>%
  {.$property = ifelse(.$property == "w", "water", "sediment"); .}

cross_ko_signif =
  gene_ko_tpm %>%
  .[.$ko %in% key_genes$KO, ] %>%
  {aggregate(.$tpm, by = .[c("ko", "Layer")], FUN = sum)} %>% {.$tpm = .$x; .$x <- NULL; .} %>%
  {.$Site = .$Layer %>% gsub("^(.+)\\.\\.(.+)$", "\\1", .); .} %>%
  merge(
    site_group %>%
      {.[c("location", "property")] = .$group %>% str_split_fixed("", 2); .} %>%
      {.$location = ifelse(.$location == "B", "Bottom", "Slope"); .} %>%
      {.$property = ifelse(.$property == "w", "water", "sediment"); .}
  ) %>%
  location.group.signif("ko", "tpm", "location", "group", TOTAL_GENE_NUM)

cross_ko_signif.property =
  gene_ko_tpm %>%
  .[.$ko %in% key_genes$KO, ] %>%
  {aggregate(.$tpm, by = .[c("ko", "Layer")], FUN = sum)} %>% {.$tpm = .$x; .$x <- NULL; .} %>%
  {.$Site = .$Layer %>% gsub("^(.+)\\.\\.(.+)$", "\\1", .); .} %>%
  merge(
    site_group %>%
      {.[c("location", "property")] = .$group %>% str_split_fixed("", 2); .} %>%
      {.$location = ifelse(.$location == "B", "Bottom", "Slope"); .} %>%
      {.$property = ifelse(.$property == "w", "water", "sediment"); .}
  ) %>%
  reshape2::dcast(formula = formula("ko + location ~ property"),
                  value.var = "tpm", fun.aggregate = mean, fill = 0) %>%
  {.$more_in_sediment = ifelse(.$sediment < .$water, -1, 1); .} %>%
  merge(cross_ko_signif) %>%
  {.}

cross_ko_signif.property %>%
  .[.$location == "Slope", ] %>%
  .[order(.$p.value.adj), ]
cross_ko_signif.property %>%
  {
    merge(.[.$location == "Bottom", ], .[.$location == "Slope", ],
          by = "ko", suffixes = c(".Bottom", ".Slope"), all = TRUE)
  } %>%
  reshape2::acast(formula("p.value.char.Bottom ~ p.value.char.Slope"),
                  fun.aggregate = length) %>%
  .[c(7, 2, 1, 3, 4, 5, 6), c(7, 2, 1, 3, 4, 5, 6)]

cross_ko_signif.property %>%
  .[.$p.value.adj != 2, ] %>%
  {
    merge(.[.$location == "Bottom", ], .[.$location == "Slope", ],
          by = "ko", suffixes = c(".Bottom", ".Slope"), all = TRUE)
  } %>%
  merge(key_genes, by.x = "ko", by.y = "KO", all.x = TRUE) %>%
  {
    ggplot(data = .,
           mapping = aes_string(
             x = "(1 - p.value.adj.Bottom) * more_in_sediment.Bottom",
             y = "(1 - p.value.adj.Slope) * more_in_sediment.Slope"
           )) +
      geom_point(mapping = aes_string(fill = "pathway", shape = "!is.na(pathway)", size = "!is.na(pathway)"),
                 alpha = 0.5) +
      scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 1)) +
      geom_label_repel(mapping = aes_string(label = "label")) +
      geom_density2d(aes_string(color = "..level.."), h = c(0.2, 0.2),
                     alpha = 0.7)
  } %>%
  {.}
