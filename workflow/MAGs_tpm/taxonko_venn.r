###
#' @Date: 2022-03-03 15:56:58
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-08-06 09:47:39
#' @FilePath: /2021_09-MT10kSW/workflow/MAGs_tpm/taxonko_venn.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, fig_out_path, keyword_args                         #####
Wtdb_abd = argv[1]
#Wtdb_abd = stringr::str_glue("Stdb.relative_abundance.tsv") %>% file_path$file_path$results() %>% as.character
fig_mag_venn = argv[2]

##### LOAD data AND transform TO basic format                              #####
genome_taxonomy = load__genome_taxonomy(load__Stdb(), load__Wtdb())
genome.relative_abundance = get_relative_abundance(Wtdb_abd, genome_taxonomy)

site_ko =
  load__gene_ko_tpm() %>%
  reshape2::acast(formula = formula("KO ~ Layer"),
                  fun.aggregate = sum, value.var = "TPM", fill = 0)
site_module =
  file_path$file_path$results("site_module.csv") %>% as.character %>%
  load_gmodule

### ######################################################################## ###
#### Define function AND Calculate data                                     ####
### ######################################################################## ###
##### Define function                                                      #####
#' DRAW Venn plot ACCORDING TO div.otu
#'
#' taxon that unannotated at given level will be regarded different across envirnment
#' e.g. if the taxon is "Bacteria;Proteobaceria;" after split,
#'      it will result in "Bacteria;Proteobaceria;env1" as it is obtained in env1
venn.grid.taxon = function(taxon.level.spec) {
  genome.relative_abundance %>%
    {.$name = .$Taxonomy %>% taxon.split(1, taxon.level.spec); .} %>%
    #.[!.$name %>% grepl(";\\([^;]+\\)", .),] %>%
    {.$Site = .$Layer %>% gsub("^(.+)\\.\\.(.+)$", "\\1", .); .} %>%
    merge(site_group) %>%
    {split(.$name, .$Group)} %>%
    lapply(unique) %>%
    {
      VennDiagram::venn.diagram(
        x = .,
        filename = NULL, imagetype = "png",
        fill = sample_meta_col[names(.)], alpha = 0.75,
        lwd = 3,
        label.col = "black",
        cex = 2, fontfamily = "Arial", fontface = "bold",
        main = taxon.level.spec, main.cex = 2, main.fontfamily = "Arial",
        main.pos = c(0.5, 0), main.just = c(0.5, 1),
        #cat.col = {names(.) %>% },
        ext.line.lty = "dotted", ext.dist = -0.1,
        disable.logging = TRUE
      )
    } %>%
    {ggpubr::as_ggplot(.)} %>%
    {. + theme(plot.margin = unit(rep(0.3, 4), "in"))}
}


#' DRAW Venn plot ACCORDING TO KO OR module
venn.grid.gene = function(x, pname) {
  x %>%
    reshape2::melt(varnames = c("KO", "Layer"), value.name = "TPM") %>%
    {.[.$TPM > 0,]} %>%
    {.$Site = .$Layer %>% gsub("^(.+)\\.\\.(.+)$", "\\1", .); .} %>%
    merge(site_group) %>%
    {split(.$KO, .$Group)} %>%
    lapply(unique) %>%
    {
      VennDiagram::venn.diagram(
        x = .,
        filename = NULL, imagetype = "png",
        fill = sample_meta_col[names(.)], alpha = 0.75,
        lwd = 3,
        label.col = "black",
        cex = 1.5, fontfamily = "Arial", fontface = "bold",
        main = pname, main.cex = 2, main.fontfamily = "Arial",
        main.pos = c(0.5, 0), main.just = c(0.5, 1),
        #cat.col = c(ME = NA, MT = NA),
        ext.line.lty = "dotted", ext.dist = -0.1,
        disable.logging = TRUE
      )
    } %>%
    {ggpubr::as_ggplot(.)} %>%
    {. + theme(plot.margin = unit(rep(0.3, 4), "in"))}
}


##### Calculate data                                                       #####

### ######################################################################## ###
#### Plot figures and OUTPUT                                                ####
### ######################################################################## ###
##### Plot figures                                                         #####
p_all <- NULL
for (taxon.level.spec in taxon.levels) {
  if (is.null(p_all)) {
    p_all = venn.grid.taxon(taxon.level.spec)
  } else {
    p_all = p_all + venn.grid.taxon(taxon.level.spec)
  }
}

##### OUTPUT                                                               #####
p =
  p_all +
  venn.grid.gene(site_ko, "KO") +
  venn.grid.gene(site_module > 0, "module") +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")
ggsave(filename = fig_mag_venn,
       plot = p,
       width = 10, height = 10, dpi = 300)
