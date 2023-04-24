###
#' @Date: 2022-06-25 10:52:06
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-08-11 11:24:49
#' @FilePath: /2021_09-MT10kSW/workflow/MAGs_tpm/relative_abundance_nmds.r
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
dist = argv[3]
method = argv[4]
taxon.level.spec = argv[4]

##### GLOBAL CONST vars                                                    #####

##### LOAD data AND transform TO basic format                              #####
genome_taxonomy = load__genome_taxonomy(
  load__Stdb(), load__Wtdb()
)
genome.relative_abundance = get_relative_abundance(Wtdb_abd, genome_taxonomy)


### ######################################################################## ###
#### Define function AND Calculate data                                     ####
### ######################################################################## ###
##### Define function                                                      #####

##### Calculate data                                                       #####
p_nmds_s =
  taxon.levels[-1] %>%
  as.list %>%
  {names(.) = unlist(.); .} %>%
  lapply(function(taxon.level.spec) {
    div.otu =
      genome.relative_abundance %>%
      {.$name = .$Taxonomy %>% taxon.split(1, taxon.level.spec); .} %>%
      {.$Abundance = get_relative_ce(., "Relative_abundance"); .} %>%
      reshape2::acast(formula = formula("name ~ Sample"),
                      value.var = "Abundance",
                      fun.aggregate = sum, fill = 0)
    list("pcoa" = "pcoa", "nmds" = "nmds") %>%
      lapply(function(method) {
        list("jaccard" = "jaccard", "bray" = "bray") %>%
          lapply(function(dist) {
            plot.beta.div(div.otu, pname = "relative abundance",
                          method = method, dist = dist, area = "polygon") +
              scale_color_manual(values = sample_meta_col) +
              scale_fill_manual(values = sample_meta_col)
          })
      })
  })


### ######################################################################## ###
#### Plot figures and OUTPUT                                                ####
### ######################################################################## ###
if (any(is.na(c(dist, method, taxon.level.spec)))) {
  p = NULL
  for (taxon.level.spec in taxon.levels[-1]) {
    for (dist in c("jaccard", "bray")) {
      for (method in c("pcoa", "nmds")) {
        if (is.null(p)) {
          p = p_nmds_s[[taxon.level.spec]][[method]][[dist]] +
            labs(x = "", y = taxon.level.spec)
        } else if ((dist == "jaccard") &
                   (method == "pcoa") &
                   (taxon.level.spec != taxon.levels[length(taxon.levels)])) {
          p = p + p_nmds_s[[taxon.level.spec]][[method]][[dist]] +
            labs(x = "", y = taxon.level.spec)
        } else if ((dist == "jaccard") &
                   (method == "pcoa") &
                   (taxon.level.spec == taxon.levels[length(taxon.levels)])) {
          p = p + p_nmds_s[[taxon.level.spec]][[method]][[dist]] +
            labs(x = paste(dist, method), y = taxon.level.spec)
        } else if ((taxon.level.spec == taxon.levels[length(taxon.levels)])) {
          p = p + p_nmds_s[[taxon.level.spec]][[method]][[dist]] +
            labs(x = paste(dist, method), y = "")
        } else {
          p = p + p_nmds_s[[taxon.level.spec]][[method]][[dist]] +
            labs(x = "", y = "")
        }
      }
    }
  }
  p = p +
    plot_layout(ncol = 4, guides = "collect") +
    plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")
  ggsave(fig_out, p, width = 22, height = 33)
} else {
  #dist = "bray"
  #taxon.level.spec = "class"
  p = p_nmds_s[[taxon.level.spec]][[method]][[dist]]
  ggsave(fig_out, p, width = 6, height = 6)
}
