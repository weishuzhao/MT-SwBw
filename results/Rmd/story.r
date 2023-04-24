###
#' @Date: 2022-07-05 15:11:27
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-08-22 09:56:05
#' @FilePath: /2021_09-MT10kSW/workflow/story.r
#' @Description:
#'      Use R to prove all my points
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


module.name = get_module.name()

Wtdb_abd = file_path$file_path$results("Stdb.relative_abundance.tsv") %>% as.character
genome_taxonomy = load__genome_taxonomy(load__Stdb(), load__Wtdb())
genome.relative_abundance = get_relative_abundance(Wtdb_abd, genome_taxonomy)
gmodule = file_path$file_path$results("gmodule.csv") %>% as.character %>% load_gmodule(symnum = FALSE)
gene_ko_tpm = load__gene_ko_tpm()

## 1.  Slope water have only a few MAGs
genome.relative_abundance %>%
  {.$name = taxon.split(.$Taxonomy, 1, 7); .} %>%
  {.[c("group", "name")]} %>%
  {table(.$group)} %>%
  {message.print(.); .} %>%
  {assertthat::assert_that(all(. == c(1178, 1275, 1323, 248)),
                           msg = "Check if Ss is not fewer than other samples")}

## 2.  MAG number recovered from single sample is not significant lower than other samples
genome.relative_abundance %>%
  {.$name = taxon.split(.$Taxonomy, 1, 7); .} %>%
  {table(.$sample)} %>%
  data.frame %>%
  {.$group = .$Var1 %>% gsub("^([^_]+)_(.+)$", "\\1", .); .} %>%
  group.signif.mark("Freq", "group") %>%
  {
    assertthat::assert_that(.["Sw", "mean"] > .["Bs", "mean"])
    assertthat::assert_that(.["Sw", "mean"] > .["Ss", "mean"])
    .
  }


TOTAL_GENE_NUM = gene_ko_tpm$KO %>% unique %>% length %>% {. * 2}
check_cross_ko_signif <- function(ko_tpm, site_group) {
  ko_tpm %>%
    {.$Site = .$Layer %>% gsub("^(.+)\\.\\.(.+)$", "\\1", .); .} %>%
    merge(site_group) %>%
    merge(sample_meta[c("location", "group")] %>% unique) %>%
    location.group.signif("KO", "tpm", "location", "group", TOTAL_GENE_NUM) %>%
    split(.$location) %>%
    {merge(.$Bottom, .$Slope, by = "KO", suffixes = c(".Bottom", ".Slope"), all = TRUE)} %>%
    reshape2::acast(formula("p.value.char.Bottom ~ p.value.char.Slope"),
                    fun.aggregate = length) %>%
    .[c(7, 2, 1, 3, 4, 5, 6), c(7, 2, 1, 3, 4, 5, 6)]
}
set.seed(404)
site_group.sample =
  site_group %>%
  split(.$group) %>%
  lapply(. %>% sample_n(3)) %>%
  bind_rows
site_group.sample %>%
  check_cross_ko_signif(gene_ko_tpm, .) %>%
  {
    assertthat::assert_that(.[7, 3] / sum(.) * (7 ** 2) > 5, msg = "many kos are only signif difference in bottom but no difference in slope")
    assertthat::assert_that(.[7, 7] / sum(.) * (7 ** 2) > 5, msg = "many kos are very signif difference in both location")
    .
  }


key_genes = load__key_genes()


Stdb = load__Stdb()

checkm = read.csv("results/checkm.tsv", sep = "\t")
checkm %>%
  {.$Site = .$Bin.Id %>% gsub("-(concoct|metabat|maxbin).+", "", .); .} %>%
  .[.$Completeness > 70 & .$Contamination > 10, ] %>%
  .[c("Bin.Id", "Site", "Marker.lineage", "Completeness", "Contamination")] %>%
  .$Site %>% table %>%
  {
    bind_rows(., Stdb$Genome %>% gsub("-(concoct|metabat|maxbin).+", "", .) %>% table)
  } %>%
  .[, !is.na(apply(., 2, sum))] %>%
  t %>%
  data.frame(sum = apply(., 1, sum))
Stdb$Genome %>% gsub("-(concoct|metabat|maxbin).+", "", .) %>% table
