###
#' @Date: 2022-07-06 13:55:31
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-09-10 12:10:21
#' @FilePath: /2021_09-MT10kSW/workflow/utils/RLib.local/R/drep.r
#' @Description:
###

load__Stdb <- function() {
  Stdb =
    file_path$file_path$results("Stdb") %>% as.character %>%
    read.csv %>% {colnames(.) = colnames(.) %>% str_to_title; .} %>%
    dplyr::select(!c("Genome_path")) %>%
    merge(site_group)
  return(Stdb)
}


load__Wtdb <- function() {
  Wtdb =
    file_path$file_path$results("Wtdb") %>% as.character %>%
    read.csv %>% {colnames(.) = colnames(.) %>% str_to_title; .} %>%
    {
      .$Taxonomy =
        {.} %>%
        {.[c("Classification", "Cluster")]} %>%
        {apply(., 1,
               function(x) {
                 x[1] %>%
                   gsub("(.)__;", str_glue("\\1__({x[2]});"), .) %>%
                   gsub("(.)__$", str_glue("\\1__({x[2]})"), .)
               })}
      .
    }
  return(Wtdb)
}


load__genome_taxonomy <- function(Stdb, Wtdb) {
  genome_taxonomy =
    merge(Wtdb[c("Cluster", "Taxonomy")], Stdb[c("Genome", "Secondary_cluster")],
          by.x = "Cluster", by.y = "Secondary_cluster") %>%
    {.$Rep = .$Genome %in% Wtdb$Genome; .} %>%
    {.$Genome = .$Genome %>% gsub(".fa$", "", .); .}
  return(genome_taxonomy)
}


get_relative_abundance <- function(genome_tdb_path, genome_taxonomy) {
  genome_tdb_path %>%
    read.csv(sep = "\t") %>%
    {.$Layer = paste0(.$Site, "..", .$Layer); .} %>%
    {.$Group = layer_2_group(.$Layer); .} %>%
    {.$Sample = paste0(.$Group, "_", .$Layer); .} %>%
    group_by(Layer) %>%
    mutate(Relative_abundance = Relative.abundance / sum(Relative.abundance) * 100) %>%
    #mutate(Relative_abundance = Relative.abundance) %>%
    merge(genome_taxonomy %>% .[.$Rep,])
}


assign_quality <- function(completeness, contamination) {
  {
    ifelse(
      (completeness < 50) | contamination >= 10, "low",
      ifelse(
        completeness <= 70, "mid",
        ifelse(
          contamination >= 5, "mid",
          ifelse(
            completeness <= 90, "good",
            ifelse(
              contamination >= 5, "good",
              "high"
            )
          )
        )
      )
    )
  } %>%
  factor(levels = c("low", "mid", "good", "high"))
}
