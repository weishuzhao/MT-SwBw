###
#' @Date: 2022-07-06 10:48:10
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-09-08 14:58:46
#' @FilePath: /2021_09-MT10kSW/workflow/others/check_relative_abundance.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


Wtdb_abd = stringr::str_glue("Wtdb.relative_abundance.tsv") %>% file_path$file_path$results() %>% as.character
Wtdb_abd_raw = stringr::str_glue("Stdb.raw.relative_abundance.tsv") %>% file_path$file_path$results() %>% as.character


### ######################################################################## ###
#### extract information of MAGs                                            ####
### ######################################################################## ###
genome_taxonomy = load__genome_taxonomy(load__Stdb(), load__Wtdb())
genome.relative_abundance = get_relative_abundance(Wtdb_abd, genome_taxonomy)
genome.raw.relative_abundance = get_relative_abundance(Wtdb_abd_raw, genome_taxonomy)

genome.relative_abundance.merge =
  merge(
    genome.relative_abundance, genome.raw.relative_abundance,
    by = c(
      "Cluster",
      "Site", "Layer", "Group", "Sample",
      "Taxonomy", "Genome", "Rep"
    ),
    all = TRUE,
    suffixes = c(".proper", ".raw")
  ) %>%
  {.[is.na(.)] = 0; .} %>%
  {.$more = .$Relative_abundance.raw - .$Relative_abundance.proper; .} %>%
  {.$more.pct = .$more / .$Relative_abundance.proper; .}

# unproper reads are mapped to many species in many samples
genome.relative_abundance.merge %>%
  .[.$Relative_abundance.proper == 0, ] %>%
  .$Taxonomy %>% unique %>% sort %>%
  {length(.) / length(genome.relative_abundance.merge$Taxonomy %>% unique %>% sort)}

# 1. cannot reject things with lower relative abundance
# 2. many MAGs have huge improvement
genome.relative_abundance.merge %>%
  split(.$Relative_abundance.proper == 0) %>%
  lapply(summary)

# Few records can pass max mapping check (4501 ~ 218)
genome.relative_abundance.merge %>%
  {
    fake.max =
      .[.$Relative_abundance.proper == 0,] %>%
      .$Relative_abundance.raw %>% max
    split(., .$Relative_abundance.raw > fake.max)
  } %>%
  lapply(summary)


# http://localhost:8787/graphics/plot_zoom_png?width=2000&height=4000
genome.relative_abundance.merge %>%
  {.$Taxon.level.spec = .$Taxonomy %>% taxon.split(1, 7); .} %>%
  group_by(
    Sample, Taxon.level.spec
  ) %>%
  summarise(
    proper = sum(Relative_abundance.proper),
    raw = sum(Relative_abundance.raw),
    more = sum(more)
  ) %>%
  {.[. == 0] = NA; .} %>%
  {
    ggplot(.) +
      geom_point(
        mapping = aes_string(
          x = "Sample", y = "Taxon.level.spec", size = "raw"),
        alpha = 1, color = "red"
      ) +
      geom_point(
        mapping = aes_string(
          x = "Sample", y = "Taxon.level.spec", size = "proper"),
        alpha = 0.5, color = "blue"
      ) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(hjust = 0)
      )
  }


genome.relative_abundance.merge %>%
  group_by(
    Sample, Cluster
  ) %>%
  summarise(
    proper = sum(Relative_abundance.proper),
    raw = sum(Relative_abundance.raw),
    more = sum(more)
  ) %>%
  {.[. == 0] = NA; .} %>%
  group_by(
    Sample
  ) %>%
  mutate(
    relative.raw = raw / sum(raw)
  ) %>%
  {
    ggplot(
      .,
      mapping = aes_string(
        x = "raw", y = "relative.raw", color = "is.na(proper)")
   ) +
      geom_point(
        shape = 20
      ) +
      scale_x_log10() + scale_y_log10() +
      stat_density2d()
  }
