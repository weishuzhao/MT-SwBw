###
#' @Date: 2022-07-20 13:43:25
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2024-05-29 15:39:26
#' @FilePath: /2021_09-MT10kSW/workflow/others/draw_fig2.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, fig_out_path, keyword_args                         #####
Wtdb_abd <- stringr::str_glue("Wtdb.relative_abundance.tsv") %>% # nolint: object_name_linter, line_length_linter.
  file_path$file_path$results() %>%
  as.character()
Wtdb_abd <- argv[1] # nolint: object_name_linter.
fig_out <- argv[2]

##### GLOBAL CONST vars                                                    #####
font_size_1 <- 13
font_size_2 <- 10
font_size_3 <- 10
axis_ticks_length <- 0.1

##### LOAD data AND transform TO basic format                              #####
genome_taxonomy <- load__genome_taxonomy(load__Stdb(), load__Wtdb())
genome_rltabd <- get_relative_abundance(Wtdb_abd, genome_taxonomy)

sample_meta_cross <- read.csv(
  "results/reads_diversity/metadata.tsv",
  sep = "\t"
) %>%
  mutate(
    X = get("Sample"),
    Layer = ifelse(
      get("Type") == "16s", get("Sample"), gsub("^[^_]+_", "", get("Sample"))
    ),
    Sample = paste0(get("Group"), "_", get("Layer"))
  ) %>%
  merge(unique(sample_meta[c("Location", "Group")]))
cross_rltabd <- read.csv("results/reads_diversity/abundance.csv") %>%
  merge(sample_meta_cross[c("X", "Sample")]) %>%
  mutate(X = NULL) %>%
  column_to_rownames("Sample")


taxon_color <- load__taxon_color()

### ######################################################################## ###
#### Define function AND Calculate data                                     ####
### ######################################################################## ###
##### Define function                                                      #####
plot.beta.div <- function(dist = "jaccard", binary = NA) {
  # >>->> argParse
  if (is.na(binary)) {
    binary <- tryCatch(
      expr = {
        match.arg(dist, c("jaccard")) == "jaccard"
      },
      error = function(e) FALSE
    )
  }
  # <<-<<                                                               <<-<<

  title_test_method <- paste0(
    ifelse(binary, "binary ", ""), dist, " distance"
  )

  # >>->> adonis2
  group <- data.frame(Location = sapply(
    rownames(cross_rltabd),
    function(x) {
      unlist(strsplit(x, "\\_"))[1]
    }
  ))
  group_adonis2 <- vegan::adonis2(formula("cross_rltabd ~ Location"), group,
    method = dist, binary = binary, by = "margin"
  )
  title_adonis_sgnf <- paste0(
    "PERMANOVA",
    " R^2=", round(group_adonis2$R2[1], 4),
    " p(Pr(>F))=", group_adonis2$`Pr(>F)`[1]
  )
  # <<-<<                                                               <<-<<

  # >>->> Dimensionality reduction
  nmds_dis <- vegan::metaMDS(
    vegan::vegdist(cross_rltabd, method = dist, binary = binary),
    trace = 0
  )
  title_test_method <- paste0(
    title_test_method, ", stress=", as.character(round(nmds_dis$stress, 4))
  )
  if (nmds_dis$stress >= 0.2) {
    warning("应力函数值 >= 0.2, 不合理")
  }
  div_otu_point_ <- data.frame(nmds_dis$points)
  xylab <- paste0("NMDS ", 1:2)
  # <<-<<                                                               <<-<<

  div_otu_point <- div_otu_point_ %>% # nolint: object_usage_linter.
    `names<-`(c("Axis.1", "Axis.2")) %>%
    rownames_to_column("Sample") %>% # nolint: object_usage_linter.
    merge(sample_meta_cross)

  p <- ggplot(data = div_otu_point) + # nolint: object_usage_linter.
    geom_point( # nolint: object_usage_linter.
      aes_string( # nolint: object_usage_linter.
        x = "Axis.1", y = "Axis.2",
        color = "Group", shape = "Type"
      ),
      size = 2, alpha = 0.65
    ) +
    scale_color_manual(values = sample_meta_col) + # nolint
    scale_fill_manual(values = sample_meta_col) + # nolint
    labs( # nolint: object_usage_linter.
      title = paste(
        title_test_method, title_adonis_sgnf,
        sep = "\n"
      ),
      x = xylab[1], y = xylab[2]
    ) +
    geom_polygon( # nolint: object_usage_linter.
      data = Reduce(
        rbind,
        lapply(
          split(div_otu_point, div_otu_point$Group),
          function(x) x[chull(x[c("Axis.1", "Axis.2")]), ]
        )
      ),
      aes_string(
        x = "Axis.1", y = "Axis.2",
        fill = "Group", color = "Group"
      ),
      alpha = 0.15, linetype = 3
    )
  # <<-<<                                                               <<-<<

  return(p)
}

##### Calculate data                                                       #####

### ######################################################################## ###
#### Plot figures and OUTPUT                                                ####
### ######################################################################## ###
##### Plot figures                                                         #####
set.seed(589)
p1s <- {
  list("jaccard" = "jaccard", "bray" = "bray") %>%
    lapply(function(dist) {
      p <- plot.beta.div(dist = dist) +
        theme_bw() +
        theme(text = element_text(family = "Arial")) +
        theme(
          axis.text = element_text(
            size = font_size_2, colour = "black", face = "bold"
          ),
          axis.title = element_text(
            size = font_size_1, face = "bold", colour = "black"
          ),
          axis.ticks.length = unit(axis_ticks_length, "cm"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
        ) +
        theme(
          legend.title = element_text(size = font_size_2, face = "bold"),
          legend.text = element_text(size = font_size_3),
          legend.key = element_blank()
        )
      p
    })
}
list("jaccard" = "jaccard", "bray" = "bray") %>%
  lapply(
    function(dist) {
      g2 <- sample_meta_cross %>%
        mutate(
          water_sediment = gsub("^[BS]", "", .data$Group)
        ) %>%
        .[c("Sample", "Group", "Location", "water_sediment", "Type")] %>%
        column_to_rownames("Sample") %>%
        .[rownames(cross_rltabd), ]
      g1 <- data.frame(Location = sapply(
        rownames(cross_rltabd),
        function(x) {
          unlist(strsplit(x, "\\_"))[1]
        }
      ))
      group_adonis2 <- vegan::adonis2(
        formula("cross_rltabd ~ Group * Type"),
        g2,
        method = dist, binary = dist == "jaccard"
      )
      paste0(
        "PERMANOVA",
        " R^2=", round(group_adonis2$R2[1], 4),
        " p(Pr(>F))=", group_adonis2$`Pr(>F)`[1]
      ) %>% print()
      group_adonis2
    }
  )
#' p1s$jaccard + p1s$bray

sample_meta_1 <- sample_meta %>%
  arrange(-rank(.data$Group), .data$Depth) %>%
  mutate(
    x = c(1:3, 4:5, 7:10, 1:7, 3:14, 16:20),
    Layer = paste0(.data$Site, "..", .data$Layers)
  ) %>%
  .[c("Layer", "Group", "x")]

p2 <- sample_meta_1 %>%
  inner_join(genome_rltabd) %>%
  mutate(
    Layer = paste0(.data$Group, .data$x),
    Taxa_label = sapply(taxon.split(.data$Taxonomy, 1, 7), get_taxon_color)
  ) %>%
  get_percent_plot("Relative_abundance", "Taxa_label",
    sample.name = "Layer",
    labs.x = "sample", labs.y = "relative abundance",
    font_size_1 = font_size_1, font_size_2 = font_size_2,
    font_size_3 = font_size_3,
    axis.ticks.length = axis_ticks_length
  )
p2_x <- p2 +
  scale_fill_manual(
    values = c(taxon_color$LEGEND_COLORS) %>%
      {
        names(.) <- taxon_color$Taxa_label
        .
      } %>%
      .[order(names(.))]
  ) +
  scale_x_discrete(
    limits = sample_meta_1 %>%
      with(paste0(get("Group"), get("x")))
  ) +
  theme(axis.text.x = element_text(color = sample_meta_col[sample_meta$Group]))


##### OUTPUT                                                               #####
pout <- p1s$jaccard + p1s$bray + p2_x +
  plot_layout(
    design = "AB\nCC",
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 18))
ggsave(filename = fig_out, plot = pout, width = 13, height = 10)
