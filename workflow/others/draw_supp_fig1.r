###
#' @Date: 2022-07-20 13:43:25
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2024-05-24 17:44:16
#' @FilePath: /2021_09-MT10kSW/workflow/others/draw_supp_fig1.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, fig_out_path, keyword_args                         #####
otu_class_abd <- argv[1]
#' otu_class_abd <- "results/reads_diversity/cross_abundance.csv"
fig_out <- argv[2]

##### GLOBAL CONST vars                                                    #####
font_size_1 <- 13
font_size_2 <- 10
font_size_3 <- 10
axis_ticks_length <- 0.1

##### LOAD data AND transform TO basic format                              #####
sample_meta_cross <-
  read.csv("results/reads_diversity/metadata.tsv", sep = "\t", as.is = TRUE) %>%
  mutate(
    X = get("Sample"),
    Layer = ifelse(
      get("Type") == "16s", get("Sample"), gsub("^[^_]+_", "", get("Sample"))
    ),
    Sample = paste0(get("Group"), "_", get("Layer"))
  ) %>%
  merge(unique(sample_meta[c("Location", "Group")])) %>%
  mutate(Group = factor(get("Group"), names(sample_meta_col))) %>%
  .[order(.$Group, .$Sample), ] %>%
  mutate(
    Group = as.character(get("Group")),
    Layer = .data$Layer %>%
      `[`(
        c(
          "TY.040" = "Sw3",
          "TY.035" = "Ss1",
          "TY.044" = "Ss2",
          "TY.048" = "Ss3",
          "TY.041" = "Ss6",
          "WQ.022" = "Bw1",
          "YW.019" = "Bw3",
          "YW.020" = "Bw5",
          "YW.023" = "Bw7",
          "TY.042" = "Bs1",
          "TY.038" = "Bs2",
          "TY.039" = "Bs15",
          "TY.046" = "Bs21"
        ),
        .
      )
  )

cross_rltabd <-
  otu_class_abd %>%
  read.csv() %>%
  merge(sample_meta_cross[c("X", "Sample")]) %>%
  mutate(X = NULL) %>%
  column_to_rownames("Sample")

### ######################################################################## ###
#### Define function AND Calculate data                                     ####
### ######################################################################## ###
##### Define function                                                      #####

##### Calculate data                                                       #####

### ######################################################################## ###
#### Plot figures and OUTPUT                                                ####
### ######################################################################## ###
##### Plot figures                                                         #####
p2 <-
  #' otu_rltabd %>%
  #' mutate(Taxa_label = taxon.split(get("Taxonomy"), 1, 3)) %>%
  cross_rltabd %>%
  rownames_to_column("Sample") %>%
  pivot_longer(
    !c("Sample"),
    names_to = "Taxa_label", values_to = "Abundance"
  ) %>%
  merge(sample_meta_cross) %>%
  filter(get("Type") == "16S") %>%
  as.data.frame() %>%
  get_percent_plot(
    "Abundance",
    fill.name = "Taxa_label",
    sample.name = "Layer",
    labs.x = "16S sample", labs.y = "relative abundance",
    TOP_N_TAXON_PER_LAYER = 10,
    font_size_1 = font_size_1, font_size_2 = font_size_2,
    font_size_3 = font_size_3, axis.ticks.length = axis_ticks_length
  )
p2_x <- p2 +
  scale_x_discrete(
    limits = filter(sample_meta_cross, get("Type") == "16S")$Layer
  ) +
  guides(fill = guide_legend(title = "Class", ncol = 3, reverse = TRUE)) +
  theme(
    axis.text.x = element_text(
      color = sample_meta_col[
        filter(sample_meta_cross, get("Type") == "16S")$Group
      ]
    )
  )

##### OUTPUT                                                               #####
ggsave(filename = fig_out, plot = p2_x, width = 9, height = 4)
