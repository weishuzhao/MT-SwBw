###
#' @Date: 2022-10-09 16:13:56
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2023-09-20 19:46:06
#' @FilePath: /2021_09-MT10kSW/workflow/MAGs_tpm/boot_summary_groups_bar.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, keyword_args, fig_out_path                         #####
fig_out <- argv[1]
icamp_file <- argv[2]
#' icamp_file = "results/MAGs/icamp.iCAMP.BootSummary.Groups.RPKM.csv"
icamp_file2 <- argv[3]
#' icamp_file2 = "results/reads_diversity/ProcessImportance_16S.csv"

##### GLOBAL CONST vars                                                    #####
assebly_factor_col <- c(
  "Heterogeneous.Selection" = "#FF9D9A",
  "Homogeneous.Selection" = "#E15759",
  "Dispersal.Limitation" = "#00CCCC",
  "Homogenizing.Dispersal" = "#006666",
  "Drift.and.Others" = "#666666"
)

##### LOAD data AND transform TO basic format                              #####
icamp <-
  icamp_file %>%
  read.csv() %>%
  filter(get("Process") != "Stochasticity") %>%
  mutate(
    Process = factor(get("Process"), levels = names(assebly_factor_col)),
    Group = get("Group") %>%
      strsplit("_vs_") %>%
      sapply(sort) %>%
      sapply(. %>% paste(collapse = " vs ")),
    Compare = get("Group") %>%
      grepl(" vs ", .) %>%
      ifelse("compare", "single") %>%
      factor(c("single", "compare"))
  ) %>%
  dplyr::select("Group", "Mean", "Process", "Compare")

icamp2 <-
  icamp_file2 %>%
  read.csv() %>%
  pivot_longer(
    !c("Group", "GroupBasedOn", "Method"),
    names_to = "Process_raw",
    values_to = "Mean"
  ) %>%
  mutate(
    Process = get("Process_raw") %>%
      factor(c("HeS", "HoS", "DL", "HD", "DR")) %>%
      {
        levels(.) <- names(assebly_factor_col)
        .
      },
    Group = get("Group") %>%
      strsplit("_vs_") %>%
      sapply(sort) %>%
      sapply(. %>% paste(collapse = " vs ")),
    Compare = get("Group") %>%
      grepl(" vs ", .) %>%
      ifelse("compare", "single") %>%
      factor(c("single", "compare"))
  ) %>%
  filter(get("GroupBasedOn") == "Group") %>%
  dplyr::select("Group", "Mean", "Process", "Compare")

p <-
  list("Metagenome" = icamp, "16S rRNA" = icamp2) %>%
  bind_rows(.id = "DataType") %>%
  filter(get("Compare") == "compare") %>%
  mutate(DataType = factor(get("DataType"), c("Metagenome", "16S rRNA"))) %>%
  ggplot(
    data = .,
    mapping = aes_string(x = "Group", y = "Mean", fill = "Process")
  ) %>%
  {
    font_size_1 <- 16
    font_size_2 <- 14
    axis_ticks_length <- 0.1

    . +
      theme(
        axis.line = element_line(colour = "black"),
        axis.text =
          element_text(size = font_size_2, colour = "black", face = "bold"),
        axis.title =
          element_text(size = font_size_1, face = "bold", colour = "black"),
        axis.ticks.length = unit(axis_ticks_length, "cm"),
        axis.text.x =
          element_text(angle = 45, vjust = 1, hjust = 1, face = "plain")
      ) +
      theme(
        legend.title = element_text(size = font_size_2, face = "bold"),
        legend.text = element_text(size = font_size_2, face = "bold")
        #' legend.position = "bottom"
      ) +
      theme(
        text = element_text(
          family = "Arial", size = font_size_1, hjust = 0.5, lineheight = 0.5
        )
      )
  } +
  geom_bar(
    stat = "identity", position = position_stack(),
    color = "black"
  ) + # stack
  scale_fill_manual(values = assebly_factor_col) +
  facet_grid(
    formula(". ~ DataType"),
    scales = "free_x", space = "free_x"
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "white", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.key = element_blank()) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside"
  )
# p
ggsave(filename = fig_out, plot = p, width = 8, height = 6, limitsize = FALSE)
