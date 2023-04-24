###
#' @Date: 2022-10-09 16:13:56
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-11-05 15:35:55
#' @FilePath: /2021_09-MT10kSW/workflow/MAGs_tpm/boot_summary_groups_bar.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, keyword_args, fig_out_path                         #####
fig_out = argv[1]
icamp_file = argv[2]
#icamp_file = "results/MAGs/icamp.iCAMP.BootSummary.Groups.RPKM.csv"

##### GLOBAL CONST vars                                                    #####
assebly_factor_col = c(
  "Heterogeneous.Selection" = "#FF9D9A",
  "Homogeneous.Selection" = "#E15759",
  "Dispersal.Limitation" = "#00CCCC",
  "Homogenizing.Dispersal" = "#006666",
  "Drift.and.Others" = "#666666"
)

##### LOAD data AND transform TO basic format                              #####
icamp =
  icamp_file %>%
  read.csv %>%
  .[.[, "Process"] != "Stochasticity", ] %>%
  {.[, "Process"] = .[, "Process"] %>% factor(levels = names(assebly_factor_col)); .} %>%
  {
    .[, "Group"] =
      .[, "Group"] %>%
      strsplit("_vs_") %>%
      sapply(sort) %>%
      sapply(. %>% paste(collapse = " vs "))
    .
  } %>%
  {
    .$Compare =
      .$Group %>%
      grepl(" vs ", .) %>% ifelse("compare", "single") %>%
      factor(c("single", "compare"))
    .
  }

#icamp %>%
#  head
  #group_by(Group) %>%   summarise(Mean = sum(Mean))
#.[, "Process"] %>% factor(levels = names(assebly_factor_col))


p =
  icamp %>%
  {
    font_size_1 = 16
    font_size_2 = 14
    axis.ticks.length = 0.1

    ggplot(data = .,
           mapping = aes_string(x = "Group", y = "Mean", fill = "Process")) +
      geom_bar(stat = "identity", position = position_stack(reverse = T), color = "black") +  # stack
      scale_fill_manual(values = assebly_factor_col) +
      facet_grid(formula(". ~ Compare"),
                 scales = "free_x", space = "free_x") +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = 'white', size = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()
      ) +
      theme(
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = font_size_2, colour = "black", face = "bold"),
        axis.title = element_text(size = font_size_1, face = "bold", colour = "black"),
        axis.ticks.length = unit(axis.ticks.length, 'cm'),

        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "plain")
      ) +
      theme(
        legend.title = element_text(size = font_size_2, face = "bold"),
        legend.text = element_text(size = font_size_2, face = "bold")
        #legend.position = "bottom"
      ) +
      theme(text = element_text(family = "Arial",
                                size = font_size_1,
                                hjust = 0.5,
                                lineheight = 0.5)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.key = element_blank()) +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            strip.text.x = element_blank())  }
#p
ggsave(filename = fig_out, plot = p, width = 10, height = 6, limitsize = FALSE)
