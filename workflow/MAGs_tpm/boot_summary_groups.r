###
#' @Date: 2022-10-09 16:13:56
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2023-04-22 16:48:30
#' @FilePath: /2021_09-MT10kSW/workflow/MAGs_tpm/boot_summary_groups.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, keyword_args, fig_out_path                         #####
fig_out <- argv[1]
icamp_file <- argv[2]
# icamp_file = "results/MAGs/icamp.iCAMP.BootSummary.Groups.RPKM.csv" # nolint

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
  .[.[, "Process"] != "Stochasticity", ] %>%
  {
    .[, "Process"] <-
      .[, "Process"] %>% factor(levels = names(assebly_factor_col))
    .
  }

icamp %>%
  head()


p <-
  icamp %>%
  {
    ggplot(
      data = .,
      mapping = aes_string(x = 1, y = "Mean", fill = "Process")
    ) +
      geom_bar(
        stat = "identity", position = position_stack(reverse = TRUE),
        color = "black"
      ) + # stack
      coord_polar(theta = "y") +
      scale_fill_manual(values = assebly_factor_col) +
      facet_grid(formula(". ~ Group")) +
      expand_limits(x = c(-0.5, 0.5)) +
      theme_void()
  }
# p
ggsave(filename = fig_out, plot = p, width = 15, height = 2, limitsize = FALSE)
