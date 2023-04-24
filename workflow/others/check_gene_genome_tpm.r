###
#' @Date: 2022-07-17 09:41:43
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-08-05 14:13:41
#' @FilePath: /2021_09-MT10kSW/workflow/others/check_gene_genome_tpm.r
#' @Description:
###
#source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)
library(tidyverse)
library(ggplot2)
library(patchwork)

site_group = list(
  "Group" = c("Sw", "Sw", "Sw", "Bw", "Bw", "Bw", "Bw", "Bw", "Bw", "Bw", "Ss",
              "Ss", "Ss", "Ss", "Ss", "Ss", "Ss", "Ss", "Ss", "Bs", "Bs", "Bs",
              "Bs"),
  "Site" = c("TY.044", "TY.041", "TY.040", "YW.021", "WQ.022", "YW.020",
             "WQ.024", "YW.019", "WQ.021", "YW.023", "MC02", "D1T1", "D1T2",
             "T1B5", "T1B3", "T1B8", "T1L6", "T1B10", "T1B11", "T3L11", "T3L8",
             "T3L14", "T1L10")
) %>% as.data.frame()


### ######################################################################## ###
#### input data                                                             ####
### ######################################################################## ###
rep_sample_tpm =
  #file_path$file_path$cache("rep_sample_tpm.csv") %>% as.character() %>%
  "results/cache/rep_sample_tpm.csv" %>%
  read.csv %>%
  {.$wrong = ifelse(.$wrong == "", "true", .$wrong); .} %>%
  {.$Group = layer_2_group(.$Layer); .}

plot_reads_denisty <- function(data, x, y) {
  ggplot(data = data,
         mapping = aes_string(x = x, y = y)) +
    geom_point(mapping = aes_string(color = "Group", shape = "wrong"),
               alpha = 0.1) +
    scale_x_log10(labels =
                    ~format(.x, scientific = TRUE) %>%
                    str_replace("^0e\\+0", "0e+0") %>%
                    str_replace("e\\+0", "%*%10^") %>%
                    parse(text = .)) +
    scale_y_log10(labels =
                    ~format(.x, scientific = TRUE) %>%
                    str_replace("^0e\\+0", "0e+0") %>%
                    str_replace("e\\+0", "%*%10^") %>%
                    parse(text = .)) +
    scale_color_manual(values = sample_meta_col) +
    ggnewscale::new_scale_color() +
    stat_density2d(mapping = aes_string(linetype = "wrong", color = "..level..")) +
    scale_shape_manual(values = c("true" = 1, "ko" = 4))
}

ps =
  rep_sample_tpm %>% split(.$Group) %>%
  lapply(function(df) plot_reads_denisty(df, x = "count", y = "tpm")) %>%
  {.$Sw + .$Ss + .$Bw + .$Bs + plot_layout(ncol = 2)} %>%
  {ggsave("count_tpm.png", plot = ., width = 22, height = 22, dpi = 300)}
ps =
  rep_sample_tpm %>% split(.$Group) %>%
  lapply(function(df) plot_reads_denisty(df, x = "count", y = "rpb")) %>%
  {.$Sw + .$Ss + .$Bw + .$Bs + plot_layout(ncol = 2)} %>%
  {ggsave("count_rpb.png", plot = ., width = 22, height = 22, dpi = 300)}
ps =
  rep_sample_tpm %>% split(.$Group) %>%
  lapply(function(df) plot_reads_denisty(df, x = "rpb", y = "tpm")) %>%
  {.$Sw + .$Ss + .$Bw + .$Bs + plot_layout(ncol = 2)} %>%
  {ggsave("rpb_tpm.png", plot = ., width = 22, height = 22, dpi = 300)}
