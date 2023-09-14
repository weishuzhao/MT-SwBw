###
#' @Date: 2022-05-04 09:42:54
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2023-04-23 16:41:17
#' @FilePath: /2021_09-MT10kSW/workflow/others/draw_supp_fig4.r
#' @Description:
###
library(ggplot2)
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, keyword_args, fig_out_path                         #####
fig_out <- argv[1]

##### GLOBAL CONST vars                                                    #####

##### LOAD data AND transform TO basic format                              #####
gene_ko_tpm <- load__gene_ko_tpm()

key_genes <- load__key_genes()
TOTAL_GENE_NUM <- key_genes$KO %>% # nolint: object_name_linter.
  unique() %>%
  length() %>%
  {
    . * 2
  }


### ######################################################################## ###
#### Define function AND Calculate data                                     ####
### ######################################################################## ###
##### Define function                                                      #####
report__tpm_key_bar_group <- function(tpm_key) {
  tpm_key.signif <-
    tpm_key %>%
    location.group.signif("KO", "TPM", "Location", "Group", TOTAL_GENE_NUM) %>%
    merge(tpm_key[c("KO", "Label", "Pathway")] %>% unique(), by = "KO") %>%
    .[.$p.value.char != "-", ]

  tpm_key_1 <-
    tpm_key %>%
    {
      merge(
        merge(
          .[c("KO", "Label", "Pathway")] %>% unique(),
          .[c("Group", "Location")] %>% unique()
        ),
        .,
        all.x = TRUE
      )
    } %>%
    {
      .[is.na(.$TPM), "Hide"] <- "TRUE"
      .[is.na(.)] <- 1
      .
    }

  #### add paint                                                            ####
  font_size_1 <- 14
  font_size_2 <- 11
  font_size_3 <- 2.5
  axis.ticks.length <- 0.1 # nolint: object_name_linter.

  p <-
    ggplot(data = tpm_key_1) +
    geom_boxplot(mapping = aes_string(
      x = "Label", y = "TPM",
      fill = "Group",
      linetype = "Hide"
    )) +
    scale_linetype_manual(values = c("1" = 1, "TRUE" = 0)) +
    guides(linetype = "none") +
    scale_fill_manual(values = sample_meta_col) +
    facet_grid(formula("Location ~ Pathway"),
      scales = "free_x", space = "free_x"
    ) +
    geom_label(
      data =
        tpm_key.signif %>%
          .[.$p.value.char != "", ],
      mapping = aes_string(
        x = "Label", y = "Max. + max(Max.) * 0.05",
        label = "p.value.char"
      ),
      label.padding = unit(0.05, "lines"), label.size = 0,
      fill = "#7f7f7f3f", size = font_size_3
    ) +
    labs(x = "")

  p1 <-
    p +
    scale_y_log10(
      labels =
        ~ format(.x, scientific = TRUE) %>%
          str_replace("^0e\\+0", "0e+0") %>%
          str_replace("e\\+0", "%*%10^") %>%
          parse(text = .)
    ) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "white", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    ) +
    theme(
      axis.line = element_line(colour = "black"),
      axis.text = element_text(
        size = font_size_2, colour = "black", face = "bold"
      ),
      axis.title = element_text(
        size = font_size_1, face = "bold", colour = "black"
      ),
      axis.ticks.length = unit(axis.ticks.length, "cm"),
      axis.text.x = element_text(
        angle = 45, vjust = 1, hjust = 1, face = "plain"
      )
    ) +
    theme(
      legend.title = element_text(size = font_size_2, face = "bold"),
      legend.text = element_text(size = font_size_2, face = "bold")
    ) +
    theme(text = element_text(
      family = "Arial",
      size = font_size_1,
      hjust = 0.5,
      lineheight = 0.5
    )) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.key = element_blank()) +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.x = element_blank()
    )

  p1
}

##### Calculate data                                                       #####
tpm_key <-
  gene_ko_tpm %>%
  group_by(KO = get("KO"), Layer = get("Layer")) %>%
  summarise(TPM = sum(get("TPM"))) %>%
  merge(key_genes, by = "KO") %>%
  mutate(Site = gsub("^(.+)\\.\\.(.+)$", "\\1", get("Layer"))) %>%
  merge(site_group) %>%
  merge(unique(sample_meta[c("Location", "Group")])) %>%
  mutate(Location = factor(get("Location"), levels = c("Slope", "Bottom")))


### ######################################################################## ###
#### Plot figures and OUTPUT                                                ####
### ######################################################################## ###
##### Plot figures                                                         #####
# key_genes %>% .$Pathway %>% {data.frame(table(.), cumsum(table(.)))}
p1 <-
  tpm_key %>%
  {
    .[.$Pathway %in% c("3HP-4HB", "rTCA", "CBB"), ]
  } %>%
  {
    .$Pathway <- factor(.$Pathway, levels = unique(.$Pathway))
    .
  } %>%
  {
    report__tpm_key_bar_group(.) + guides(fill = "none")
  }
p2 <-
  tpm_key %>%
  {
    .[.$Pathway %in% c("D-AA", "L-sugar"), ]
  } %>%
  {
    .$Pathway <- factor(.$Pathway, levels = unique(.$Pathway))
    .
  } %>%
  {
    report__tpm_key_bar_group(.) + guides(fill = "none")
  }
p3 <-
  tpm_key %>%
  {
    .[.$Pathway %in% c("Alkane", "Aromatic", "Complex sugar"), ]
  } %>%
  {
    .$Pathway <- factor(.$Pathway, levels = unique(.$Pathway))
    .
  } %>%
  {
    report__tpm_key_bar_group(.) + guides(fill = "none")
  }
p4 <-
  tpm_key %>%
  filter(
    get("Pathway") %in%
      c("S", "As", "Se", "ROS", "Glycine betaine", "TMAO", "DMSP")
  ) %>%
  mutate(Pathway = factor(get("Pathway"), levels = unique(get("Pathway")))) %>%
  report__tpm_key_bar_group() +
  guides(fill = "none")



##### OUTPUT                                                               #####
c() %>%
  {
    p1 + p2 + p3 + p4 +
      plot_layout(design = "AAAB
                            CCCC
                            DDDD") +
      plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
      theme(plot.tag = element_text(size = 18))
  } %>%
  {
    ggsave(
      filename = fig_out,
      plot = .,
      width = 18, height = 18, dpi = 300
    )
  }
