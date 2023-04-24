###
#' @Date: 2022-05-04 09:42:54
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2023-04-23 10:31:58
#' @FilePath: /2021_09-MT10kSW/workflow/MAGs_tpm/ko_tpm_MAG_N_pie.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, keyword_args, fig_out_path                         #####
fig_out <- argv[1]

##### LOAD data AND transform TO basic format                              #####
gene_ko_tpm <- load__gene_ko_tpm()

key_genes <- load__key_genes()
key_genes_n <-
  key_genes %>%
  filter(get("Pathway") == "N") %>%
  mutate(Pathway = get("Arrow") %>% factor(levels = unique(.)))

nxr_genomes <-
  "data/nxrAB/final_nxrAB.csv" %>%
  read.csv() %>%
  mutate(Genome = get("genome"), KO = get("ko")) %>%
  split(.$KO) %>%
  lapply(. %>% .$Genome %>% unique())
### ######################################################################## ###
#### Define function AND Calculate data                                     ####
### ######################################################################## ###
##### Define function                                                      #####

##### modify data for analysis                                             #####
tpm_key_n_sum_mean <-
  gene_ko_tpm %>%
  merge(key_genes_n, by = "KO") %>%
  mutate(
    Label = ifelse(
      get("KO") == "K00370",
      ifelse(
        get("Genome") %in% nxr_genomes[["K00370"]],
        "nxrA", "narG"
      ),
      ifelse(
        get("KO") == "K00371",
        ifelse(
          get("Genome") %in% nxr_genomes[["K00371"]],
          "nxrB", "narH"
        ),
        get("Label") %>% as.character()
      )
    ) %>%
      factor(levels = unique(.)),
    Pathway = ifelse(
      get("Label") %in% c("narG", "narH"),
      "Dissimilatory Nitrate -> Nitrite",
      ifelse(
        get("Label") %in% c("nxrA", "nxrB"),
        "Nitrite -> Nitraite",
        as.character(get("Pathway"))
      )
    ) %>%
      factor(levels = unique(.))
  ) %>%
  mutate(Site = gsub("^(.+)\\.\\.(.+)$", "\\1", get("Layer"))) %>%
  merge(site_group) %>%
  merge(unique(sample_meta[c("Location", "Group")])) %>%
  mutate(Group = factor(get("Group"), levels = names(sample_meta_col))) %>%
  split(paste(.$Label, .$Layer)) %>%
  lapply(function(df) {
    data.frame(
      df[1, c(
        "Group", "Site", "KO", "Layer", "Label", "Arrow", "Pathway", "Location"
      )],
      TPM = sum(df$TPM)
    )
  }) %>%
  bind_rows() %>%
  split(paste(.$Group, .$Label)) %>%
  lapply(function(df) {
    data.frame(
      df[1, c("Group", "KO", "Label", "Arrow", "Pathway", "Location")],
      TPM = mean(df$TPM)
    )
  }) %>%
  bind_rows() %>%
  .[order(.$Group, .$Arrow, -.$TPM), ]
tpm_key_n_sum_mean_width <-
  tpm_key_n_sum_mean %>%
  split(.$Label) %>%
  lapply(
    . %>% data.frame(Width = sum(.$TPM))
  ) %>%
  bind_rows() %>%
  mutate(
    width = get("Width") %>%
      log10() %>%
      ceiling() %>%
      symnum(
        cutpoints = c(0, 2, 20, 200, 2000, 5000),
        symbols = c(1, 2, 4, 8, 12)
      ) %>%
      as.integer(),
    Width = get("width") / get("Width")
  ) %>%
  split(.$Pathway) %>%
  lapply(
    . %>%
      {
        a <- unique(.$Label)
        .$iLabel <-
          .$Label %>%
          sapply(. %>%
            {
              . == a
            } %>%
            which())
        .
      }
  ) %>%
  bind_rows()

p <-
  tpm_key_n_sum_mean_width %>%
  ggplot(
    data = .,
    mapping = aes_string(
      x = "Group", y = "log10(TPM + 1)", # "TPM * Width",  #
      color = "Group", fill = "Group"
    )
  ) +
  geom_bar(
    stat = "identity", position = "dodge",
    width = 1, alpha = 0.3
  ) + # stack
  geom_text(
    aes_string(label = "Label", y = 0),
    color = "black"
  ) + # stack
  scale_color_manual(values = sample_meta_col) +
  scale_fill_manual(values = sample_meta_col) +
  coord_polar(theta = "x") +
  facet_grid(formula("Pathway ~ iLabel")) +
  guides(fill = "none", color = "none") +
  theme_void()
# p

ggsave(filename = fig_out, plot = p, width = 5, height = 12, limitsize = FALSE)
