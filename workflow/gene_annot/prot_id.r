###
#* @Date: 2022-03-20 16:47:38
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-05-20 10:45:57
#' @FilePath: /2022_05-ZFMG-release/workflow/gene_annot/prot_id.r
#* @Description:
###
source("Scripts/utils/RLib.r")

get_blastp_tsv <- function(file) {
  MT2ME = read.csv(
    file,
    sep = "\t",
    header = FALSE,
    col.names = unlist(
      strsplit(
        "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "\\s"
      )
    )
  )

  MT2ME_max =
    MT2ME[,c("qseqid", "pident")] %>%
    group_by(get("qseqid")) %>%
    filter(get("pident") == max(get("pident"))) %>%
    ungroup()

  MT2ME_max
}
cross_max = bind_rows(get_blastp_tsv("workflow/gene_annot/faa/MT.faa.blastp.tsv"),
                      get_blastp_tsv("workflow/gene_annot/faa/ME.faa.blastp.tsv"))


p = ggplot(data = cross_max) +
  geom_histogram(mapping = aes_string(x = "pident"),
                 breaks = seq(26, 100, 1)) +
  scale_y_continuous(labels = ~format(.x, scientific = TRUE) %>%
                       str_replace("e\\+0", "%*%10^") %>%
                       parse(text = .)) +
  labs(x = "highest identity") +

  scale_x_continuous(breaks = seq(30, 100, 10)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5
  ))

p = p +
  geom_vline(xintercept = c(30, 50, 90), linetype = 2:4,
             size = 0.7)

font_size_1 = 15
font_size_2 = 13
axis.ticks.length = 0.3

p = p +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black",
                             face = "bold", size = font_size_2),
    axis.title = element_text(
      size = font_size_1, face = "bold", colour = "black"),
    legend.title = element_text(size = font_size_2, face = "bold"),
    legend.text = element_text(size = font_size_2, face = "bold")
    #legend.position = "bottom"
  )


ggsave(filename = "figs/S5A.svg",
       plot = p, width = 7.5, height = 6, dpi = 300)

################################################################################

rm(list = ls(all.names = TRUE))

get_blastn_tsv <- function(file) {
  print(file)
  MT2ME = read.csv(
    file,
    sep = "\t",
    header = FALSE,
    col.names = unlist(
      strsplit(
        "qaccver saccver pident length mismatch gapopen qstart qend sstart send
   evalue bitscore",
        "\\s+"
      )
    )
  )

  MT2ME_max =
    MT2ME[,c("qaccver", "pident")] %>%
    as_tibble() %>%
    group_by(get("qaccver")) %>%
    filter(get("pident") == max(get("pident"))) %>%
    ungroup()

  MT2ME_max
}


cross_max_n = bind_rows(get_blastn_tsv("workflow/gene_annot/fna/MT.fna.blastn.tsv"),
                        get_blastn_tsv("workflow/gene_annot/fna/ME.fna.blastn.tsv"))

cross_max_n %>%
  split(x = .$pident,
        f = sapply(.$qaccver,
                   function(x) gsub("^([^|]+)\\|(.+)$", "\\1", x))) %>%
  lapply(., function(x) length(x)) %>%
  as.data.frame(.) %>%
  t(.)


p = ggplot(data = cross_max_n) +
  geom_histogram(mapping = aes_string(x = "pident"),
                 breaks = seq(68, 100, 0.5)) +
  scale_y_continuous(labels = ~format(.x, scientific = TRUE) %>%
                       str_replace("e\\+0", "%*%10^") %>%
                       parse(text = .)) +
  labs(x = "highest identity") +

  scale_x_continuous(breaks = seq(70, 100, 5)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5
  ))

p = p +
  geom_vline(xintercept = c(90, 95), linetype = 2:3,
             size = 0.7)

font_size_1 = 15
font_size_2 = 13
axis.ticks.length = 0.3

p = p +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black",
                             face = "bold", size = font_size_2),
    axis.title = element_text(
      size = font_size_1, face = "bold", colour = "black"),
    legend.title = element_text(size = font_size_2, face = "bold"),
    legend.text = element_text(size = font_size_2, face = "bold")
    #legend.position = "bottom"
  )


ggsave(filename = "figs/S5B.svg",
       plot = p, width = 7.5, height = 6, dpi = 300)
