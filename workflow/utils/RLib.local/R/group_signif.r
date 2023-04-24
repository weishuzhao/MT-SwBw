###
#' @Date: 2022-07-19 10:49:43
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-07-19 10:51:12
#' @FilePath: /2021_09-MT10kSW/workflow/utils/RLib.local/R/group_signif.r
#' @Description:
###

location.group.signif <-
  function(ko_tpm_group,
           ko, tpm, location = "location", group = "group",
           TOTAL_GENE_NUM) {
    ko_tpm_group %>%
      split(.[ko]) %>%
      lapply(
        function(x) {
          x %>%
            split(.[location]) %>%
            lapply(
              function(x1) {
                p.value = ifelse(
                  length(unique(x1[, group])) == 2,
                  wilcox.test(formula(paste(tpm, "~", group)), data = x1,
                              paired = FALSE)$p.value,
                  2
                )
                c("p.value" = p.value, summary(x1[, tpm]))
              }
            ) %>%
            bind_rows(.id = location)
        }
      ) %>%
      bind_rows(.id = ko) %>%
      {
        .$p.value.adj =
          .$p.value %>%
          {ifelse(. == 2, 2, p.adjust(., "BH", n = TOTAL_GENE_NUM * 2))}
        .
      } %>%
      {
        .$p.value.char =
          as.character(symnum(.$p.value.adj,
                              cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1, 2),
                              symbols = c("***", "**", "*", ".", "", "-")))
        .
      }
  }

get_ko.diff.annots <- function(tpm.key_or_tpm.key.signif, TOTAL_GENE_NUM) {
  if ("p.value.char" %in% tpm.key_or_tpm.key.signif) {
    tpm.key.signif = tpm.key_or_tpm.key.signif
  } else {
    tpm.key.signif = location.group.signif(tpm.key_or_tpm.key.signif,
                                           "ko", "tpm", "location", "group",
                                           TOTAL_GENE_NUM)
  }
  tpm.key.signif[c("ko", "location", "p.value.char")] %>%
    merge(key_genes, by.x = "ko", by.y = "KO") %>%
    reshape2::acast(., formula("ko ~ location"), value.var = "p.value.char") %>%
    apply(
      1,
      function(x) ifelse(
        all(grepl("\\*{3}", x[1]), grepl("\\*{3}", x[2])), "all.signif",
        ifelse(
          all(grepl("\\*+", x[1]), grepl("^$", x[2])), "mixed",
          "others"
        )
      )
    ) %>%
    merge(key_genes, by.x = 0, by.y = "KO") %>%
    {
      data.frame(
        pathway = .$pathway, label = .$label, new = .$new, diff = .$x,
        row.names = .$Row.names
      )
    }
}
