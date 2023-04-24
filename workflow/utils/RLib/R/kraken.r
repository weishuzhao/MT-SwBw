###
#* @Date: 2022-02-26 22:31:37
#* @LastEditors: Hwrn
#* @LastEditTime: 2022-02-27 13:40:03
#* @FilePath: /metaSC/R/RLib/kraken.r
#* @Description:
###

report__kraken_rare <- function(kraken_rare, report_path_to_sample){
  kraken_rare$sample = report_path_to_sample(kraken_rare$sample)

  kraken_rare$location = sapply(
    kraken_rare$sample,
    function(x) {unlist(strsplit(x, "\\_"))[1]})

  kraken_rare_stat = as.data.frame(t(sapply(
    split(kraken_rare[c("sample_size", "otus", "sample")],
          kraken_rare$sample),
    function(x) x[which.max(unlist(x["sample_size"])),])))
  kraken_rare_stat$sample_size = unlist(kraken_rare_stat$sample_size)
  kraken_rare_stat$otus = unlist(kraken_rare_stat$otus)
  kraken_rare_stat$sample = unlist(kraken_rare_stat$sample)
  colnames(kraken_rare)
  ggplot(data = kraken_rare,
         mapping = aes_string(x = "sample_size", y = "otus",
                              color = "sample")) +
    stat_summary(geom = "line", fun = "mean") +
    stat_summary(mapping = aes_string(color = "sample"),
                 geom = "errorbar",
                 fun.min = function(x) mean(x - sd(x)/sqrt(length(x))),
                 fun.max = function(x) mean(x + sd(x)/sqrt(length(x)))) +
    stat_summary(geom = "point", fun = "mean") +
    geom_text_repel(
      data = kraken_rare_stat,
      mapping = aes_string(x = "sample_size", y = "otus",
                           label = "sample")
    ) +
    #scale_x_continuous(
    #  expand = expansion(),
    #  labels = ~format(.x, scientific = TRUE) %>%
    #    str_replace("1e\\+0", "10^") %>%
    #    #str_replace("e\\+0", "*10^") %>%
    #    #str_replace("e\\-0", "*10^-") %>%
    #    parse(text = .)
    #) +
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
    )) +
    guides(fill = guide_legend(reverse = TRUE))

}


exec_kraken_rare <- function(file, step=2e6) {
  raw = read.csv(file, sep = "\t",
                 col.names = c("percents", "all_reads", "exact_reads",
                               "level_type", "level_id", "name"))
  kraken_counts = raw[sapply(raw["level_type"],
                             function(x) gsub("(\\w).*", "\\1", x) != "R"),
                      "exact_reads"]
  cols = c(seq(step, sum(kraken_counts), step), sum(kraken_counts))
  data.frame(
    sample_size = cols,
    otus = sapply(cols, function(x) as.integer(rarefy(t(kraken_counts), x))),
    sample = file
  )

}

# demo
if (FALSE) {
  kraken_report_files = paste0(
    "kraken-",
    c(
      "S1", "S2", "S2B", "S3",
    ),
    ".report"
  )
  report_path_to_sample <- function(x) {
    x = gsub("kraken.S(.+)B.report", "B_\\1", x)
    x = gsub("kraken.S(.+).megahit.report", "A_\\1", x)
    return(x)
  }

}

report__kraken_rare_from_file <- function(
  kraken_report_files,
  report_path_to_sample
) {
  #kraken_rare = read.csv("kraken_rare.csv",
  #                       stringsAsFactors = FALSE)
  kraken_rare = data.frame()
  for (file in kraken_report_files) {
    kraken_rare = rbind(kraken_rare, exec_kraken_rare(file))
  }

  report__kraken_rare(kraken_rare, report_path_to_sample)

}
