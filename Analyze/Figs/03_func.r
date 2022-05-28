###
#* @Date: 2021-08-22 21:53:56
#* @LastEditors: Hwrn
#* @LastEditTime: 2021-10-27 11:26:22
#* @FilePath: /2021_09-MT10kSW/Analyze/Figs/03_func.r
#* @Description:
###
source('Analyze/Figs/02_func.r')
source('Analyze/Figs/02_load.r')
RELOAD = FALSE
#RELOAD = TRUE


summarize__Wtdb <- function(Wtdb) {
  p = ggplot(data = Wtdb) +
    geom_point(mapping = aes_string(x = "contamination", y = "completeness",
                                    color = "annot_level", size = "score"),
               alpha = 0.5) +
    scale_color_hue(h=c(0, 270) + 90, c=100, l=70, direction = 1) +
    theme_classic(base_line_size = 0) +
    scale_x_continuous(limits = c(0, 10)) +
    scale_y_continuous(limits = c(50, 100))
  p +
    annotate("segment", x =  0, xend = 10, y = 75, yend =  75) +
    annotate("segment", x =  0, xend = 10, y = 90, yend =  90) +
    annotate("segment", x =  5, xend =  5, y = 90, yend = 100) +
    annotate("segment", x =  0, xend =  5, y = 95, yend =  95)

  sumqual = cumsum(rev(summary(Wtdb$comp_cont)))
  sumqual = data.frame(
    sumqual,
    sumqual.rate = sumqual/length(Wtdb$genome))

  sumphylum = summary(taxon.split(as.character(Wtdb$taxonomy), 2))
  sumphylum = data.frame(
    sumphylum,
    sumphylum.rate = sumphylum/length(Wtdb$genome)
  )

  return(list(p, sumqual, sumphylum))
}


summarize__Wtdb.sample <- function(Wtdb, annotation_row) {
  Wtdb.t = t(Wtdb[apply(Wtdb[, rownames(annotation_row)], 1, sum) > 1,
                  rownames(annotation_row)])
  #Wtdb.t[Wtdb.t == 0] = NA
  p = pheatmap(Wtdb.t,
               show_colnames = FALSE,
               annotation_row = annotation_row)

  sumsingle = data.frame(row.names = rownames(annotation_row))
  sumsingle$sumsingle = apply(Wtdb[apply(Wtdb[, rownames(sumsingle)
                                              ], 1, sum) == 1,
                                   rownames(sumsingle)],
                              2, sum)
  sumsingle$sum = apply(Wtdb[, rownames(annotation_row)], 2, sum)
  sumsingle$sumsingle.rate = sumsingle$sumsingle / sumsingle$sum

  tax = Wtdb[apply(Wtdb[, rownames(annotation_row)], 1, function(x) {
    ifco = 0
    for (site in levels(as.factor(annotation_row$site))) {
      ifco = ifco + (sum(x[annotation_row$site == site]) > 0)}
    return(ifco > 1)
  }), "taxonomy"]

  return(list(p, sumsingle, tax))
}


collapse__Wtdb <- function(Wtdb, annotation_row,
                           taxon.level = 7, FUN = sum) {
  # fold Wtdb by taonomy
  Wtdb$taxonomy = as.character(taxon.split(Wtdb$taxonomy, 1, taxon.level))
  Wtdb = Wtdb[strsplit(Wtdb$taxonomy, "__$") == Wtdb$taxonomy, ]

  Wtdb.t = data.frame(
    Wtdb %>%
      group_by(get("taxonomy")) %>%
      mutate(nculs = n()) %>%
      filter(get("nculs") > 1) %>%
      filter(score == max(get("score"))) %>%
      arrange(get("taxonomy"))
  )

  if (dim(Wtdb.t)[1] == 0) {
    warning("Nothing left after filtering")
    return()
  }

  Wtdb.t[, rownames(annotation_row)] = t(sapply(
    Wtdb.t$taxonomy, function(x) apply(
      Wtdb[Wtdb$taxonomy == x, rownames(annotation_row)], 2, FUN))
  )[, rownames(annotation_row)]

  return(Wtdb.t)
}


report__div.Wtdb <- function(Wtdb, value = "score", FUN = sum) {
  order.Wtdb = unique(sort(as.character(Wtdb$dpc)))
  div.16s.Wtdb = data.frame(
    dpc = rep(order.Wtdb, length(levels(div.16s$sample))),
    sample = unlist(lapply(levels(div.16s$sample),
                           function(x) rep(x, length(order.Wtdb))))
  )
  div.16s.Wtdb$name = taxon.split(div.16s.Wtdb$dpc, 3)

  # check taxonomy
  div.16s.Wtdb$name.16s = as.character(div.16s.Wtdb$name)
  unique(sort(Wtdb$taxonomy[
    taxon.split(Wtdb$taxonomy, 3) %in% div.16s.Wtdb$name.16s[
      !div.16s.Wtdb$name.16s %in% taxon.order$name]]))
  # https://doi.org/10.1099/ijsem.0.003920
  div.16s.Wtdb$name.16s[div.16s.Wtdb$name.16s == "Actinomycetia"
  ] = "Actinobacteria"

  # map div.16s to Bin diversity
  div.16s.g = taxon.pct.annot(div.16s, taxon = "name",
                              taxon.sort = div.16s.Wtdb$name.16s)
  div.16s.Wtdb$precent.16s = apply(
    div.16s.Wtdb[, c("sample", "name.16s")], 1,
    function(x) sum(div.16s.g[div.16s.g$sample == x[[1]] &
                                div.16s.g$name == x[[2]],
                              "annot.percent"]))
  div.16s.Wtdb$precent.16s[is.na(div.16s.Wtdb$precent.16s)] = 0

  Stdb.dpc = taxon.split(Stdb$taxonomy, 1, 3)  # speed up
  div.16s.Wtdb$value = apply(
    div.16s.Wtdb[, c("sample", "dpc")], 1,
    function(x) FUN(Stdb[Stdb$sample == x[1] &
                           Stdb.dpc == x[2],
                         value]))

  ggplot(data = div.16s.Wtdb[!div.16s.Wtdb$dpc == "", ]) +
    #geom_bar(mapping = aes(x = order, fill = domain), y = "TY.040", width = 1) +
    geom_point(mapping = aes_string(x = "name.16s", y = "sample",
                                    size = "value",
                                    color = "log(precent.16s + 1)"),
               alpha = 0.5) +

    # considering the smallest value must be 0
    scale_size_continuous(range = c(-1, 10)) +
    scale_color_distiller(palette = "PuBuGn", direction = 1) +
    guides(size = guide_legend(title = paste0(deparse(substitute(FUN)),
                                              "(", value, ")")),
           color = guide_colorbar(
             title = "abundance\n(log(reads %) from 16s)")) +

    labs(title = "", x = "class", y = "sample") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

}


report__div.Stdb <- function(Wtdb.cluster) {
  # give a subset of Wtdb, use cluster
  # use environment Stdb
  # A Dot-Matrix plot: 2nd_Cluster (taxon) * sample, size = depth, color = comp_cont
  Stdb.c = Stdb[Stdb$secondary_cluster %in% Wtdb.cluster,]
  Stdb.c$taxonomy = sapply(Stdb.c$secondary_cluster, function(x) {
    Wtdb[Wtdb$cluster == x, ]$taxonomy
  })
  Stdb.c$comp_cont = "low quality"
  for (i in comp_cont) {
    Stdb.c$comp_cont[Stdb.c$Completeness >= i[1] & Stdb.c$Contamination <= i[2]
                   ] = paste0(">", i[1], ", <", i[2])}
  Stdb.c$comp_cont = factor(Stdb.c$comp_cont, levels = sapply(
    comp_cont, function(i) paste0(">", i[1], ", <", i[2])
  ))
  Stdb.c$species = taxon.split(Stdb.c$taxonomy, 7)
  colnames(Stdb.c)

  tmp = table(Stdb.c[Stdb.c$comp_cont %in% c(">90, <5", ">95, <5"), "taxonomy"]) >= 3
  Stdb.c[Stdb.c$taxonomy %in% names(tmp[tmp]),]
  
  ggplot(data =Stdb.c[Stdb.c$taxonomy %in% names(tmp[tmp]),]) +
    #geom_bar(mapping = aes(x = order, fill = domain), y = "TY.040", width = 1) +
    geom_point(mapping = aes_string(x = "sample", y = "taxon.split(taxonomy, 7)",
                                    size = "AvgDepth",
                                    color = "comp_cont"),
               alpha = 0.5) +
    scale_color_hue(h=c(0, 270) + 90, c=100, l=70, direction = -1) +

    theme_classic() +
    theme(axis.text.y = element_text(hjust = 0),
          axis.text.x = element_text(angle = 90, hjust = 1))

}
