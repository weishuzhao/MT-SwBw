###
#* @Date: 2021-09-20 17:32:57
#* @LastEditors: Hwrn
#* @LastEditTime: 2021-11-10 18:56:53
#* @FilePath: /2021_09-MT10kSW/Analyze/Figs/04_func.r
#* @Description:
###


report__abd.heat <- function() {
  annot.C = unique(module.name[,c('B', 'A', 'C')])
  row.names(annot.C) = annot.C$C; annot.C$C = NULL
  pheatmap(annot_for_heat(existance, module.name,
                          function(x){sum(x)/length(x)}),
           cluster_rows = FALSE, cutree_cols = 2,
           display_numbers = annot_for_heat(existance, module.name,
                                            function(x){paste(sum(x), length(x),
                                                              sep = "/")}),
           annotation_row = annot.C,
           scale = 'none')

  pheatmap(annot_for_heat(abundance, module.name,
                          function(x){sum(x)/length(x)}),
           cluster_rows = FALSE, cutree_cols = 2,
           display_numbers = t(abd.sgf), annotation_row = annot.C,
           scale = 'none')

  pheatmap(scale(annot_for_heat(abundance, module.name,
                                function(x){sum(x)/length(x)})),
           cluster_rows = FALSE, cutree_cols = 2,
           display_numbers = T, annotation_row = annot.C,
           scale = 'none')

  return(list(annot.C = annot.C))
}


report__pcoa.KO <- function(div.table) {
  pcoa.16s = pcoa(vegdist(div.table, method = "bray"))
  pcoa.16s.point = data.frame(pcoa.16s$vectors[, 1:2],
                              location = sapply(
                                rownames(pcoa.16s$vectors),
                                function(x){unlist(strsplit(x, "\\."))[1]}))
  xylab = paste0("PCo", 1:2, " [",
                 round(pcoa.16s$values$Relative_eig[1:2]*100, 2), "%]")

  p = ggplot() +
    geom_point(data = pcoa.16s.point,
               aes_string(x = "Axis.1", y = "Axis.2",
                          color = "location"),
               size = 2) +
    geom_text_repel(data = pcoa.16s.point,
                    aes_string(x = "Axis.1", y = "Axis.2",
                               label = "rownames(pcoa.16s.point)")) +
    #stat_ellipse(data = pcoa.16s.point,
    #             aes_string(x = "Axis.1", y = "Axis.2",
    #                        color = "location"),
    #             level = 0.95, show.legend = FALSE, size = 1) +
    labs(title = "", x = xylab[1], y = xylab[2]) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))

  return(p)
}


load__pca.KO <- function(div.table) {
  pca.16s = prcomp(div.table)
  pca.16s.importance = summary(pca.16s)$importance

  plot(ceiling(1:length(pca.16s.importance)/3),
       pca.16s.importance,
       col=1:3, pch=19)
  legend('topright', c(rownames(pca.16s.importance)), col=1:3, pch=19)

  return(pca.16s)
}
report__pca <- function(pca.16s, PCx = 1, PCy = 2,
                        annot.point = data.frame(),
                        annot.arrow = data.frame(),
                        .scale = 5) {
  .id = rownames(pca.16s$x)
  pca.points = data.frame(
    x = pca.16s$x[, PCx],
    y = pca.16s$x[, PCy],
    .id = .id,
    annot.point[.id,]
  )
  pca.points$area = samples.log[rownames(pca.points), "area"]

  .id = rownames(pca.16s$rotation)
  if (! "col" %in% colnames(annot.arrow)) {
    annot.arrow[.id, "col"] = .id
  }
  pca.arrow = data.frame(
    x = (pca.16s$rotation[, PCx] * sum(abs(pca.16s$x)) /
           sum(abs(pca.16s$rotation)) * .scale),
    y = (pca.16s$rotation[, PCy] * sum(abs(pca.16s$x)) /
           sum(abs(pca.16s$rotation)) * .scale),
    .id = .id,
    annot.arrow[.id,]
  )

  p = ggplot() +
    geom_point(data = pca.points,
               aes_string(x = "x", y = "y")) +  # "occur_group"
    geom_segment(data = pca.arrow,
                 aes_string(x = "0", y = "0", xend = "x", yend = "y",
                            col = "col"),
                 arrow = arrow(angle = 22.5, length = unit(2, "mm"),
                               type = "closed")) +
    geom_text_repel(data = pca.points,
                    aes_string(x = "x", y = "y", label = ".id")) +
    geom_text_repel(data = pca.arrow,
                    aes_string(x = "x", y = "y", label = ".id")) +
    labs(title = "", x = paste0("PC", PCx), y = paste0("PC", PCy)) +
    theme_classic(base_line_size = 0)

  return(p)
}


report__enrich <- function(div.table) {
  warning("cannot run")
  dds = DESeqDataSetFromMatrix(
    countData = t(div.table),
    colData = data.frame(row.names = rownames(div.table),
                         area = samples.log[rownames(div.table), "area"]),
    design = formula("~area"))

  dds2 <- DESeq(dds[rowSums(counts(dds)) > 1,])
  res_0.05 = results(dds2, alpha = 0.05)
  res_0.05 = res_0.05[order(res_0.05$padj),]
  resdata_0.05 = merge(as.data.frame(res_0.05),
                       as.data.frame(counts(dds2, normalized = TRUE)),
                       by = "row.names", sort = FALSE)
  res.diff_0.05 = subset(resdata_0.05, padj < 0.05 & abs(log2FoldChange) > 1)

  res.diff_0.05$significance = ifelse(
    res.diff_0.05$padj < 0.0001, "p < 0.0001",
    ifelse(res.diff_0.05$padj < 0.001, "p < 0.001",
           ifelse(res.diff_0.05$padj < 0.01, "p < 0.01",
                  ifelse(res.diff_0.05$padj < 0.05, "p < 0.05", ""))))

  res.diff_0.05$Row.names = factor(res.diff_0.05$Row.names,
                                   res.diff_0.05$Row.names[order(
                                     -res.diff_0.05$log2FoldChange)])

  p = ggplot() +
    geom_point(data = res.diff_0.05,
               aes_string(x = "log2FoldChange", y = "Row.names",
                          col = "Row.names",
                          size = "significance",
                          shape = "log2FoldChange > 0")) +
    scale_size_manual(values = c("p < 0.05" = 1,
                                 "p < 0.01" = 3,
                                 "p < 0.001" = 5,
                                 "p < 0.0001" = 7)) +
    scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 16)) +
    scale_color_manual(values = order.color) +
    guides(color = "none") +
    theme_bw()
  p

  return(p)
}
#report__enrich(abundance[, c(M.e.all, M.e.any)])


report__abd.sgf.TY <- function(metabolisms = list(),
                               alternative = "two.sided") {
  KO_sample_abd = load_if_not(KO_sample_abd, load__KO_sample_abd, 1)
  entry.KO = read.csv("Analyze/pathway/entry_KO.tsv",
                      sep = "\t", header = 1)
  module.KO = unique(merge(entry.KO, module.name,
                           all.x = TRUE, by = "entry")
  )[, c("A", "B", "C", "entry", "KO")]

  m.name = names(metabolisms)[1]
  abd.sgf.TY = sapply(metabolisms[[m.name]], function(x) {
    tmp = KO_sample_abd[KO_sample_abd$KO %in%
                          module.KO$KO[module.KO[, m.name] == x], ]
    if (dim(tmp)[1] > 0) {
      tmp = aggregate(tmp$RPb, by = list(tmp$sample), sum)
      tmp.ks = ks.test(
        tmp[samples.log[tmp$Group.1, "area"] == "southern slope", "x"],
        tmp[samples.log[tmp$Group.1, "area"] == "central axis", "x"],
        alternative = alternative)
      return(tmp.ks$p.value)
    } else {
      return(1)
    }
  })
  return(abd.sgf.TY)
}


report__sgf.TY <- function() {
  de = data.frame(
    C = sapply(c(M.e.both, M.e.TY, M.e.nTY),
               function(entry) module.name[module.name$entry == entry, "C"]))
  sgf.TY = aggregate(row.names(de),
                     by = list(as.character(de$C)),
                     function(x) length(x))
  rownames(sgf.TY) = sgf.TY$Group.1; sgf.TY$Group.1 <- NULL

  report__sgf.TY.gl <- function(alternative = "two.sided") {
    abd.sgf.C = report__abd.sgf.TY(
      list(C = as.character(unique(
        sapply(M.e.both, function(entry) module.name[module.name$entry ==
                                                       entry, "C"])))),
      alternative = alternative)
    message("metabolisms ", alternative,
            " in slope with confident of p<0.02 includes:\n",
            paste(names(abd.sgf.C[abd.sgf.C < 0.02]), collapse = ", "), "\n",
            "with same p-value '", unique(abd.sgf.C[abd.sgf.C < 0.02]), "'")
    message("metabolisms with 0.02<=p<0.05 includes:", "\n",
            paste(names(abd.sgf.C[abd.sgf.C < 0.05 & abd.sgf.C >= 0.02]),
                  collapse = ", "), "\n",
            "with same p-value '", unique(abd.sgf.C[abd.sgf.C < 0.05 &
                                                      abd.sgf.C >= 0.02]), "'")
    message(length(abd.sgf.C[abd.sgf.C < 0.05]), " of ",
            length(abd.sgf.C), " metabolism is detected as signifcant different")
    
    abd.sgf.entry = report__abd.sgf.TY(list(entry = M.e.both),
                                       alternative = alternative)
    length(abd.sgf.entry)
    length(abd.sgf.entry[abd.sgf.entry < 0.05])

    return(abd.sgf.entry[abd.sgf.entry < 0.05])
  }

  for (entries in list(list("axis rich", names(report__sgf.TY.gl("less"))),
                       list("slope rich", names(report__sgf.TY.gl("greater"))),
                       list("slope compl", M.e.TY),
                       list("axis compl", M.e.nTY))) {
    de = data.frame(
      C = sapply(entries[[2]],
                 function(entry) module.name[module.name$entry == entry, "C"]))
    tmp = aggregate(row.names(de),
                    by = list(as.character(de$C)),
                    function(x) paste(x, collapse = ", "))
    sgf.TY[tmp$Group.1, entries[[1]]] = tmp$x
  }
  sgf.TY[is.na(sgf.TY)] = ""

  return(sgf.TY)
}


report__KO.abd.rank <- function(entry = "M00531") {
  # either entry or KO should be given.
  entry.KO = load_if_not(entry.KO)
  KO_sample_abd = load_if_not(KO_sample_abd)
  module.name = load_if_not(module.name)

  if (class(entry) == "character") {
    col.title = module.name[module.name$entry %in% c(entry), "name"]
    col.legend = entry
    col.KO = entry.KO[entry.KO$entry %in% c(entry), "KO"]
  } else if (class(entry) == "list") {
    col.title = names(entry)
    col.legend = "KO"
    col.KO = entry[[col.title]]
  }
  # add rank to KO_sample_abd
  tmp = unsplit(lapply(split(KO_sample_abd, KO_sample_abd$sample),
                       function(x) {
                         x$rank = order(order(x[, "RPb"],
                                               decreasing = TRUE))
                         return(x)}),
                KO_sample_abd$sample)
  tmp = tmp[tmp$KO != "",]
  p = ggplot() +
    geom_line(data = tmp,
              mapping = aes_string(x = "rank", y = "RPb",
                                   linetype = "samples.log[sample, ]$area",
                                   group = "sample"),
              col = "darkgray") +
    geom_point(data = tmp[tmp$KO %in%
                            col.KO,],
               mapping = aes_string(x = "rank", y = "RPb",
                                    col = "KO"),
               alpha = 0.7) +
    scale_y_log10(expand=expansion(),
                  labels=~ format(.x, scientific = TRUE)
                  %>% str_replace('1.*e\\+0', '10^') %>% parse(text = .)) +
    geom_polygon(
      data = as.data.table(tmp[tmp$KO %in% col.KO,]
                           )[, .SD[chull(rank, log10(RPb))], by = KO],
      mapping = aes_string(x = "rank", y = "RPb",
                           col = "KO"),
      fill = NA) +
    guides(col = guide_legend(title = col.legend),
           linetype = guide_legend(title = "area")) +
    labs(title = col.title) +
    #facet_grid(.~samples.log[sample, ]$area) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
  return(p)
}


pathview.path <- function(ko,
                          ko.data = as.matrix(
                            data.frame(genome1 = rep(1, 3),
                                       row.names = c("K00001", "K00004",
                                                     "K00008"))),
                          path = "./",
                          out.suffix = ""
) {
  last_path = getwd()
  setwd(path)
  p = pathview(gene.data = ko.data, species = 'ko',
               pathway.id = ko,
               out.suffix = out.suffix,
               kegg.native = T, kegg.dir = '~/Data/Database2/KEGG/pathview')
  setwd(last_path)
  return(p)
}
