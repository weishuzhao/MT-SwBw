###
#* @Date: 2021-09-16 14:12:32
#* @LastEditors: Hwrn
#* @LastEditTime: 2021-11-17 15:00:55
#* @FilePath: /2021_09-MT10kSW/Analyze/Figs/02_load.r
#* @Description:
#   vars loading
###
source("Analyze/Figs/00_load.r")
#RELOAD = TRUE
RELOAD = FALSE


load__taxon.order.1 <- function() {
  df1 = data.frame()
  for (i in samples.log$X.Sample) {
    samplen = i
    df1 = rbind(df1, data.frame(
      read.csv(paste0("Analyze/alphadiv/phyloflash-", samplen, "-megahit.csv"),
               sep = ",", col.names = c("taxon", "reads"), as.is = T),
      sample = samplen))
  }
  taxon = sort(unique(df1$taxon))

  taxon.order = data.frame(
    domain = taxon.split(taxon, 1),
    phylum = taxon.split(taxon, 2),
    class  = taxon.split(taxon, 3),
    order  = taxon.split(taxon, 4),
    name.16s = taxon.split(taxon, 1, 4),
    dpco = taxon.split(taxon, 1, 4)
  )
  taxon.order$name.kaiju = as.character(taxon.split(taxon, 3))
  taxon.order = taxon.order[
    sapply(taxon.order$domain,
           function(x) strsplit(as.character(x), "")[[1]][1]) != "E", ]
  taxon.order = taxon.order[!duplicated(taxon.order$dpco), ]
  rownames(taxon.order) = taxon.order$dpco

  dpc = unique(taxon.split(taxon, 1, 3))
  class.dup = unique(taxon.split(dpc, 3)[duplicated(taxon.split(dpc, 3))])

  taxon.order$name = ifelse(  # name by mixed order and class
    taxon.order$class %in% c(
    #  "Gammaproteobacteria", "Alphaproteobacteria"
    ),
    taxon.order$order, as.character(taxon.order$class)
  )
  taxon.order$name = ifelse(
    sapply(as.character(taxon.order$class),  # unknown class
           function(x) strsplit(x, "")[[1]][1]) == "("
    | taxon.order$class %in% class.dup,
    "unassigned", as.character(taxon.order$name))

  return(taxon.order)
}
update__taxon.order <- function(taxon.order) {
  ## fix bug
  # https://doi.org/10.1099/ijsem.0.004213
  taxon.order[taxon.order$name %in% c("Bdellovibrionia"),
              c("name")] = c("Oligoflexia")
  taxon.order[taxon.order$name %in%
                c("Desulfuromonadia","Myxococcia", "Polyangia"),
              c("name")] = c("Deltaproteobacteria")

  taxon.order["Bacteria;Planctomycetota;Planctomycetes;kaiju",
              c("name.kaiju", "name")] = c(
                "Planctomycetia", "Planctomycetes"
              )  # https://www.ncbi.nlm.nih.gov/pubmed/26654112
  taxon.order["Bacteria;Proteobacteria;Gammaproteobacteria;kaiju",
              c("name.kaiju", "name")] = c(
                "Betaproteobacteria", "Gammaproteobacteria"
              )  # https://www.arb-silva.de/documentation/faqs/

  # https://doi.org/10.3389/fmicb.2017.00682
  taxon.order[taxon.order$name %in% c("Campylobacteria"),
              c("name")] = c("Epsilonproteobacteria")
  # https://doi.org/10.1007/s10482-015-0532-1
  taxon.order[taxon.order$name %in% c("Chlamydiae"),
              c("name")] = c("Chlamydiia")
  taxon.order[taxon.order$name == "Acidobacteriae",
              c("name.kaiju")] = c("Acidobacteriia"
              )  # Acidobacteriales
  taxon.order[taxon.order$name == "Babeliae",
              c("name.kaiju")] = c("Candidatus Babeliae"
              )  # Candidatus
  taxon.order[taxon.order$name == "Nitrospiria",
              c("name.kaiju")] = c("Candidatus Babeliae"
              )  # Nitrospirales

  # https://doi.org/10.3389/fmicb.2019.02083
  taxon.order["Bacteria;Bacteroidota;Rhodothermia;Balneolales",
              c("name.kaiju", "name")] = c("Balneolia", "Rhodothermia")
  for (name.kaiju in c("Cytophagia", "Flavobacteriia", "Sphingobacteriia")) {
    taxon.order[paste("Bacteria;Bacteroidota",
                      "Bacteroidia", name.kaiju, sep = ";"),
                c("name.kaiju", "name")] = c(name.kaiju, "Bacteroidia")
  }
  # https://doi.org/10.1007/s00792-019-01150-3
  for (name.kaiju in c("Acidimicrobiia", "Actinobacteria", "Coriobacteriia",
                       "Nitriliruptoria", "Rubrobacteria", "Thermoleophilia")) {
    taxon.order[paste("Bacteria;Actinobacteriota",
                      "Actinobacteria", name.kaiju, sep = ";"),
                c("name.kaiju", "name")] = c(name.kaiju, "Actinobacteria")
  }

  return(taxon.order)
}
load__taxon.order <- function() {
  return(update__taxon.order(load__taxon.order.1()))
}
taxon.order = load_if_not(taxon.order, .force.reload = RELOAD)


load__div.16s <- function() {
  df1 = data.frame()
  for (i in samples.log$X.Sample) {
    samplen = i
    df1 = rbind(df1, data.frame(
      read.csv(paste0("Analyze/alphadiv/phyloflash-", samplen, "-megahit.csv"),
               sep = ",", col.names = c("taxon", "reads"), as.is = T),
      sample = samplen))
  }
  df1$dpco = taxon.split(df1$taxon, 1, 4)
  # filter bad taxon  # Eukaryota
  df1 = df1[sapply(df1$dpco, function(x) strsplit(as.character(x), ""
                                                  )[[1]][1]) != "E", ]
  df1$name = taxon.order$name[
    sapply(as.character(df1$dpco),
           function(x) which(x == taxon.order$name.16s))]

  df1$name[sapply(df1$name, function(x) "Marinimicrobia" %in% unlist(strsplit(x, " ")))
           ] = "Marinimicrobia"

  div.16s = df1[df1$name %in% taxon.div.top(df1, "name", topi = 0), ]
  return(div.16s)
}
div.16s = load_if_not(div.16s, .force.reload = RELOAD)


load__div.kaiju <- function() {
  df2 = read.csv("Analyze/alphadiv/kaiju.order.tsv", sep = "\t",
                 col.names = c("sample", "precent", "reads",
                               "taxon_id", "taxon_name"),
                 as.is = TRUE)
  df2$sample = sapply(df2$sample, function(x) {
    unlist(strsplit(unlist(strsplit(x, "kaiju-"))[2], "-megahit"))[1]
  })
  df2 = df2[!(df2$taxon_name %in% c("Viruses", "unclassified")), ]
  df2 = df2[sapply(df2$taxon_name,
                   function(x) strsplit(x, "")[[1]][1]) != "E", ]

  # NA in kaiju.taxon_name
  kaiju.name.na = sort(unique(df2$taxon_name[sapply(
    as.character(df2$taxon_name),
    function(x)length(unlist(strsplit(x, "NA"))) > 1)]))

  df2$taxon_name[df2$taxon_name == "Archaea;Candidatus Altiarchaeota;NA;Candidatus Altiarchaeales;"
                 ] = "Archaea;Altiarchaeota;Altiarchaeia;Altiarchaeales;"
  df2$taxon_name[df2$taxon_name == "Archaea;Nanoarchaeota;NA;Nanoarchaeales;"
                 ] = "Archaea;Nanoarchaeota;Nanoarchaeia;Nanoarchaeales;"
  # https://doi.org/10.1007/s10482-010-9488-3
  df2$taxon_name[df2$taxon_name == "Archaea;Thaumarchaeota;NA;Cenarchaeales;"
                 ] = "Archaea;Crenarchaeota;Marine Benthic Group A;Cenarchaeales;"
  df2$taxon_name[df2$taxon_name == "Archaea;Thaumarchaeota;NA;Nitrosopumilales;"
                 ] = "Archaea;Crenarchaeota;Nitrososphaeria;Nitrosopumilales;"
  df2$taxon_name[df2$taxon_name == "Bacteria;Bacteroidetes;NA;Bacteroidetes Order II. Incertae sedis;"
                 ] = "cannot be assigned to a (non-viral) class"
  df2$taxon_name[df2$taxon_name == "Bacteria;Candidatus Eremiobacteraeota;NA;Candidatus Eremiobacterales;"
                 ] = "cannot be assigned to a (non-viral) class"
  df2$taxon_name[df2$taxon_name == "Bacteria;Candidatus Melainabacteria;NA;Candidatus Gastranaerophilales;"
                 ] = "Bacteria;Cyanobacteria;Vampirivibrionia;Gastranaerophilales;"
  df2$taxon_name[df2$taxon_name %in% c("Bacteria;Cyanobacteria;NA;Chroococcales;",
                                       "Bacteria;Cyanobacteria;NA;Chroococcidiopsidales;",
                                       "Bacteria;Cyanobacteria;NA;Gloeoemargaritales;",
                                       "Bacteria;Cyanobacteria;NA;Nostocales;",
                                       "Bacteria;Cyanobacteria;NA;Oscillatoriales;",
                                       "Bacteria;Cyanobacteria;NA;Pleurocapsales;",
                                       "Bacteria;Cyanobacteria;NA;Spirulinales;")
                 ] = "cannot be assigned to a (non-viral) class"
  df2$taxon_name[df2$taxon_name == "Bacteria;Cyanobacteria;NA;Synechococcales;"
                 ] = "Bacteria;Cyanobacteria;Cyanobacteriia;Synechococcales;"
  df2$taxon_name[df2$taxon_name == "Bacteria;NA;Candidatus Babeliae;Candidatus Babeliales;"
                 ] = "cannot be assigned to a (non-viral) class"
  df2$taxon_name[df2$taxon_name == "Bacteria;NA;NA;Haloplasmatales;"
                 ] = "Bacteria;Tenericutes;Mollicutes;Haloplasmatales"

  df2$name = ifelse(df2$taxon_name ==
                       "cannot be assigned to a (non-viral) class",
                    "unassigned", as.character(taxon.split(df2$taxon_name, 3)))
  df2$name = sapply(df2$name, function(x) ifelse(
    x %in% taxon.order$name.kaiju,
    taxon.order$name[x == taxon.order$name.kaiju], x
  ))

  div.kaiju = df2[df2$name %in% taxon.div.top(df2, "name", topi = 0), ]
  return(div.kaiju)
}
div.kaiju = load_if_not(div.kaiju, .force.reload = RELOAD)


compare__div.16s.kaiju <- function(topi = 10) {
  name.16s = sort(unique(div.16s$name))
  name.kaiju = sort(unique(div.kaiju$name))
  name.16s.top = taxon.div.top(div.16s, "name", topi)
  name.kaiju.top = taxon.div.top(div.kaiju, "name", topi)

  rt1 = list("these are *name.16s.top* not in *name.kaiju*",
             name.16s.top[!name.16s.top %in% name.kaiju])
  rt2 = list("these are *name.kaiju.top* not in *name.16s*",
             name.kaiju.top[!name.kaiju.top %in% name.16s])

  rt3 = list("the problemable name of SILVA",
             name.16s[!name.16s %in% name.kaiju])
  rt4 = list("the problemable name of NCBI",
             name.kaiju[!name.kaiju %in% name.16s])

  return(list(rt1, rt2, rt3, rt4))
}
#compare__div.16s.kaiju(7)
#taxon.order[taxon.order$name == "Bdellovibrionia",]
#"Oligoflexia" %in% sort(unique(div.kaiju$name))


load__name.top <- function(taxon = "name", topi = 7) {
  name.16s = unique(div.16s[, taxon])
  name.kaiju = unique(div.kaiju[, taxon])

  name.16s.top = taxon.div.top(div.16s, "name", topi = topi)
  name.kaiju.top = taxon.div.top(div.kaiju, "name", topi = topi)

  name.top = taxon.factor.sort(unique(c(
    name.16s.top, name.kaiju.top[name.kaiju.top %in% name.16s], "others")))

  return(name.top)
}
load__name.top <- function(taxon = "name", topi = 7) {
  name.16s = unique(div.16s[, taxon])
  name.kaiju = unique(div.kaiju[, taxon])
  
  name.16s.top = taxon.div.top(div.16s, "name", topi = topi)
  name.kaiju.top = taxon.div.top(div.kaiju, "name", topi = topi)
  
  name.top = taxon.factor.sort(unique(c(
    name.16s.top, #name.kaiju.top[name.kaiju.top %in% name.16s],
    "others")))
  
  return(name.top)
}
name.top = load_if_not(name.top, .force.reload = RELOAD)


trans__div2table <- function(div.any,
                             .formula = formula("sample~name"),
                             .value.var = "reads") {
  # sample~name, value.var = "reads
  div.table = reshape2::dcast(div.any, formula = .formula,
                              value.var = .value.var, fun.aggregate = sum)
  rownames(div.table) = div.table$sample
  div.table[, c("sample", "others", "unassigned")] <- NULL

  return(div.table)
}


report__pcoa.16s <- function() {
  div.16s.table = trans__div2table(load_if_not(div.16s, load__div.16s, 1))

  pcoa.16s = pcoa(vegdist(div.16s.table, method = "bray"))
  pcoa.16s.point = data.frame(-pcoa.16s$vectors[, 1:2],
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
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

  return(p)
}


report__enrich.16s <- function(order.color = vector()) {
  div.16s.table = trans__div2table(load_if_not(div.16s, .force.reload = TRUE))

  dds = DESeqDataSetFromMatrix(
    countData = t(div.16s.table),
    colData = data.frame(row.names = rownames(div.16s.table),
                         location = ifelse(
                           sapply(
                             rownames(div.16s.table),
                             function(x){
                               unlist(strsplit(x, "\\."))[1]} == "TY"),
                           "slope", "axis"
                         )),
    design = formula("~location"))

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


print(c("taxon.order",
        "div.16s", "div.kaiju",
        "name.top",
        "clust.16s"))
