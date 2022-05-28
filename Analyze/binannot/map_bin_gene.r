###
#* @Date: 2021-11-17 10:45:14
#* @LastEditors: Hwrn
#* @LastEditTime: 2022-01-18 21:27:18
#* @FilePath: /2021_12-MT10kSW/data/2021_09-MT10kSW/Analyze/binannot/map_bin_gene.r
#* @Description:
###
source("Analyze/Figs/03_func.r")
source("Analyze/Figs/03_load.r")
source("Analyze/Figs/04_load.r")

sample_name = "TY.040"
load__bin_gene = function() {
  bin_gene = data.frame()
  for (i in rownames(samples.log)) {
    tmp = read.csv(paste0("Analyze/binannot/", i, "-PF_bin.tsv"),
                   sep = "\t", header = 1, row.names = 1) /
      mrk_stat[i, "value"]
    bin_gene[rownames(tmp), colnames(tmp)] = tmp

  }
  colnames(bin_gene) = gsub(".megahit", "-megahit", colnames(bin_gene))
  bin_gene[is.na(bin_gene)] = 0

  return(bin_gene)
}
bin_gene = load_if_not(bin_gene)  #, .force.reload = TRUE)


report__pf_bin <- function(pf_names = c("PF00189.15"),  # rpS3
                           cutoff = 1) {
  #taxon.sort = data.frame(cluster = Wtdb$cluster,
  #                        taxonomy = Wtdb$taxonomy)
  #taxon.sort$cluster = factor(taxon.sort$cluster,
  #                            levels = taxon.sort$cluster[
  #                              order(taxon.sort$taxonomy)])
  tab = reshape2::acast(KO_sample_abd, formula("KO~sample"),
                        value.var = "RPb", fun.aggregate = sum) * 100
  tmp = data.frame()
  for (pf_name in pf_names) {  # pf_name = pf_names[1]
    if (!pf_name %in% rownames(bin_gene)) {
      warning(pf_name, " isn't found in any bin in any sample!")
      next
    }
    div.Stdb = data.frame(
      Stdb,
      protein = pf_name,
      annot.percent = c(t(bin_gene[pf_name,
                                   gsub("^(\\d)", "X\\1",
                                        gsub("-|;| ", ".",
                                             Stdb$genome))]) * 100))

    if (pf_name %in% rownames(tab)) {
      div.others = data.frame(
        genome = "others", score = 0, secondary_cluster = "others",
        Completeness = 0, Contamination = 100, Strain.heterogeneity = 100,
        length = 0, sample = colnames(tab), AvgDepth = 0, taxonomy = "others",
        comp_cont = "low quality", protein = pf_name,
        annot.percent = tab[pf_name, ] - reshape2::acast(
          div.Stdb, formula("protein~sample"),
          value.var = "annot.percent", fun.aggregate = sum)[,colnames(tab)])
    } else {div.others = data.frame()}

    tmp = rbind(tmp, div.Stdb, div.others)
  }
  tmp$text = ifelse(taxon.split(tmp$taxonomy, 6) == "",
                    as.character(tmp$secondary_cluster),
                    as.character(taxon.split(tmp$taxonomy, 6)))
  tmp$text = as.character(tmp$secondary_cluster)
  tmp$label.text = ifelse(tmp$annot.percent >= cutoff,
                          #as.character(tmp$secondary_cluster),
                          as.character(tmp$text),
                          NA)

  tmp$index = taxon.split(tmp$taxonomy, 1, 6)

  p = ggplot(data = tmp) +
    geom_bar(aes_string(x = "sample", y = "annot.percent",
                        fill = "index",
                        alpha = "comp_cont %in% c('>90, <5', '>95, <5')"),
             position = position_stack(reverse = TRUE),
             stat = "identity", col = "black") +
    #scale_fill_manual(values = order.color) +
    #scale_fill_manual(values = brewer_pal(palette = "BrBG")(17)) +
    scale_y_continuous(breaks = seq(0, 100, 10)) +
    labs(title = "", x = "sample", y = "precent %") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    facet_grid(formula(".~protein"))

  p$data["label.y"] = (ggplot_build(p)$data[[1]]["y"] -
                         p$data["annot.percent"] / 2)

  return(p)
}
p = report__pf_bin(
  c(unlist(strsplit(paste(
    "K00372 K00360 K00366",
    ""), split = " "))
  )) +
  guides(fill = "none", alpha = "none") +
  geom_text(aes_string(x = "sample", y = "label.y",
                       label = "label.text"),
            position = "identity")


sort(rownames(bin_gene), decreasing = TRUE)[1:56]


for (pf_names in list(
  #"K00372 K00360 K00367",
  #"K00370 K00371 K00374",
  #"K02567 K02568",
  #"K00362 K00363",
  #"K00366 K03385 K15876",
  #"K10944 K10945 K10946 K10535",
  #"K00821 K00611 K01940 PF00189.15",  #"K00611 K01940 K01755 PF00189.15",
  #"K00931 K00147 K00286",
  #"K00260 K00261 K00262",
  #"K13788 K00625 K00925",
  "K01895 K01907",
  #"K00486 K00452 K03392 K10217",
  #"K03381 K01856 K03464", "K01055 K14727",
  #"K11258 K09011 K01702",
  #"K00166 K00167 K11381 K09699",
  #"K01758 K01697",
  #"K01825 K01782 K00632",
  #"K00249 K06445 K01692 K07516",
  #"K01595 K00024",
  #"K01647 K05942 K01681 K01682",
  #"K00031 K00030",
  #"K00036 K01057 K07404",
  #"K19243 K01810 K15916",
  #"K00134 K00927 K00131",
  #"K03737 K00169 K00170 K00171 K00172",
  #"K00720 K07553 K04718 K01634",
  #"K11780 K11781 K14941 K11212 K12234",
  #"K17222 K17223 K17224", "K17225 K22622 K17226 K17227",
  #"K02111 K02112 K02113 K02114", "K02115 K02108 K02109 K02110",
  #"K02117 K02118 K02119 K02120", "K02121 K02122 K02107 K02123 K02124",
  #"K02172 K02171 K17836", "K07644 K07665",
  #"K18301 K18302 K18303 K18139",
  #sort(rownames(bin_gene), decreasing = TRUE)[51:53],  # marker
  ""
)) {
  if (pf_names == "") next

  #pf_names = "K01548 K02119 K02128"
  p = report__pf_bin(
    c(unlist(strsplit(paste(
      pf_names,
      ""), split = " "))
    ), cutoff = 1) +
    guides(fill = "none", alpha = "none") +
    geom_text(aes_string(x = "sample", y = "label.y",
                         label = "label.text"),
              position = "identity")
  #print(p); next

  lapply(unique(p$data[p$data$comp_cont %in% c(">90, <5", ">95, <5") &
                         #p$data$protein == "K00260" &
                         p$data$annot.percent > 0,
                       c("text")]),
         function(x) list(unique(p$data[p$data$text == x,
                                        c("text", "taxonomy")]),
                          p$data[#p$data$protein == "K00260" &
                                   p$data$text == x,
                                 c("sample", "comp_cont",
                                   "annot.percent", "label.y", "protein")]))

  ggsave(plot = p,
         path = "tmp/", filename = paste0(gsub(" ", "-", pf_names), ".jpg"),
         device = "jpg", dpi = 300, width = 15, height = 9)

}


fc.genomes = Stdb[Stdb$secondary_cluster %in% ws.clusters, ]$genome
fc.genes = list(
  NO3.red.Ass    = "K00372 K00360 K00367",  # NasA, NasB, NarB
  NO3.red.Dis    = "K00370 K00371 K00374",  # NarG, NarH, NarI
  NO3.red.Nap    = "K02567 K02568",  # NapA, NapB
  NO2.red.NirABD = "K00366 K00362 K00363",  # NirA, NirB, NirD
  NO2.red.NrfAH  = "K03385 K15876",  # NrfA, NrfH
  NH3.oxd.AOA    = "K10944 K10945 K10946",  # AOA, K10535 miss
  Glu.Arg        = "K00821 K00611 K01940",  # Glu -> Arg
  Arg.M00038     = "K00486 K00452 K03392 K10217",  # M00038
  PhOH.M00568.1  = "K03381 K01856 K03464",
  PhOH.M00568.2  = "K01055 K14727"
  #sort(rownames(bin_gene), decreasing = TRUE)[51:53],  # marker
)
fc.genes = lapply(fc.genes, function(x) unlist(strsplit(x, " ")))
fc.bin_gene = bin_gene[unlist(strsplit(unlist(fc.genes), " ")),
                       as.character(fc.genomes)]
fc.bin_gene = fc.bin_gene[apply(fc.bin_gene, 1, sum) > 0,
                          apply(fc.bin_gene, 2, sum) > 0]
#fc.bin_gene[fc.bin_gene == 0] = NA
annotation_row = melt(as.data.frame(sapply(
  fc.genes, "[", i = 1:max(sapply(fc.genes, length)))), id.vars = c())
annotation_row = annotation_row[!is.na(annotation_row["value"]), ]
p = pheatmap(fc.bin_gene, display_numbers = ifelse(fc.bin_gene == 0, "0", ""),
         annotation_col = data.frame(
           species = sapply(fc.genomes, function(x)
             taxon.split(Stdb[Stdb$genome == x, "taxonomy"], 7)),
           sample = sapply(as.character(fc.genomes), function(x)
             unlist(strsplit(x, "\\."))[1]),
           row.names = fc.genomes),
         annotation_row = data.frame(
           catalyze = annotation_row[["variable"]],
           row.names = annotation_row[["value"]]),
         clustering_method = "average",
         show_colnames = FALSE,
         silent = TRUE)
svg(filename = "Analyze/binannot/fc.bin_gene.heat.svg",
    width = 10, height = 6)
p
dev.off()


# 2022-01-16 20:21:32
ws.clusters = c("126_1", "133_1", "140_1", "157_1")
fc.genomes = Stdb[Stdb$secondary_cluster %in% ws.clusters &
                    as.integer(Stdb$comp_cont) > 4, ]$genome
fc.genes = rbind(
  #data.frame(cpds = "C00064, C00169", kos = c('K01954', 'K01955', 'K01956', 'K11540', 'K11541')),
  data.frame(rct = "NO2- <=> NH3",   cpds = "C00014, C00088", kos = c('K00362', 'K00363', 'K03385')),
  data.frame(rct = "NO3- <=> NO2-" , cpds = "C00088, C00244", kos = c('K00360', 'K00372')),
  data.frame(rct = "NH3+HCOOH",      cpds = "C00014, C00488", kos = c('K01455')),
  data.frame(rct = "Gln.syn",        cpds = "C00014, C00064", kos = c('K01915', 'K23265')),
  data.frame(rct = "Gly.syn",        cpds = "C00014, C00037", kos = c('K00281', 'K00282', 'K00283', 'K00382', 'K00605', 'K02437')),
  data.frame(rct = "R-CN.deg",       cpds = "C00014, C00726", kos = c('K01501')),
  data.frame(rct = "P-CONH2.syn",    cpds = "C00014, C00169", kos = c('K00926')),
  data.frame(rct = "NO2- ==> NO",    cpds = "C00088, C00533", kos = c('K00368')),
  data.frame(rct = "NH2OH <=> NH3",  cpds = "C00014, C00192", kos = c('K10944', 'K10945', 'K10946')),
  data.frame()
)
fc.genes$name = c("nirB", "nirD", "nrfA",
                  "nasB", "nasA",
                  "formamidase",
                  "glnA, GLUL",
                  "purQ", "GLDC, gcvP", "gcvPA", "gcvPB", "DLD, lpd, pdhD", "gcvT, AMT", "gcvH, GCSH",
                  "nitrilase",
                  "arcC",
                  "nirK",
                  "pmoA-amoA", "pmoB-amoB", "pmoC-amoC")
fc.genes$desc = c("nitrite reductase (NADH) large subunit",
                  "nitrite reductase (NADH) small subunit",
                  "nitrite reductase (cytochrome c-552)",
                  "assimilatory nitrate reductase electron transfer subunit",
                  "assimilatory nitrate reductase catalytic subunit",
                  "formamidase",
                  "glutamine synthetase",
                  "phosphoribosylformylglycinamidine synthase subunit PurQ / glutaminase",
                  "glycine dehydrogenase",
                  "glycine dehydrogenase subunit 1",
                  "glycine dehydrogenase subunit 2",
                  "dihydrolipoamide dehydrogenase",
                  "aminomethyltransferase",
                  "glycine cleavage system H protein",
                  "nitrilase",
                  "carbamate kinase",
                  "nitrite reductase (NO-forming)",
                  "methane/ammonia monooxygenase subunit A",
                  "methane/ammonia monooxygenase subunit B",
                  "methane/ammonia monooxygenase subunit C")
rownames(fc.genes) = as.character(fc.genes$kos)
fc.bin_gene = bin_gene[rownames(fc.genes),
                       as.character(fc.genomes)]
fc.bin_gene = fc.bin_gene[apply(fc.bin_gene, 1, sum) > 0,
                          apply(fc.bin_gene, 2, sum) > 0]
#fc.bin_gene[fc.bin_gene > 0] = 1
#fc.bin_gene[fc.bin_gene == 0] = NA
p = pheatmap(fc.bin_gene, display_numbers = ifelse(fc.bin_gene == 0, "0", ""),
             annotation_col = data.frame(
               species = sapply(fc.genomes, function(x)
                 taxon.split(Stdb[Stdb$genome == x, "taxonomy"], 7)),
               sample = sapply(as.character(fc.genomes), function(x)
                 unlist(strsplit(x, "\\."))[1]),
               row.names = fc.genomes),
             annotation_row = fc.genes[c("rct")],
             clustering_method = "average",
             show_colnames = FALSE,
             silent = TRUE)
svg(filename = "Analyze/binannot/fc.bin_gene.heat.svg",
    width = 10, height = 6)
p
dev.off()



#sapply(ws.clusters, function(x) Wtdb$name[Wtdb$cluster == x])
sapply(as.character(ws.clusters), function(x) {
  genomes = Stdb[Stdb$secondary_cluster == x, ]$genome
  apply(bin_gene[p$tree_row$labels, as.character(genomes)], 1, sum)
})[c("K00611", "K00363", "K00452", "K03392", "K00486", "K10217", "K01940",
     "K10944", "K10945", "K10946", "K03385", "K00362", "K00372", "K00360",
     "K00821", "K14727", "K01055", "K03381", "K03464"),] > 0

report__MAGs_KOs = function(MAGs) {
  KOs = apply(bin_gene[, as.character(MAGs)], 1, sum) > 0
  KOs = KOs & sapply(names(KOs), function(x)
    length(unlist(strsplit(x, "PF|TIGR"))) == 1)
  return(paste(names(KOs[KOs]), collapse = " "))
}

bin_gene[c("K00368", "K15864"), Stdb[Stdb$secondary_cluster == "140_1", "genome"]]

for (pf_names in fc.genes) {
  if (pf_names == "") next

  #pf_names = "K01548 K02119 K02128"
  p = report__pf_bin(
    c(unlist(strsplit(paste(
      pf_names,
      ""), split = " "))
    ), cutoff = 1) +
    guides(fill = "none", alpha = "none") +
    geom_text(aes_string(x = "sample", y = "label.y",
                         label = "label.text"),
              position = "identity")
  #print(p); next

  lapply(unique(p$data[p$data$comp_cont %in% c(">90, <5", ">95, <5") &
                         #p$data$protein == "K00260" &
                         p$data$annot.percent > 0,
                       c("text")]),
         function(x) list(unique(p$data[p$data$text == x,
                                        c("text", "taxonomy")]),
                          p$data[#p$data$protein == "K00260" &
                            p$data$text == x,
                            c("sample", "comp_cont",
                              "annot.percent", "label.y", "protein")]))

  ggsave(plot = p,
         path = "tmp/", filename = paste0(gsub(" ", "-", pf_names), ".jpg"),
         device = "jpg", dpi = 300, width = 15, height = 9)

}
