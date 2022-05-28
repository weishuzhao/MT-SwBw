###
#* @Date: 2021-09-20 17:32:57
#* @LastEditors: Hwrn
#* @LastEditTime: 2021-11-10 17:54:25
#* @FilePath: /2021_09-MT10kSW/Analyze/Figs/04_draw.r
#* @Description:
###
source("Analyze/Figs/04_load.r")


annot.C = report__abd.heat()[["annot.C"]]
dev.off()
png("Analyze/Figs/03A.png", units = "in", width = 11, height = 6,res = 300)
print()
dev.off()


abundance[, M.e.both]
pheatmap(annot_for_heat(abundance[, c(M.e.both, "sample")], module.name,
                        function(x){sum(x)/length(x)}),
         cluster_rows = FALSE, cutree_cols = 2,
         annotation_row = annot.C,
         scale = 'none')


report__pcoa.KO(abundance[, c(M.e.both, M.e.nTY, M.e.TY)])


annot.arrow = data.frame(
  row.names = module.name$entry,
  col = module.name$entry %in% M.e.both,  # module.name$C,  #
  name = module.name$entry
)
report__pca(load__pca.KO(abundance[, c(M.e.both, M.e.nTY, M.e.TY)]), 1, 2,
            annot.arrow = annot.arrow, .scale = 10)


#> KO.mono
# [1] "K01870" "K02904" "K02528" "K01872" "K02892" "K02961" "K02965" "K02878"
# [9] "K03177" "K02874" "K02926" "K02890" "K01889" "K03076" "K02881" "K02994"
#[17] "K02933" "K02982" "K01890" "K03703" "K02906" "K02948" "K03106" "K02519"
#[25] "K02600" "K02996" "K02992" "K02871" "K01409" "K02946" "K00927" "K02967"
#[33] "K02956" "K02864" "K02863" "K02952" "K02950"

message("To check whether a KO family is rare or not, we should first ",
        "understand the abundance of all markers in the sample")
report__KO.abd.rank(entry = "M00531")
report__KO.abd.rank(list(mono.marker = KO.mono))


dev.off()
pdf("Analyze/pathway/entries.pdf", width = 12, height = 8)
for (entry in unlist(list(module.name[module.name$entry %in%
                                      c(M.e.both), "entry"],
                          module.name[module.name$entry %in%
                                      c(M.e.TY), "entry"],
                          module.name[module.name$entry %in%
                                      c(M.e.nTY), "entry"]))) {
  print(report__KO.abd.rank(entry = entry))
}
dev.off()


message("The total abundance of gene in the complete pathway is comparable, ",
        "while those in different pathway is incomparable. ",
        "We next first try to compare genes involved in the same pathway ",
        "in either TY or non-TY. ",
        "Warning: some gene with low abundance (by error annotation or play ",
        "minar role) may influent our results")
sgf.TY = report__sgf.TY()
write.table(sgf.TY, file = "Analyze/binannot/sgf.tsv",
            sep = "\t", quote = FALSE)

# use scripts in `Analyze/binannot/README.md#map-entry-to-pathway`
#  to re-product
tmp = read.csv("tmp.tsv", sep = "\t", header = 1, row.names = 1)
cover = apply(tmp, 1, function(x) sum(x) == 0)
tmp = tmp[!rownames(tmp) %in% names(cover[cover]),]
for (map in c('map01100',  # Metabolic pathways
              'map01110',  # Biosynthesis of secondary metabolites
              'map01120',  # Microbial metabolism in diverse environments
              'map01240',  # Biosynthesis of cofactors
              'map01230',  # Biosynthesis of amino acids
              'map01220',  # Degradation of aromatic compounds
              'map01212',  # Fatty acid metabolism
              'map01210',  # 2-Oxocarboxylic acid metabolism
              'map01200',  # Carbon metabolism
              'map00660', 'map00053',  # cover back
              '')) {
  if (!map %in% colnames(tmp))
    next
  # check if all entries of given `map` is covered by other maps
  cover = apply(tmp[tmp[map] > 0, !colnames(tmp) == map],
                1, function(x) sum(x) > 0)
  cover
  if (all(cover)) {
    tmp = tmp[!colnames(tmp) == map]
    print(paste(map, "can be safely removed"))
  } else {
    print(paste(map, "cannot be safely removed"))
  }
}
for (map in c('map00190',  # cover 光合磷酸化 (map00195)
              'map00250',  # M00027  	GABA 分流, 亦在丁酸代谢中出现
              'map00541',  # M00793  	dTDP-L-鼠李糖生物合成, 亦可能在链霉素和聚酮合成通路中
              'map00620', 'map00362',
              'map00630',  # M00741 丙酰 CoA 代谢, 亦在 map00640 和 map00280
              '')) {
  if (!map %in% colnames(tmp))
    next
  tmp[which(tmp[map] == 0), -which(colnames(tmp) == map)]
  # find those maps that already covered by this `map`
  # that is, 仅保留包含未被 `map` 覆盖的 entry 对应的其他 maps
  cover = apply(tmp[which(tmp[map] == 0), -which(colnames(tmp) == map)],
                2, function(x) sum(x) == 0)
  if (any(cover)) {
    #pheatmap(t(tmp[c(map, names(cover[cover]))]))
    print(paste(paste(names(cover[cover]), collapse = ","), "can be removed",
                "as covered by", map))
    tmp = tmp[!colnames(tmp) %in% names(cover[cover])]
  }
}
pheatmap(t(tmp[apply(tmp, 1, sum) > 1, ]),
         annotation_row = data.frame(
           uniq = apply(tmp[apply(tmp, 1, sum) == 1,],
                        2, function(x) ifelse(sum(x), "uniq", NA))))
paste0("'", paste(sort(colnames(tmp)), collapse = "', '"), "'")


tab = reshape2::acast(KO_sample_abd, formula("KO~sample"),
                      value.var = "RPb", fun.aggregate = sum)
for (path in colnames(tmp)) {
  p = pathview.path(gsub("^map(\\d{5})$", "\\1", path), tab, 
                    path = "Analyze/binannot/pathview/")
}
