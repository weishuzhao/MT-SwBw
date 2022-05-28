###
#* @Date: 2021-09-12 23:20:00
#* @LastEditors: Hwrn
#* @LastEditTime: 2021-09-15 20:23:12
#* @FilePath: /2021_09-MT10kSW/Analyze/drep/01_draw.r
#* @Description:
###


comp_cont = list(c(50, 10), c(75, 10), c(90, 10),
                 c(90, 5), c(95, 5))


load__Wtdb <- function () {
  Wtdb = read.csv("Analyze/drep/Wtdb.csv", as.is = TRUE)
  Wtdb$name = sapply(as.character(Wtdb$taxonomy),
                     function(x) {
                       s = strsplit(x, ';.__')[[1]]
                       i = length(s)
                       while (s[i] == "") i = i - 1
                       s[i]})
  Wtdb$annot_level = factor(
    sapply(Wtdb$taxonomy, function (x) {ifelse(
      taxon.split(x, 1) == "d__", "root", ifelse(
        taxon.split(x, 2) == "p__", "domain", ifelse(
          taxon.split(x, 3) == "c__", "phylum", ifelse(
            taxon.split(x, 4) == "o__", "class", ifelse(
              taxon.split(x, 5) == "f__", "order", ifelse(
                taxon.split(x, 6) == "g__", "family", ifelse(
                  taxon.split(x, 7) == "s__", "genus",
                  "species")))))))}),
    levels = c("fine classified", "species", "genus", "family",
               "order", "class", "phylum", "domain"))

  for (i in comp_cont) {
    Wtdb$comp_cont[Wtdb$contamination < i[1] & Wtdb$completeness > i[2]
                   ] = paste0(">", i[1], ", <", i[2])}
  Wtdb$comp_cont = factor(Wtdb$comp_cont, levels = sapply(
    comp_cont, function (i) paste0(">", i[1], ", <", i[2])
  ))

  return(Wtdb)
}
Wtdb = load__Wtdb()


summarize__Wtdb <- function() {
  p = ggplot(data = Wtdb) +
    geom_point(mapping = aes(x=contamination, y=completeness,
                             color=annot_level, size=score),
               alpha=0.5) +
    theme_classic(base_line_size = 0) +
    scale_x_continuous(limits = c(0, 10)) + scale_y_continuous(limits = c(50, 100))
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
s = summarize__Wtdb()
s[[1]] +
  annotate("segment", x =  0, xend = 10, y = 75, yend =  75) +
  annotate("segment", x =  0, xend = 10, y = 90, yend =  90) +
  annotate("segment", x =  5, xend =  5, y = 90, yend = 100) +
  annotate("segment", x =  0, xend =  5, y = 95, yend =  95)

s[[2]]
s[[3]]


summarize__Wtdb.sample <- function (sample.index) {
  annotation_row = data.frame(
    row.names = colnames(Wtdb)[sample.index],
    site = sapply(
      colnames(Wtdb)[sample.index],
      function (x) unlist(strsplit(x, "\\."))[1]
    ))
  p = pheatmap(t(Wtdb[apply(Wtdb[, sample.index], 1, sum) > 1, sample.index]),
               show_colnames = FALSE,
               annotation_row = annotation_row)
  
  sumsingle = apply(Wtdb[apply(Wtdb[, sample.index], 1, sum) == 1, 
                         sample.index], 2, sum)
  sumsingle = data.frame(
    sumsingle, 
    sumsingle.rate = sumsingle / apply(Wtdb[, sample.index], 2, sum)
  )

  tax = Wtdb[apply(Wtdb[, sample.index], 1, function (x) {
    (sum(x[c(1:3)]) > 0) +
      (sum(x[c(4:12)]) > 0) > 1
  }), "taxonomy"]
  
  return(list(colnames(Wtdb)[sample.index], p,
              sumsingle, tax))
}
s = summarize__Wtdb.sample(16:27)
"
> s[[4]]
d__Archaea;p__Thermoplasmatota;c__Poseidoniia_A;o__Poseidoniales;f__Thalassarchaeaceae;g__;s__
d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae;g__Nitrosopumilus;s__Nitrosopumilus sp013390905
d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Flavobacteriaceae;g__PRS1;s__PRS1 sp003709165
d__Bacteria;p__Chloroflexota;c__Dehalococcoidia;o__UBA3495;f__UBA3495;g__UBA9611;s__
d__Bacteria;p__Marinisomatota;c__Marinisomatia;o__Marinisomatales;f__TCS55;g__UBA2126;s__UBA2126 sp004124385
d__Bacteria;p__Planctomycetota;c__Phycisphaerae;o__Phycisphaerales;f__SM1A02;g__GCA-002718515;s__
d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__PS1;f__Thioglobaceae;g__DUCF01;s__
"
