###
#* @Date: 2021-09-20 17:32:57
#* @LastEditors: Hwrn
#* @LastEditTime: 2021-09-20 20:48:59
#* @FilePath: /2021_09-MT10kSW/Analyze/Figs/03_load.r
#* @Description:
###


comp_cont = list(c(50, 10), c(75, 10), c(90, 10),
                 c(90, 5), c(95, 5))


load__Wtdb <- function(usePrimary = FALSE) {
  if (usePrimary) {
    Wtdb = read.csv("Analyze/drep/Wtdb_1.csv", as.is = TRUE)
    Wtdb$cluster = as.character(Wtdb$primary_cluster)
    Wtdb$primary_cluster <- NULL
  } else {
    Wtdb = read.csv("Analyze/drep/Wtdb.csv", as.is = TRUE)
  }
  Wtdb$taxonomy[Wtdb$taxonomy == ""] = paste(
    "d", "p", "c", "o", "f", "g", "s__", sep = "__;")
  Wtdb$name = sapply(Wtdb$taxonomy,
                     function(x) {
                       s = strsplit(x, ';.__')[[1]]
                       i = length(s)
                       while (s[i] == "") {i = i - 1}
                       s[i]})
  Wtdb$annot_level = factor(
    sapply(Wtdb$taxonomy, function(x) {ifelse(
      taxon.split(x, 1) == "", "root", ifelse(
        taxon.split(x, 2) == "", "domain", ifelse(
          taxon.split(x, 3) == "", "phylum", ifelse(
            taxon.split(x, 4) == "", "class", ifelse(
              taxon.split(x, 5) == "", "order", ifelse(
                taxon.split(x, 6) == "", "family", ifelse(
                  taxon.split(x, 7) == "", "genus", "species")))))))}),
    levels = c("species", "genus", "family", "order",
               "class", "phylum", "domain", "root"))

  Wtdb$comp_cont = "low quality"
  for (i in comp_cont) {
    Wtdb$comp_cont[Wtdb$completeness >= i[1] & Wtdb$contamination <= i[2]
    ] = paste0(">", i[1], ", <", i[2])}
  Wtdb$comp_cont = factor(Wtdb$comp_cont, levels = c(
    "low quality",
    sapply(comp_cont, function(i) paste0(">", i[1], ", <", i[2]))))
  Wtdb$dpc = taxon.split(Wtdb$taxonomy, 1, 3)

  return(Wtdb)
}
Wtdb = load_if_not(Wtdb)


load__Stdb <- function() {
  Stdb = read.csv("Analyze/drep/Stdb.csv")
  Stdb$taxonomy = sapply(Stdb$secondary_cluster, function(x) {
    Wtdb[Wtdb$cluster == x, ]$taxonomy
  })
  Stdb$comp_cont = "low quality"
  for (i in comp_cont) {
    Stdb$comp_cont[Stdb$Completeness >= i[1] & Stdb$Contamination <= i[2]
    ] = paste0(">", i[1], ", <", i[2])}
  Stdb$comp_cont = factor(Stdb$comp_cont, levels = c(
    "low quality",
    sapply(comp_cont, function(i) paste0(">", i[1], ", <", i[2]))))

  return(Stdb)
}
Stdb = load_if_not(Stdb)  #, .force.reload = TRUE)

