###
#* @Date: 2021-09-20 17:32:57
#* @LastEditors: Hwrn
#* @LastEditTime: 2021-11-08 21:00:12
#* @FilePath: /2021_09-MT10kSW/Analyze/Figs/04_load.r
#* @Description:
###
source("Analyze/Figs/00_load.r")
source("Analyze/Figs/04_func.r")


load__mrk_stat <- function() {
  mrk_stat = read.table("Analyze/annot/marker/mkr_stat.tsv")
  mrk_stat$value = mrk_stat$mrk.med

  return(mrk_stat[rownames(samples.log), ])
}
mrk_stat = load_if_not(mrk_stat)


load__module.name <- function() {
  module.name = read.csv("Analyze/pathway/module_name.tsv",
                         sep = "\t", header = 1)
  return(module.name)
}
module.name = load_if_not(module.name)


load__KO_sample_abd <- function() {
  # already corrected by marker gene abundance
  KO_sample_RBp = read.csv("Analyze/pathway/KO_sample_RPb.tsv",
                           sep = "\t", header = 1)
  KO_sample_RBp.KO = KO_sample_RBp$KO; KO_sample_RBp$KO <- NULL

  KO_sample_RBp = as.data.frame(t(t(as.matrix(KO_sample_RBp)) /
                                    mrk_stat[colnames(KO_sample_RBp),
                                             "value"]))
  KO_sample_RBp$KO = KO_sample_RBp.KO
  KO_sample_abd = reshape2::melt(KO_sample_RBp,
                                 id.vars = "KO",
                                 variable.name = "sample",
                                 value.name = "RPb")
  KO_sample_abd = KO_sample_abd[KO_sample_abd$RPb > 0, ]

  return(KO_sample_abd)
}
KO_sample_abd = load_if_not(KO_sample_abd)


load__abd.sgf <- function() {
  KO_sample_abd = load_if_not(KO_sample_abd, load__KO_sample_abd)
  entry.KO = read.csv("Analyze/pathway/entry_KO.tsv",
                      sep = "\t", header = 1)
  module.KO = unique(merge(entry.KO, module.name,
                           all.x = TRUE, by = "entry")
                     )[, c("A", "B", "C", "entry", "KO")]

  abd.sgf = data.frame()
  for (metabolism in unique(module.name$C)) {  # metabolism = "Carbon fixation"
    tmp = KO_sample_abd[KO_sample_abd$KO %in%
                          module.KO$KO[module.KO$C == metabolism], ]
    if (dim(tmp)[1] > 0) {
      tmp = as.data.frame(
        compare_means(formula = formula("RPb~sample"), tmp,
                      method = "wilcox.test", ref.group = ".all."))
      abd.sgf[tmp$group2, metabolism] = tmp$p.signif
    } else {
      abd.sgf[, metabolism] = NA
    }
  }
  abd.sgf[abd.sgf == "ns" | is.na(abd.sgf)] = ""

  return(abd.sgf)
}
abd.sgf = load_if_not(abd.sgf)


load__ex_ab <- function() {
  sample_pathway = read.csv("Analyze/pathway/abd_cpl.tsv",
                            sep = "\t", header = 1, row.names = 1)
  sample_pathway = sample_pathway / mrk_stat[rownames(sample_pathway), "value"]

  existance = as.data.frame(apply(sample_pathway, 2, function(x) x > 0))
  existance$sample = rownames(sample_pathway)

  abundance = as.data.frame(apply(sample_pathway, 2, abs))
  abundance$sample = rownames(sample_pathway)

  return(list(ex = existance, ab = abundance))
}
ex_ab = load_if_not(ex_ab)
existance = ex_ab$ex
abundance = ex_ab$ab


annot_for_heat <- function(df_in, module.name, FUN) {
  df = melt(df_in, id.vars = 'sample',
            variable.name = 'entry', value.name = 'value')
  df1 = merge(df, module.name, by = 'entry', all.x = T)
  df2 = aggregate(formula("value ~ sample + C"), df1, FUN = FUN)
  head(df2)
  tab = dcast(df2, formula("C ~ sample"))
  tab = tab[order(factor(tab$C, levels = unique(module.name$C))),]
  row.names(tab) = factor(tab$C, levels = unique(module.name$C)); tab$C = NULL
  return(tab[,rownames(df_in)])
}


load__M.e.TY <- function(existance) {
  M.e.TY = suppressWarnings(colnames(existance)[
    apply(existance[samples.log[rownames(existance),
                                "area"] == "southern slope",], 2, all)])
  M.e.TY = M.e.TY[!is.na(M.e.TY)]
  M.e.nTY = suppressWarnings(colnames(existance)[
    apply(existance[samples.log[rownames(existance),
                                "area"] == "central axis",], 2, all)])
  M.e.nTY = M.e.nTY[!is.na(M.e.nTY)]
  
  M.e.both = M.e.TY[M.e.TY %in% M.e.nTY]
  M.e.TY = M.e.TY[!M.e.TY %in% M.e.both]
  M.e.nTY = M.e.nTY[!M.e.nTY %in% M.e.both]
  
  return(list(both = M.e.both, TY = M.e.TY, nTY = M.e.nTY))
}
M.e.TY = load_if_not(M.e.TY, existance, .force.reload = TRUE)
M.e.both = M.e.TY[["both"]]
M.e.nTY = M.e.TY[["nTY"]]
M.e.TY = M.e.TY[["TY"]]
