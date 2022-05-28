###
#* @Date: 2021-08-22 21:53:56
#* @LastEditors: Hwrn
#* @LastEditTime: 2021-11-17 16:51:41
#* @FilePath: /Work/2021_09-MT10kSW/Analyze/Figs/02_func.r
#* @Description:
###
source("Analyze/Figs/00_func.r")


## function define
#' Get top div
#'
#' @param formula {formula} formula of tow columns in `div` to reshape
#'                          for example, sample ~ order
#' @param topi {int} top i `taxon` to keep in each sample
#' @param min.reads {int} min number of reads to keep
#' @param min.sample {int} min number of sample that a `taxon` in
#'
#' @return div.top {factor} name from `taxon` in `div`
#' @export
#'
#' @examples
taxon.div.top <- function(div,
                          taxon = "order",
                          value.var = "reads",
                          topi = 10,
                          min.value = 10,
                          min.sample = 2) {
  div.otu = reshape2::acast(div,
                            formula = formula(paste0("sample~", taxon)),
                            value.var = value.var,
                            fun.aggregate = sum)
  div.otu = div.otu[, apply(div.otu, 2, function(x) (sum(x) >= 0 &
                                                       sum(x > 0) >= 0))]


  if (topi < 1 || topi > length(colnames(div.otu)))
    topi = length(colnames(div.otu))

  div.top = unique(as.vector(
    apply(div.otu, 1, function(x)
      names(x[order(x, decreasing = T)][1:topi]))))
  return(div.top)
}


#' Demo of group all rows from `div` by `taxon` together
#'
#' @param div
#' @param taxon
#'
#' @return
#' @export
#'
#' @examples
taxon.group <- function(div, taxon = "order") {
  div.otu = reshape2::dcast(div,
                            formula(paste0("sample~", taxon)),
                            value.var = "reads",
                            fun.aggregate = sum)
  div.grouped = melt(
    div.otu,
    id.vars = "sample",
    variable.name = "order",
    value.name = "reads"
  )
  return(div.grouped)

}


#' Sort a vector of taxon according to `taxon.order`
#' WARNING: taxon.order must be provided with correct `dpco` and `name`
#'
#' @param names
#'
#' @return
#' @export
#'
#' @examples
taxon.factor.sort <- function(names) {
  name.short = unique(names)
  name.long = sapply(name.short,
                     function(x)
                       ifelse(
                         x %in% taxon.order$name,
                         as.character(taxon.order[x == taxon.order$name,
                                                  ]$dpco[1]),
                         "-------"
                       ))
  # set order here
  taxon.sort = factor(name.short, levels = name.short[order(name.long)])
  names = sapply(names, function(x)
    taxon.sort[which(x == taxon.sort)])

  return(names)
}


#' Change `reads` in `div` to "100%" style of each `sample`
#'
#' @param div {dataframe} include cols: sample, `taxon`, reads
#' @param taxon taxon to focus on, other will be ignored as "others"
#' @param taxon.sort sorted taxon (to `index`)
#' @param cutoff {float} between (0, 1)
#'
#' @return div {dataframe} include cols: sample, name, annot.percent,
#'                                       index, annot.percent,
#'                                       label.y, label.text
#'
#' @export
#'
#' @examples
taxon.pct.annot <- function(div,
                            taxon = "order",
                            value.var = "reads",
                            taxon.sort,
                            cutoff = 0.5,
                            div_is_pct = FALSE) {
  taxon.sort.g = sapply(sort(taxon.sort),
                        function(x) {
                          x = gsub("^(\\d)", "X\\1", gsub("-|;| ", ".", x))
                          return(x)})

  div[, taxon] = ifelse(div[, taxon] %in% taxon.sort,
                        as.character(div[, taxon]), "others")
  div.otu = reshape2::acast(div,
                            formula = formula(paste0("sample~", taxon)),
                            value.var = value.var,
                            fun.aggregate = sum,
                            check.names = FALSE)

  if (div_is_pct) { div.otu.pct = div.otu * 100
  } else {div.otu.pct = div.otu / apply(div.otu, 1, sum) * 100}

  div.grouped = reshape2::melt(
    data.frame(div.otu.pct,
               sample = rownames(div.otu.pct)),
    id.vars = "sample",
    variable.name = "name",
    value.name = "annot.percent"
  )
  div.grouped$name = sapply(
    div.grouped$name,
    function(x) sort(taxon.sort)[which(x == taxon.sort.g)[1]])
  div.grouped$name[is.na(div.grouped$name)] = "others"
  div.grouped = div.grouped[! is.na(div.grouped$annot.percent) &
                              div.grouped$annot.percent > 0, ]
  div.grouped$index = sapply(div.grouped$name,
                             function(x) which(x == sort(taxon.sort))[1])

  ce = arrange(div.grouped, sample, index, annot.percent)
  ce = ddply(ce, "sample", transform,
             label.y = cumsum(annot.percent) - annot.percent / 2)
  ce$label.text = ifelse(ce$annot.percent >= cutoff,
                         as.character(ce$name), NA)

  return(ce)
}


print(c("taxon.div.top",
        "taxon.group",
        "taxon.factor.sort",
        "taxon.pct.annot"))


#计算多种 Alpha 多样性指数，结果返回至向量
alpha_index <- function(x, method = 'richness', base = exp(1)) {
  if (method == 'richness') {  #丰富度指数 
    result = rowSums(x > 0)}
  else if (method == 'chao1') {  # Chao1 指数
    result = estimateR(x)[3, ]}
  else if (method == 'ace') {  # ACE 指数
    result = estimateR(x)[5, ]}
  else if (method == 'shannon') {  # Shannon 指数 
    result = diversity(x, index = 'shannon', base = base)}
  else if (method == 'simpson') {  # Gini-Simpson 指数
    result = diversity(x, index = 'simpson')}
  else if (method == 'pielou') {  # Pielou 均匀度
    result = diversity(x, index = 'shannon', 
                       base = base) / log(estimateR(x)[1, ], base)}
  else if (method == 'gc') {  # goods_coverage
    result = 1 - rowSums(x == 1) / rowSums(x)}
  return(result)
}


#根据抽样步长（step），统计每个稀释梯度下的 Alpha 多样性指数，结果返回至列表
alpha_curves <- function(x, step, method = 'richness',
                         rare = NULL, base = exp(1)) {
  x_nrow = nrow(x)
  rare = ifelse(is.null(rare), rowSums(x), rep(rare, x_nrow))
  alpha_rare = list()
  
  for (i in 1:x_nrow) {
    step_num = seq(0, rare[i], step)
    if (max(step_num) < rare[i]) {
      step_num = c(step_num, rare[i])}

    alpha_rare_i = NULL
    for (step_num_n in step_num) 
      alpha_rare_i = c(alpha_rare_i, alpha_index(
        x = rrarefy(x[i, ], step_num_n), method = method, base = base))
    names(alpha_rare_i) = step_num
    alpha_rare = c(alpha_rare, list(alpha_rare_i))
  }

  names(alpha_rare) = rownames(x)
  return(alpha_rare)
}

