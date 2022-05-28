###
#* @Date: 2021-08-22 21:56:28
#* @LastEditors: Hwrn
#* @LastEditTime: 2021-09-20 19:49:38
#* @FilePath: /2021_09-MT10kSW/Analyze/Figs/00_func.r
#* @Description:
###
#rm(list = ls())

library(data.table)

library(reshape2)
library(plyr)
library(dplyr)
library(vegan)
library(ape)
library(stringr)

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(scales)
library(ggpubr)
suppressMessages(library(pathview))


is.defined <- function(var) {
  env = parent.frame()
  sym = deparse(substitute(var))
  return(exists(sym, env))
}


load_if_not <- function(var, ..., .force.reload=FALSE) {
  env = parent.frame()
  sym = deparse(substitute(var))
  if (!exists(sym, env)) {
    message(paste0("var '", sym, "' don't exists, loading"))
  } else if (.force.reload) {
    message(paste0("var '", sym, "' reloading"))
  } else {
    return(var)
  }
  env[[sym]] = env[[paste0("load__", sym)]](...)
  return(env[[sym]])
}


setcolor <- function(vec) {
  if (!class(vec) == "factor") {
    vec = as.factor(vec)
  }
  vec.color = data.frame(level = as.character(levels(vec)),
                         color = rainbow(length(levels(vec))),
                         stringsAsFactors = F)
  vec.color = sapply(vec, function(x) {
    vec.color$color[which(x == vec.color$level)]})

  return(vec.color)
}


taxon.as.num <- function(taxon) {
  if (class(taxon) == "character") {
    taxon = which(tolower(strsplit(taxon, "")[[1]][1]) == c("d", "p", "c", "o",
                                                            "f", "g", "s"))
  }
  return(taxon)
}


taxon.split <- function(taxon.full, start, end = NULL) {
  if (is.null(end))
    end = start

  start = taxon.as.num(start)
  end = taxon.as.num(end)

  taxon.new = sapply(as.character(taxon.full), function(x) {
    if (is.na(x)) return(x)
    if (strsplit(x, "^.__") == x) {
      tax_order = unlist(strsplit(as.character(x), ";"))
    } else {
      tax_order = unlist(strsplit(unlist(strsplit(x, "^.__"))[2], ";.__"))
    }
    return(paste(c(tax_order, rep("", end))[start:end], collapse = ";"))
  })

  taxon.new.factor = factor(taxon.new,
                            levels = unique(taxon.new[order(taxon.full)]))
  return(taxon.new.factor)
}
