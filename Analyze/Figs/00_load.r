###
#* @Date: 2021-11-06 21:39:41
#* @LastEditors: Hwrn
#* @LastEditTime: 2021-11-08 21:00:26
#* @FilePath: /2021_09-MT10kSW/Analyze/Figs/00_load.r
#* @Description:
###
source("Analyze/Figs/00_func.r")


load__samples.log <- function() {
  samples.log_1 = read.csv("~/Work/2021_09-ZFMG_MG/00_data/sample_meta.tsv",
                           na.strings = "---",
                           sep = "\t", header = 1, as.is = TRUE)
  samples.log_2 = read.csv("00_data/sample_meta.tsv",
                           na.strings = "——",
                           sep = "\t", header = 1, as.is = TRUE)
  samples.log = merge(x = samples.log_1, by.x = "X.Sample", all = FALSE,
                      y = samples.log_2, by.y = "Station"
                      )[c(1:13, 19:23)]

  colnames(samples.log)
  rownames(samples.log) = samples.log$X.Sample

  message("We only use the water samples")
  samples.log = samples.log[samples.log$Sample.Type == "water", ]
  samples.log$area = sapply(samples.log$X.Sample, function(x)
    ifelse(unlist(strsplit(x, "\\."))[1] == "TY",
           "southern slope", "central axis"))
  
  return(samples.log)
}
samples.log = load_if_not(samples.log)


load__samples_info <- function() {
  samples_info = read.csv("~/Work/2021_09-ZFMG_MG/Oerr/sample_info.tsv",
                         sep = "\t", header = 1, as.is = TRUE)
  rownames(samples_info) = samples_info$key
  samples_info = samples_info[rownames(samples.log), ]

  return(samples_info)
}
samples_info = load_if_not(samples_info)
