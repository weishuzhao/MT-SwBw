###
#* @Date: 2021-08-22 21:53:45
#* @LastEditors: Hwrn
#* @LastEditTime: 2021-09-22 09:38:42
#* @FilePath: /2021_09-MT10kSW/Analyze/Figs/03_draw.r
#* @Description:
###
source('Analyze/Figs/03_func.r')
source('Analyze/Figs/03_load.r')
report = FALSE


report__90_ZFMG_S <- function() {
  Wtdb = load__Wtdb(usePrimary = TRUE)
  dropSedimant = which(!annotation_row$site == "Trench Water")
  s = summarize__Wtdb.sample(Wtdb, data.frame(
    row.names = rownames(annotation_row)[dropSedimant],
    site = annotation_row[dropSedimant, ]
  ))
  return(s)
}#; report__90_ZFMG_S()


report__95_ZFMG_W <- function() {
  Wtdb = load__Wtdb()
  dropSedimant = which(!annotation_row$site == "Trench Sediment")
  s = summarize__Wtdb.sample(Wtdb, data.frame(
    row.names = rownames(annotation_row)[dropSedimant],
    site = annotation_row[dropSedimant, ]
  ))
  return(s)
}#; report__95_ZFMG_W()


report__90_ZFMG_W <- function() {
  Wtdb = load__Wtdb(usePrimary = TRUE)
  dropSedimant = which(!annotation_row$site == "Trench Sediment")
  s = summarize__Wtdb.sample(Wtdb, data.frame(
    row.names = rownames(annotation_row)[dropSedimant],
    site = annotation_row[dropSedimant, ]
  ))
  return(s)
}#; report__90_ZFMG_W()


colnames(Wtdb)
annotation_row = data.frame(
  row.names = c(colnames(Wtdb)[16:27]),
  site = c(rep("TY", 3), rep("WQ", 5), rep("YW", 4)),
  stringsAsFactors = FALSE)


report = FALSE
#report = TRUE
if (report) {dev.off()}
if (report) {pdf("Analyze/Figs/03.pdf", width = 11, height = 6)}
s95 = summarize__Wtdb(Wtdb)
s95[[1]] +
  annotate("segment", x =  0, xend = 10, y = 75, yend =  75) +
  annotate("segment", x =  0, xend = 10, y = 90, yend =  90) +
  annotate("segment", x =  5, xend =  5, y = 90, yend = 100) +
  annotate("segment", x =  0, xend =  5, y = 95, yend =  95)

s95[[2]]
s95[[3]]


s95.s = summarize__Wtdb.sample(Wtdb, annotation_row)
s95.s[[2]]
s95.s[[3]]


report__div.Wtdb(Wtdb, value = "score", FUN = max)
report__div.Wtdb(Wtdb, value = "score", FUN = length) +
  guides(size = guide_legend(title = "Bin number"))


Wtdb.cluster = Wtdb$cluster[Wtdb$comp_cont %in% c(">90, <5", ">95, <5") &
                              Wtdb$annot_level == "species"]
report__div.Stdb(Wtdb.cluster) +
  theme(axis.text.y = element_text(hjust = 0),
        text = element_text(size = 6))

if (report) {dev.off()}

ws.clusters = sapply(c(
  "Zunongwangia profunda", "Nocardioides sp001627335", "Brevibacterium casei",
  "NORP165 sp002400895", "PRS1 sp003709165", "CR02bin9 sp004356555",
  "Nitrosopumilus sp013390905", "UBA2126 sp004124385"
  ), function(x) Wtdb[Wtdb$name == x, "cluster"])


Stdb.rep = data.frame(
  Stdb[c("genome", "score", "comp_cont", "sample")],
  raw_bases = Stdb$AvgDepth * Stdb$length
) %>% 
  group_by(sample) %>% 
  arrange(sample, -score) %>% 
  mutate(cum_num = cumsum(raw_bases)) %>% 
  ungroup()
Stdb.rep$pct_bases = Stdb.rep$raw_bases / 
  data.frame(samples_info, 
             adj_bases = samples_info$depth * samples_info$GenomeSize
             )[Stdb.rep$sample, "adj_bases"]

Stdb.rep$pct_bases.adj = Stdb.rep$pct_bases *
  data.frame(samples_info, 
             adj_bases = samples_info$Mapped.reads / samples_info$Reads
  )[Stdb.rep$sample, "adj_bases"]


ggplot(Stdb.rep[as.integer(Stdb.rep$comp_cont) > 4, ]) +
  geom_bar(
    mapping = aes_string(
      x = "sample", y = "pct_bases.adj",
      group = "order(rev(order(score)))",
      fill = "comp_cont"),
    position = position_stack(reverse = T), stat = "identity") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("Analyze/Figs/04S2.png",
       device = "png", dpi = 300,
       width = 8, height = 5)


report__div.Wtdb(Wtdb, value = "score", FUN = max)
ggsave(filename = "Analyze/Figs/02A.png", device = "png", dpi = 300,
       width = 7, height = 4)

report__div.Wtdb(Wtdb, value = "score", FUN = length) +
  guides(size = guide_legend(title = "Bin number"))
ggsave(filename = "Analyze/Figs/02B.png", device = "png", dpi = 300,
       width = 7, height = 4)


Wtdb.cluster = Wtdb$cluster[Wtdb$comp_cont %in% c(">90, <5", ">95, <5") &
                              Wtdb$annot_level == "species"]
report__div.Stdb(Wtdb.cluster) +
  theme(axis.text.y = element_text(hjust = 0),
        text = element_text(size = 6))
ggsave(filename = "Analyze/Figs/02C.png", device = "png", dpi = 300,
       width = 7, height = 4)
