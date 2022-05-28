###
#* @Date: 2021-08-22 21:53:45
#* @LastEditors: Hwrn
#* @LastEditTime: 2021-11-17 15:00:02
#* @FilePath: /2021_09-MT10kSW/Analyze/Figs/02_draw.r
#* @Description:
###
#rm(list = ls())
source('Analyze/Figs/02_func.r')
source('Analyze/Figs/02_load.r')
report = FALSE
#report = TRUE


load__div.all <- function(name.top) {
  div.16s = load__div.16s()
  div.kaiju = load__div.kaiju()

  div.all = rbind(
    data.frame(
      taxon.pct.annot(div.16s, taxon.sort = name.top, taxon = "name"),
      source = "16s"
    ),
    data.frame(
      taxon.pct.annot(div.kaiju, taxon.sort = name.top, taxon = "name"),
      source = "kaiju"
    )
  )
  div.all$name = sapply(div.all$name,
                        function(x) name.top[which(x == name.top)])

  return(div.all)
}


## 1. feature of diversity
### 1.1. compare difference and get consistent order
if (report) {
  compare__div.16s.kaiju(7)
}


### 1.1. things in each diversity estimating method

if (report) {
  order.16s = unique(div.16s$order)
  order.kaiju = unique(div.kaiju$order)
  order.rps3 = unique(div.rps3$order)
  order.common = order.16s[order.16s %in% order.kaiju]

  order.16s.top = taxon.div.top(div.16s, taxon = "name")
  order.kaiju.top = taxon.div.top(div.kaiju, taxon = "name")

  print("too few items in 'div.rps3', skip")
  print(unique(order.rps3))
  print("rps3 calculating method is also different")
  print("")

  print("only a few order share between '16s' and 'kaiju'")
  print(c("order.16s" = length(order.16s),
          "order.kaiju" = length(order.kaiju),
          "order.common" = length(order.common)))

  print(paste("if consider rps3, length", length(order.rps3)))
  print(order.rps3[order.rps3 %in% order.common])
  print("")

  sort(taxon.order[unlist(
    sapply(name.top, function(x) which(as.character(x) ==
                                         taxon.order$order))),
  ]$dpco)
  print(name.top == taxon.order$order)
}

### 1.summary.
print("Mainly compare div.16s and div.kaiju")

## 2. melt div to adjust
div.all = load__div.all(name.top)
print(names(div.all))

load__order.color <- function() {
  set.seed(0)
  order.names = levels(div.all$name)
  ndx = expand.grid(2:3/3, 0.5, 1:9/9)[1:15, ]
  mycolor = hsv(ndx[, 3], ndx[, 2], ndx[, 1])
  order.color = sapply(order.names,
                       function(x) mycolor[x == order.names])

  return(order.color)
}
order.color = load_if_not(order.color, .force.reload = TRUE)

p = ggplot(data = div.all) +
  geom_bar(aes_string(x = "sample", y = "annot.percent", fill = "name"),
           position = position_stack(reverse = T),
           stat = "identity", col = "black") +

  scale_fill_manual(values = order.color) +
  #scale_fill_manual(values = brewer_pal(palette = "BrBG")(17)) +

  scale_y_continuous(breaks = seq(0, 100, 10)) +
  labs(title = "", x = "sample", y = "precent %") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.~source)


p.1_2 = report__pcoa.16s()
p.1_3 = report__enrich.16s(order.color)
## 4. plot
#report = TRUE
if (report) {dev.off()}
if (report) {pdf("Analyze/Figs/02.pdf", width = 8, height = 6)}

### 2.1. match label and color to bar
p +
  guides(fill = "none") +
  geom_text(aes_string(x = "sample", y = "label.y",
                       label = "label.text"),
            position = 'identity', col = 'white', size = 2)

### 2.2. formal plot
p +
  guides(fill = guide_legend(reverse = TRUE))

## 3. PCA of 16s
p.1_3
if (report) {dev.off()}
if (report) {pdf("Analyze/Figs/02_1.pdf", width = 6, height = 6)}
p.1_2
if (report) {dev.off()}


p = ggplot(data = div.all[div.all$source == "16s",]) +
  geom_bar(aes_string(x = "sample", y = "annot.percent", fill = "name"),
           position = position_stack(reverse = T),
           stat = "identity", col = "black") +
  
  scale_fill_manual(values = order.color) +
  #scale_fill_manual(values = brewer_pal(palette = "BrBG")(17)) +
  
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  labs(title = "", x = "sample", y = "precent %") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.~source)

p +
  guides(fill = guide_legend(reverse = TRUE))
ggsave(filename = "Analyze/Figs/01B.png", device = "png", dpi = 300,
       width = 5, height = 5)
ggsave(filename = "Analyze/Figs/01C.png", plot = p.1_2,
       device = "png", dpi = 300,
       width = 5, height = 5)
ggsave(filename = "Analyze/Figs/01D.png", plot = p.1_3,
       device = "png", dpi = 300,
       width = 7.5, height = 5)


kraken.report = data.frame()
for (i in rownames(samples_info)) {
  kraken.report.i = read.csv(paste0("Pipe/", i, "-megahit/01_alpha_div/",
                                    "kraken-", i, "-megahit.report"),
                             sep = "\t", col.names = c(
                               "pct", "all.reads", "this.reads",
                               "tax.level", "uid", "name"))
  rownames(kraken.report) = kraken.report$uid
  kraken.report[kraken.report.i$uid, i] = kraken.report.i$this.reads
}
kraken.report[is.na(kraken.report)] = 0
kraken.report = kraken.report[apply(kraken.report, 1, sum) > 0,]
apply(kraken.report, 2, sum)

rarecurve(t(kraken.report), 2000000)

kraken.report.TY = data.frame()
for (i in c("TY.040", "TY.041", "TY.044")) {
  kraken.report.i = read.csv(paste0("Pipe/", i, "-megahit/01_alpha_div/",
                                    "kraken-", i, "-megahit.report"),
                             sep = "\t", col.names = c(
                               "pct", "all.reads", "this.reads",
                               "tax.level", "uid", "name"))
  kraken.report.TY[kraken.report.i$uid, i] = kraken.report.i$this.reads

  if ("TY" %in% names(kraken.report.TY)) {
    kraken.report.TY[kraken.report.i$uid, "TY"] = (
      kraken.report.i$this.reads + ifelse(
        is.na(kraken.report.TY[kraken.report.i$uid, "TY"]), 
        0, kraken.report.TY[kraken.report.i$uid, "TY"]))
  } else {
    kraken.report.TY[kraken.report.i$uid, "TY"] = (
      kraken.report.i$this.reads)
  }
}
kraken.report.TY[is.na(kraken.report.TY)] = 0

rarecurve(t(kraken.report.TY), 2000000)
