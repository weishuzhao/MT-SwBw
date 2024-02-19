###
#' @Date: 2022-06-28 20:35:49
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2023-11-14 19:50:44
#' @FilePath: /2021_09-MT10kSW/workflow/others/draw_map.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


asc_path <- argv[1]
# asc_path = "data/GEBCO_28_Jun_2022_d37548b61e36/gebco_2021_n12.5_s10.0_w140.5_e143.5.asc"
fig_out1 <- argv[2]
# fig_out1 = file_path$file_path$figs("sample_site_map") %>% as.character()
fig_out2 <- argv[2]
# fig_out2 = file_path$file_path$figs("sample_site_map_1") %>% as.character()


co2 <- sp::read.asciigrid(asc_path)


## 绘图
par(mar = rep(0, 4)) %>%
  {
    sgdf2xyzc <- function(co3) {
      z <- co3$z

      # 用最小值和最大值, 把 x, y 包裹起来, 对应 z0
      x <- (seq_len(nrow(z)) * 10) %>%
        c(min(.) - 1e-10, ., max(.) + 1e-10)
      y <- (seq_len(ncol(z)) * 10) %>%
        c(min(.) - 1e-10, ., max(.) + 1e-10)

      # 用z0的颜色，把整个栅格包裹起来
      z0 <- min(z) - 20
      z1 <-
        z %>%
        cbind(z0, ., z0) %>%
        rbind(z0, ., z0) %>%
        `rownames<-`(NULL) %>%
        `colnames<-`(NULL)

      ## 创建用于显示颜色的矩阵
      # 默认全部使用绿色
      fcol <- matrix("green3", nrow = nrow(z1) - 1, ncol = ncol(z1) - 1)

      # 用灰色把所有的绿色都包裹起来,即设置四周的边界值
      i2 <- c(1, ncol(fcol))
      i1 <- c(1, nrow(fcol))
      fcol[, i2] <- "gray"
      fcol[i1, ] <- "gray"

      ## Take average of four neighboring values for palette
      ## 将上面设置的默认色，用都取相邻的四个格网颜色的平均值进行替换
      zi <- (
        co3$z[-1, -1] + co3$z[-1, -ncol(co3$z)] +
          co3$z[-nrow(co3$z), -1] + co3$z[-nrow(co3$z), -ncol(co3$z)]
      ) / 4
      pal <- terrain.colors(40, alpha = NULL)[
        cut(zi, quantile(zi, seq(0, 1, len = 41)), include.lowest = TRUE)
      ]

      fcol[-i1, -i2] <- pal
      return(list(x = x, y = y, z = z1, col = fcol))
    }
    co3 <- sp::as.image.SpatialGridDataFrame(co2)
    vars <- sgdf2xyzc(co3)
    # res = persp(
    res <- list(
      vars$x, vars$y, vars$z,
      col = vars$col,
      theta = 0, phi = 75, shade = 0.7, border = NA
    )
    par(.)
  }

sf <-
  as.data.frame(co2) %>%
  `colnames<-`(c("Depth", "Longitude", "Latitude"))


sample_site_meta <- sample_meta %>%
  list(
    metagenome = .,
    "16S" = read.csv("data/sample_meta_16S.tsv", sep = "\t") %>%
      `colnames<-`(
        c("Site", "Latitude", "Longitude", "Depth", "Group", "Type")
      ) %>%
      mutate(
        Location = ifelse(grepl("^S", get("Group")), "Slope", "Bottom"),
        Layers = ifelse(grepl("s$", get("Group")), "sediment", "water")
      )
  ) %>%
  bind_rows(.id = "Type") %>%
  group_by(Site, Group, Depth, Latitude, Longitude, Type) %>%
  summarise(Layers = paste(unique(Layers), collapse = ", ")) %>%
  ungroup() %>%
  mutate(
    LayersType = paste0(Layers, " (", Type, ")"), Layers = NULL, Type = NULL
  ) %>%
  group_by(
    Site, Group, Depth,
    Latitude = round(Latitude, 3), Longitude = round(Longitude, 3)
  ) %>%
  summarise(LayerTypes = paste(unique(LayersType), collapse = "; ")) %>%
  data.frame()

p <- ggplot(
  data = sf,
  mapping = aes_string(
    x = "Longitude", y = "Latitude",
    fill = "Depth",
    z = "Depth"
  )
) +
  geom_raster(interpolate = TRUE) +
  geom_contour(breaks = seq(-12000, -6000, 2000), color = "black") +
  scale_fill_gradientn(
    colors = paletteer::paletteer_c("grDevices::Plasma", 30),
    limits = c(-11000, -1000), n.breaks = 11
  ) +
  ggnewscale::new_scale_fill() +
  geom_point(
    data = sample_site_meta,
    mapping = aes_string(
      x = "Longitude", y = "Latitude", shape = "Type"
    ),
    color = "#BCAAA4",
    size = 2.5, alpha = 0.6
  ) +
  scale_color_gradientn(
    colors = paletteer::paletteer_c("ggthemes::Classic Orange-Blue", 30)
  ) +
  scale_fill_manual(values = sample_meta_col) +
  guides(color = guide_colorbar(title = "Sample Depth")) +
  labs(x = "Longitude (E)", y = "Latitude (N)")

font_size_1 <- 15
font_size_2 <- 13
font_size_3 <- 10
axis_ticks_length <- 0.1

p1 <-
  p +
  theme_bw() +
  theme(
    axis.text = element_text(
      size = font_size_2, colour = "black", face = "bold"
    ),
    axis.title = element_text(
      size = font_size_1, face = "bold", colour = "black"
    ),
    axis.ticks.length = unit(axis_ticks_length, "cm"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  theme(
    legend.title = element_text(size = font_size_2, face = "bold"),
    legend.text = element_text(size = font_size_3),
    legend.key = element_blank() #' element_rect(fill = "gray")
    #' legend.position = "bottom"
  ) +
  theme(text = element_text(
    family = "Arial",
    size = font_size_1,
    hjust = 0.5,
    lineheight = 0.5
  )) +
  theme(plot.title = element_text(hjust = 0.5))

p1_all <-
  p1 +
  geom_label_repel(
    data = sample_site_meta %>%
      filter(get("Latitude") < 11.30 | get("Longitude") > 142.5),
    mapping = aes_string(
      x = "Longitude", y = "Latitude",
      label = "Site", fill = "Group"
    ),
    min.segment.length = 0.3,
    max.overlaps = 40
  ) +
  geom_rect(
    data = . %>% filter(get("Depth") > 0),
    xmin = 142.17, xmax = 142.23, ymin = 11.31, ymax = 11.35,
    fill = NA, color = "white"
  ) +
  coord_fixed(xlim = c(141, 143), ylim = c(10.8, 11.8))
p1_bottom <-
  p1 +
  geom_contour(breaks = seq(-10800, -10200, 200), color = "black") +
  geom_label_repel(
    data = sample_site_meta %>%
      filter(get("Latitude") >= 11.30 & get("Longitude") <= 142.5),
    mapping = aes_string(
      x = "Longitude", y = "Latitude",
      label = "Site", fill = "Group"
    ),
    min.segment.length = 0.3,
    max.overlaps = 40
  ) +
  coord_fixed(xlim = c(142.17, 142.23), ylim = c(11.31, 11.35))


ggsave(fig_out1, p1_all, width = 12, height = 8)
ggsave(fig_out2, p1_bottom, width = 8, height = 5)
