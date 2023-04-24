###
#' @Date: 2022-06-28 20:35:49
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-08-05 12:55:00
#' @FilePath: /2021_09-MT10kSW/workflow/others/draw_map.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


asc_path = argv[1]
#asc_path = "data/GEBCO_28_Jun_2022_d37548b61e36/gebco_2021_n12.5_s10.0_w140.5_e143.5.asc"
fig_out = argv[2]
#fig_out = file_path$file_path$figs("sample_site_map") %>% as.character()


co2 = sp::read.asciigrid(asc_path)
co3 = sp::as.image.SpatialGridDataFrame(co2)

SGdf_to_xyzc <- function(co3) {
  z = co3$z

  #用最小值和最大值, 把 x, y 包裹起来, 对应 z0
  x = 1:nrow(z) %>% {10 * .} %>% {c(min(.) - 1e-10, ., max(.) + 1e-10)}
  y = 1:ncol(z) %>% {10 * .} %>% {c(min(.) - 1e-10, ., max(.) + 1e-10)}

  # 用z0的颜色，把整个栅格包裹起来
  z0 = min(z) - 20
  z = rbind(z0, cbind(z0, z, z0), z0)

  ## 创建用于显示颜色的矩阵
  #默认全部使用绿色
  fcol <- matrix("green3", nrow = nrow(z) - 1, ncol = ncol(z) - 1)

  #用灰色把所有的绿色都包裹起来,即设置四周的边界值
  fcol[, i2 <- c(1,ncol(fcol))] <- "gray"
  fcol[i1 <- c(1,nrow(fcol)), ] <- "gray"

  ## Take average of four neighboring values for palette
  ## 将上面设置的默认色，用都取相邻的四个格网颜色的平均值进行替换
  zi <- (co3$z[-1,-1] + co3$z[-1,-ncol(co3$z)] + co3$z[-nrow(co3$z),-1] + co3$z[-nrow(co3$z),-ncol(co3$z)])/4
  pal <- terrain.colors(40, alpha = NULL)[cut(zi, quantile(zi, seq(0,1, len = 41)), include.lowest = TRUE)]

  fcol[-i1,-i2] <- pal
  return(list(x = x, y = y, z = z, col = fcol))
}

## 绘图
#par(mar=rep(0,4))
#vars = SGdf_to_xyzc(co3)
#res = persp(vars$x, vars$y, vars$z, col = vars$col,
#            theta = 0, phi = 75, shade = 0.7, border = NA)


sf =
  as.data.frame(co2) %>%
  {colnames(.) = c("Depth", "Longitude", "Latitude"); .} %>%
  {message.print(summary(.)); .}
sample_site_meta =
  sample_meta %>%
  {.[!duplicated(.$Site),]}

p =
  ggplot(data = sf,
         mapping = aes_string(x = "Longitude", y = "Latitude",
                              fill = "Depth",
                              z = "Depth")) +
  geom_raster(interpolate = TRUE) +
  geom_contour(breaks = seq(-11000, -5000, 1000), color = "black") +
  scale_x_continuous(limits = c(141, 143)) +
  scale_y_continuous(limits = c(10.5, 12)) +
  scale_fill_gradientn(colors = paletteer::paletteer_c("grDevices::Plasma", 30),
                       limits = c(-11000, -1000), n.breaks = 11) +
  geom_point(data = sample_site_meta,
             mapping = aes_string(x = "Longitude", y = "Latitude",
                                  color = "-Depth", shape = "Group"),
             size = 2.5, alpha = 0.6) +
  scale_color_gradientn(colors = paletteer::paletteer_c("ggthemes::Classic Orange-Blue", 30)) +
  scale_shape_manual(values = c("Bs" = 15, "Bw" = 0, "Ss" = 17, "Sw" = 2)) +
  ggnewscale::new_scale_fill() +
  geom_label_repel(data = sample_site_meta,
                   mapping = aes_string(x = "Longitude", y = "Latitude",
                                        label = "Site", fill = "Group"),
                   min.segment.length = 0.3,
                   max.overlaps = 40) +
  scale_fill_manual(values = sample_meta_col) +
  guides(color = guide_colorbar(title = "Sample Depth")) +
  labs(x = "Longitude (E)", y = "Latitude (N)")


font_size_1 = 15
font_size_2 = 13
font_size_3 = 10
axis.ticks.length = 0.1

p1 =
  p +
  theme_bw() +
  theme(
    axis.text = element_text(size = font_size_2, colour = "black", face = "bold"),
    axis.title = element_text(size = font_size_1, face = "bold", colour = "black"),
    axis.ticks.length = unit(axis.ticks.length, 'cm'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  theme(
    legend.title = element_text(size = font_size_2, face = "bold"),
    legend.text = element_text(size = font_size_3),
    legend.key = element_blank()  # element_rect(fill = "gray")
    #legend.position = "bottom"
  ) +
  theme(text = element_text(family = "Arial",
                            size = font_size_1,
                            hjust = 0.5,
                            lineheight = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5))


ggsave(fig_out, p1, width = 12, height = 8)
