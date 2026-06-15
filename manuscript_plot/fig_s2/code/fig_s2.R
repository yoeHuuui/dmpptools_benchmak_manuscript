library(ggplot2)
library(tidyverse)
library(scales)

data <- read.table("fig_s2/data/benchmark_promoter_scale.time.tsv",
                   header = T)
query_time2 <- data %>%
  mutate(
    elapsed_sec = as.numeric(elapsed_sec),
    sample_size = as.numeric(sample_size),
    sample_number = factor(
      sample_size,
      levels = c(1, 5, 10, 50, 100, 500, 1000, 2000, 5000)),
    tool = factor(
      tool,
      levels = c("dmpp", "wgbstools","ballc"),
      labels = c("dmpp",  "wgbstools", "BAllCools"))) %>%
  filter(exit_status == 0, !is.na(elapsed_sec), elapsed_sec > 0)

p <- ggplot()+
  stat_summary(data = query_time2, 
               aes(x = sample_number, y = elapsed_sec, 
                   color = tool, group = tool),
               fun = mean, geom = "line", linewidth = 0.6)+
  stat_summary(data = query_time2,
               aes(x = sample_number, y = elapsed_sec, color = tool, group = tool),
               fun = mean, geom = "point", size = 1.8)+
  stat_summary(data = query_time2,
               aes(x = sample_number, y = elapsed_sec, color = tool, group = tool),
               fun.data = mean_sdl,
               fun.args = list(mult = 1),
               geom = "errorbar", width = 0, linewidth = 0.4, 
               show.legend = FALSE)+
  scale_y_log10(
    breaks = c(0.1,1,10,100, 100),
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.05, 0.12)))+
  labs(x = "Number of queried regions",
       y = "Wall-clock time (s)", fill = NULL, color = NULL) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 13, face = "bold", color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.45, color = "black"),
    axis.ticks = element_line(linewidth = 0.45, color = "black"),
    legend.position = "top",
    panel.spacing.x = unit(1.0, "lines"),
    plot.margin = margin(8, 12, 8, 8)
  )

pdf("fig_s2/fig_s2.pdf", width = 6, height = 5)
p
dev.off()
