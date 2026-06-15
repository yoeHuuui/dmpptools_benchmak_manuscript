library(ggplot2)
library(tidyverse)
library(scales)

file_list <-list.files("./fig_2/data/")
file_list <- file_list[grep("time.tsv", file_list)]

data_list <- lapply(file_list, function(x){
  path <- paste0("./fig_2/data/", x)
  tmp <- read.table(path, header = T)
  
})

data <- do.call(rbind, data_list)

plot_df <- data %>%
  mutate(
    elapsed_sec = as.numeric(elapsed_sec),
    exit_status = as.character(exit_status),
    status = as.character(status),
    benchmark = case_when(
      benchmark %in% c("benchmark_sine") ~ "SINE",
      benchmark %in% c("benchmark_line") ~ "LINE",
      benchmark %in% c("benchmark_ccre") ~ "cCRE",
      benchmark %in% c("benchmark_promoter") ~ "Promoter",
      TRUE ~ benchmark),
    benchmark = factor(benchmark,
                       levels = c("Promoter", "SINE", "LINE", "cCRE")),
    tool = case_when(
      tool %in% c("dmpp", "DMPP") ~ "dmpp",
      tool %in% c("wgbstools", "wgbs", "WGBS") ~ "wgbstools",
      TRUE ~ tool),
    tool = factor(tool, levels = c("dmpp", "wgbstools"))) %>%
  filter(status == "OK", exit_status == "0", !is.na(elapsed_sec),
         elapsed_sec > 0, !is.na(benchmark), !is.na(tool))

paired_df <- plot_df %>%
  group_by(benchmark, sample, tool) %>%
  summarise(elapsed_sec = median(elapsed_sec, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = tool, values_from = elapsed_sec) %>%
  filter(!is.na(dmpp), !is.na(wgbstools), dmpp > 0, wgbstools > 0) %>%
  mutate(speedup = wgbstools / dmpp)

speedup_summary <- paired_df %>%
  group_by(benchmark) %>%
  summarise(
    n_samples = n(),
    median_speedup = median(speedup, na.rm = TRUE),
    p25_speedup = quantile(speedup, 0.25, na.rm = TRUE),
    p75_speedup = quantile(speedup, 0.75, na.rm = TRUE),
    median_dmpp_sec = median(dmpp, na.rm = TRUE),
    median_wgbstools_sec = median(wgbstools, na.rm = TRUE),
    p95_dmpp_sec = quantile(dmpp, 0.95, na.rm = TRUE),
    p95_wgbstools_sec = quantile(wgbstools, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(round(median_speedup, 1), "\u00d7")
  )

print(speedup_summary)

label_df <- plot_df %>%
  group_by(benchmark) %>%
  summarise(
    y = max(elapsed_sec, na.rm = TRUE) * 1.8, .groups = "drop"
  ) %>%
  left_join(speedup_summary, by = "benchmark") %>%
  mutate(
    label = paste0(label)
  )

p <- ggplot(plot_df, aes(x = benchmark, y = elapsed_sec, fill = tool)) +
  geom_boxplot(width = 0.62, outlier.shape = NA,  alpha = 0.85,
               position = position_dodge(width = 0.72)) +
  geom_point(aes(color = tool),
             position = position_jitterdodge(jitter.width = 0.12,
                                             jitter.height = 0,
                                             dodge.width = 0.72),
             size = 1.1, alpha = 0.35, stroke = 0)+
  geom_text(data = label_df, aes(x = benchmark, y = y, label = label),
            inherit.aes = FALSE, size = 4.0, fontface = "bold",
            lineheight = 0.9)+
  scale_y_log10(breaks = c(0.3, 1, 3, 10, 30, 100, 300, 1000),
                labels = label_number(accuracy = 0.1),
                expand = expansion(mult = c(0.04, 0.22)))+
  scale_fill_manual(values = c("dmpp" = "#4C78A8", "wgbstools" = "#F58518"))+
  scale_color_manual(values = c("dmpp" = "#4C78A8", "wgbstools" = "#F58518"))+
  labs(x = NULL, y = "Wall-clock time (s)",
       fill = NULL, color = NULL) +
  theme_classic()+
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    legend.position = "top",
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    plot.margin = margin(8, 12, 6, 8))

pdf("fig_2/fig_2b.pdf", width = 5, height = 4)
p
dev.off()


## 
data <- read.table("fig_2/data/benchmark_4_overlap_sweep_overlap.tsv",
                   header = T)
data <- data[data$tool %in% c("dmpp", "wgbstools"), ]

data %>% group_by(benchmark, sample, tool) %>% 
  summarise(elapsed_sec=mean(elapsed_sec)) -> data

data$benchmark <- str_split(data$benchmark, "_",simplify = T)[,6]
data$benchmark <- str_replace(data$benchmark, "pct", "%")

p <- ggplot()+
  stat_summary(data = data, 
               aes(x = benchmark, y = elapsed_sec, 
                   color = tool, group = tool),
               fun = mean, geom = "line", linewidth = 0.6)+
  stat_summary(data = data,
               aes(x = benchmark, y = elapsed_sec, color = tool, group = tool),
               fun = mean, geom = "point", size = 1.8)+
  stat_summary(data = data,
               aes(x = benchmark, y = elapsed_sec, color = tool, group = tool),
               fun.data = mean_sdl,
               fun.args = list(mult = 1),
               geom = "errorbar", width = 0, linewidth = 0.4, 
               show.legend = FALSE)+
  scale_color_manual(values = c("dmpp" = "#4C78A8", "wgbstools" = "#F58518"))+
  labs(x = NULL, y = "Wall-clock time (s)",
       fill = NULL, color = NULL) +
  theme_classic()+
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    legend.position = "top",
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    plot.margin = margin(8, 12, 6, 8))

pdf("fig_2/fig_2c.pdf", width = 5, height = 4)
p
dev.off()
