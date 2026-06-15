library(ggplot2)
library(tidyverse)
library(scales)

query_time <- read.table("fig_4/data/query_time.all.tsv",
                         header = T) %>%
  filter(exit_status == 0) %>%
  mutate(elapsed_sec = as.numeric(elapsed_sec),
         sample_number = as.character(sample_size),
         tool = factor(tool, levels = c("dmpp_list", "ballc")))

query_time2 <- query_time %>%
  mutate(
    elapsed_sec = as.numeric(elapsed_sec),
    sample_size = as.numeric(sample_size),
    sample_number = factor(
      sample_size,
      levels = c(10, 20, 50, 100, 200, 500, 1000)
    ),
    context = factor(context, levels = c("CG", "CHG", "CHH")),
    tool = factor(
      tool,
      levels = c("dmpp_list", "ballc"),
      labels = c("dmpp", "BAllCools"))) %>%
  filter(exit_status == 0, !is.na(elapsed_sec), elapsed_sec > 0)


plot_df <- query_time2 %>%
  mutate(sample_size = as.numeric(sample_size),
         tool = factor(tool, levels = c("dmpp", "BAllCools"))) %>%
  group_by(context, sample_size, tool) %>%
  summarise(n = n(),
            p5 = quantile(elapsed_sec, 0.05, na.rm = TRUE),
            median = median(elapsed_sec, na.rm = TRUE),
            p95 = quantile(elapsed_sec, 0.95, na.rm = TRUE),
            .groups = "drop")

p <- ggplot(plot_df, aes(x = sample_size, y = median, 
                         color = tool, group = tool)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = p5, ymax = p95), linewidth = 0.8,
                 alpha = 0.75) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  scale_x_log10(breaks = sort(unique(plot_df$sample_size)))+
  facet_wrap(. ~context, nrow=3) +
  scale_color_manual(values = c("dmpp" = "#4C78A8","BAllCools" = "#F58518"))+
  labs(x = "Number of single cells",
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

pdf("fig_4/fig_4a.pdf", height = 9,  width = 4)
p
dev.off()

p95.df <- query_time2 %>%
  group_by(tool, context, sample_size) %>%
  summarise(P95 = quantile(elapsed_sec, 0.95, names = FALSE), .groups = "drop")

speedup.tbl <- p95.df %>%
  mutate(dmpp_P95 = P95[tool == "dmpp"], .by = c(context, sample_size)) %>%
  filter(tool != "dmpp") %>%
  mutate(
    P95 = round(P95, 2),
    dmpp_P95 = round(dmpp_P95, 2),
    speedup_vs_dmpp = round(P95 / dmpp_P95, 1)
  ) %>%
  arrange(context, sample_size, tool)

write.table(speedup.tbl, "fig_4/query_context_summary.tsv",
            quote = F, sep = "\t", row.names = F, col.names = T)


query_df <- read.table("fig_4/data/gene_query_n1000_c1_20260607_134038.tsv",
                       sep = "\t", header = T)
query_df <- subset(query_df, phase != "warmup")
p95 <- quantile(query_df$seconds, 0.95, names = FALSE)
p <- ggplot(query_df, aes(x=phase, y=seconds))+
  geom_violin(width=0.5, color="#4C78A8")+
  geom_boxplot(width=0.2, color="#4C78A8")+
  geom_hline(yintercept = p95, linetype="dashed")+
  geom_text(aes(x=phase, y=p95+3),label=paste0("p95 = ", round(p95,2)))+
  ylim(2,12)+
  xlab("")+ylab("Query time (s)")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pdf("fig_4/fig_4b.pdf", width = 2.5, height = 3)
p
dev.off()

