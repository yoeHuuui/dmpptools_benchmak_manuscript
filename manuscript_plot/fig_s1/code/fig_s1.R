library(tidyverse)
library(scales)
library(ggplot2)

summary <- read.table("fig_s1/data/hg38.transition.summary.tsv",
                      sep = "\t", header = TRUE)

summary <- subset(summary, byte_class == "1 byte")

summary$delta <- 1 - summary$fraction

eps <- 1e-12
summary$delta_plot <- pmax(summary$delta, eps)

pdf("fig_s1/fig.s1a.pdf", width = 4, height = 3)
ggplot(summary,
       aes(x = context_prev, y = context_next, fill = delta_plot))+
  geom_tile(color = "white",  linewidth = 0.7)+
  geom_text(aes(label = ifelse(delta == 0, "0", 
                               scales::label_scientific(digits = 2)(delta))),
            size = 3.5)+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  scale_fill_gradient(low = "#dbe9f6",high = "#2b6cb0",
                      trans = "log10",
                      labels = scientific_format(digits = 2),
                      name = "> 1 bytes"
  )+
  theme_bw(base_size = 11) +
  labs(x = "Previous context", y = "Next context")+
  theme(axis.text = element_text(color = "black"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

dev.off()

anno <- c("allc" = "ALLC", "ballc" = "BALLC", 
          "dmpp" = "dmpp", "wgbstools" = "wgbstools",
          "dmpp (single)" = "dmpp (single)")

format_levels <- c("ALLC", "BALLC", "dmpp",
                   "wgbstools","dmpp (single)")

format_cols <- c("ALLC" = "#dbe9f6", "BALLC" = "#a6c8e0",
                 "dmpp" = "#3f88c5", "wgbstools" = "#6f73c8",
                 "dmpp (single)" = "#b59ad9")

filesize <- read.table("fig_1/data/filesize_20260512.txt")

filesize$group <- str_split(filesize$V1, "\\/", simplify = TRUE)[,3]
filesize$sample <- str_split(filesize$V1, "\\/", simplify = TRUE)[,2]

filesize[grep("nostrand", filesize$group), "group"] <- "dmpp (single)"

filesize.df <- filesize %>%
  group_by(sample, group) %>%
  summarise(size = sum(V2), .groups = "drop") %>%
  filter(!group %in% c("yame", "CGbz", "cgmap")) %>%
  mutate(
    group = unname(anno[group]),
    group = factor(group, levels = format_levels),
    size_gb = size / 1024^3
  ) %>%
  filter(!is.na(group))


pdf("fig_s1/fig.s1b.pdf", width = 4.5, height = 3.5)
ggplot(filesize.df, aes(x=group, y=size/1024/1024, fill = group))+
  geom_boxplot()+
  scale_fill_manual(values = format_cols)+
  xlab("")+ylab("File size (MB)")+
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.y = element_text(color = "black", size = 11),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.45)
  )
dev.off()
