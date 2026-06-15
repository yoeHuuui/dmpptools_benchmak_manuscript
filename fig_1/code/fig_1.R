library(tidyverse)
library(scales)

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

total.df <- filesize.df %>%
  group_by(group) %>%
  summarise(total_gb = sum(size_gb), .groups = "drop") %>%
  mutate(
    group = factor(group, levels = format_levels)
  )

allc_size <- total.df %>%
  filter(group == "ALLC") %>%
  pull(total_gb)

total.df <- total.df %>%
  mutate(
    fold = allc_size / total_gb,
    label_size = ifelse(
      total_gb >= 1024,
      paste0(round(total_gb / 1024, 1), " TB"),
      paste0(round(total_gb, 1), " GB")
    )
  )

compare_groups <- c("dmpp", "wgbstools", "dmpp (single)")

anno.df <- tibble(
  group_ref = c("ALLC", "wgbstools"),
  group_cmp = c("dmpp", "dmpp (single)")
) %>%
  left_join(
    total.df %>% select(group, ref_size = total_gb),
    by = c("group_ref" = "group")
  ) %>%
  left_join(
    total.df %>% select(group, cmp_size = total_gb),
    by = c("group_cmp" = "group")
  ) %>%
  mutate(
    xstart = as.numeric(factor(group_ref, levels = format_levels)),
    xend = as.numeric(factor(group_cmp, levels = format_levels)),
    fold = ref_size / cmp_size,
    y = pmax(ref_size, cmp_size) * c(1.10, 1.25),
    fold_label = paste0(round(fold, 1), "-fold reduction")
  )

p <- ggplot(total.df, aes(x = group, y = total_gb, fill = group))+
  geom_col(width = 0.62, color = "black", linewidth = 0.35)+
  geom_text(aes(label = label_size), vjust = -0.5, size = 3.3)+
  geom_segment(data = anno.df,
               aes(x = xstart, xend = xend, y = y, yend = y),
               inherit.aes = FALSE,
               linewidth = 0.4)+
  geom_text(data = anno.df,
            aes(x = (xstart + xend) / 2, y = y * 1.08, label = fold_label),
            inherit.aes = FALSE, size = 3.1, fontface = "bold") +
  scale_fill_manual(values = format_cols)+
  scale_y_continuous(
    labels = label_number(scale_cut = cut_short_scale()),
    expand = expansion(mult = c(0, 0.22)))+
  labs(x = NULL, y = "Total file size (GB)")+
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.y = element_text(color = "black", size = 11),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.45)
  )

pdf("fig_1/fig.1d.pdf", width = 4.5, height = 3.5)
p
dev.off()


allc <- read.table("fig_1/data/allc_filesize_20260512.txt")
dmpp <- read.table("fig_1/data/dmpp_filesize_20260512.txt")

allc$sample <- str_split(allc$V1, "\\/", simplify = T)[,3]
allc <- allc[grep("allc", allc$sample), ]
allc$group <- "ALLC"
dmpp$sample <- str_split(dmpp$V1, "\\/", simplify = T)[,3]
dmpp <- dmpp[grep("dmpp", dmpp$sample), ]
dmpp$group <- "dmpp"
filesize.df <- rbind(allc, dmpp)

total.df <- filesize.df %>%
  group_by(group) %>%
  summarise(total = sum(V2), .groups = "drop")

total.df$total_Tb <-  total.df$total/1024/1024/1024/1024
total.df$label_size <- paste0(round(total.df$total_Tb,1), " TB")

format_cols <- c("ALLC" = "#dbe9f6",
                 "dmpp" = "#3f88c5")

anno.df <- tibble(
  group_ref = c("ALLC"),
  group_cmp = c("dmpp")
) %>%
  left_join(
    total.df %>% select(group, ref_size = total_Tb),
    by = c("group_ref" = "group")
  ) %>%
  left_join(
    total.df %>% select(group, cmp_size = total_Tb),
    by = c("group_cmp" = "group")
  ) %>%
  mutate(
    xstart = as.numeric(factor(group_ref, levels = c("ALLC", "dmpp"))),
    xend = as.numeric(factor(group_cmp, levels = c("ALLC", "dmpp"))),
    fold = ref_size / cmp_size,
    y = pmax(ref_size, cmp_size) * c(1.10),
    fold_label = paste0(round(fold, 1), "-fold reduction")
  )

p <- ggplot(total.df, aes(x = group, y = total_Tb, fill = group))+
  geom_col(width = 0.62, color = "black", linewidth = 0.35)+
  geom_text(aes(label = label_size), vjust = -0.5, size = 3.3)+
  geom_segment(data = anno.df,
               aes(x = xstart, xend = xend, y = y, yend = y),
               inherit.aes = FALSE,
               linewidth = 0.4)+
  geom_text(data = anno.df,
            aes(x = (xstart + xend) / 2, y = y * 1.08, label = fold_label),
            inherit.aes = FALSE, size = 3.1, fontface = "bold") +
  scale_fill_manual(values = format_cols)+
  scale_y_continuous(
    labels = label_number(scale_cut = cut_short_scale()),
    expand = expansion(mult = c(0, 0.22)))+
  labs(x = NULL, y = "Total file size (TB)")+
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.y = element_text(color = "black", size = 11),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.45)
  )

pdf("fig_1/fig.1e.pdf", width = 3.5, height = 3.5)
p
dev.off()