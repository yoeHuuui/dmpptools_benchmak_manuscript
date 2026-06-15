library(tidyverse)
library(readxl)
library(scales)

data <- read_excel("fig_s7/data/on_fly_train.xlsx")
colnames(data)[1] <- "Number"

data <- data %>%
  mutate(
    Number = factor(Number, levels = c("1K", "5K", "10K", "50K", "93K")),
    memmap_prepare = make_dmpp_bed + bed_to_parquet + parquet_to_memmap,
    memmap_total = memmap_prepare + train_memmap,
    dmpp_total = train_dmpp,
    speedup = memmap_total / dmpp_total
  )

plot.df <- bind_rows(
  data %>%
    transmute(Number, workflow = "Memmap workflow", 
              component = "Preparation", time_sec = memmap_prepare),
  data %>%
    transmute(Number, workflow = "Memmap workflow",
              component = "Training", time_sec = train_memmap),
  data %>%
    transmute(Number, workflow = "dmpp on-the-fly",
              component = "On-the-fly training", time_sec = train_dmpp)) %>%
  mutate(
    time_hour = time_sec / 3600,
    workflow = factor(workflow, levels = c("Memmap workflow", 
                                           "dmpp on-the-fly")),
    component = factor(component,
                       levels = c("Preparation", "Training", "On-the-fly training")),
    sample_id = as.numeric(Number),
    x_pos = case_when(
      workflow == "Memmap workflow" ~ sample_id - 0.18,
      workflow == "dmpp on-the-fly" ~ sample_id + 0.18))

label.df <- data %>%
  mutate(sample_id = as.numeric(Number),
         xstart = sample_id - 0.18, xend = sample_id + 0.18,
         y_top = pmax(memmap_total, dmpp_total) / 3600,
         label = paste0(round(speedup, 1), "x reduction"))

label_samples <- c("10K", "50K", "93K")

label.df <- label.df %>%
  filter(Number %in% label_samples)

y_offset <- max(label.df$y_top) * 0.045

label.df <- label.df %>%
  mutate(y_bracket = y_top + y_offset,
         y_tick = y_bracket - y_offset * 0.25,
         y_text = y_bracket + y_offset * 0.85)

component_cols <- c("Preparation" = "#b59ad9", "Training" = "#6f73c8",
                    "On-the-fly training" = "#3f88c5")

p <- ggplot(plot.df, aes(x = x_pos, y = time_hour, fill = component)) +
  geom_col(width = 0.32, color = "black", linewidth = 0.25)+
  
  geom_segment(data = label.df,
               aes(x = xstart, xend = xend, y = y_bracket, yend = y_bracket),
               inherit.aes = FALSE, linewidth = 0.45)+
  geom_text(data = label.df,
            aes(x = (xstart + xend) / 2, y = y_text, label = label),
            inherit.aes = FALSE, size = 3.1, fontface = "bold")+
  scale_fill_manual(values = component_cols) +
  scale_x_continuous(
    breaks = seq_along(levels(data$Number)),
    labels = levels(data$Number),
    expand = expansion(mult = c(0.04, 0.04)))+
  scale_y_continuous(
    labels = label_number(accuracy = 0.5),
    expand = expansion(mult = c(0, 0.22))
  )+
  labs(x = "Number of training samples",
       y = "Time to complete the first 1K training steps (h)",
       fill = NULL)+
  theme_classic(base_size = 11) +
  theme(legend.position = "top", legend.justification = "center",
        legend.key.width = unit(0.8, "cm"),
        legend.text = element_text(color = "black", size = 10),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11),
        axis.line = element_line(linewidth = 0.45),
        axis.ticks = element_line(linewidth = 0.45))

pdf("fig_s7/fig_s7a.pdf", width = 6, height = 4)
p
dev.off()

speedup.df <- data %>%
  transmute(Number,
            memmap_workflow_h = memmap_total / 3600,
            dmpp_on_fly_h = dmpp_total / 3600,
            speedup = memmap_total / dmpp_total,
            training_overhead_percent = (train_dmpp / train_memmap - 1) * 100)

speedup.df

###
data <- read.table("fig_s7/data/gpu_wait_long.tsv", header = T)

long.df <- data %>%
  mutate(
    sample = factor(sample, levels = c("1K", "5K", "10K", "50K", "93K")),
    method = factor(method, levels = c("memmap", "dmpp")),
    util_bin = case_when(
      gpu_util < 50 ~ "<50%",
      gpu_util < 95 ~ "50–95%",
      TRUE ~ "≥95%"),
    util_bin = factor(util_bin, levels = rev(c("<50%", "50–95%", "≥95%")))
  )

util.df <- long.df %>%
  count(sample, method, util_bin) %>%
  group_by(sample, method) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup()

util_cols <- c("<50%" = "#d9d9d9", "50–95%" = "#9ecae1", "≥95%" = "#3182bd")

util.df <- subset(util.df, util_bin !="≥95%")

p_util_bin <- ggplot(util.df,
                     aes(x = sample, y = fraction, fill = util_bin)) +
  geom_col(width = 0.68, color = "black", linewidth = 0.25) +
  facet_wrap(~method, nrow = 1) +
  scale_fill_manual(values = util_cols) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(x = "Number of training samples",
       y = "Fraction of steps with GPU utilization <95%",
       fill = "GPU utilization") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 10),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 11),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.45)
  )

pdf("fig_s7/fig_s7b.pdf", width = 6, height = 4)
p_util_bin
dev.off()

