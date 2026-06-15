library(tidyverse)
library(ggplot2)
library(ggrepel)
library(readxl)

root <- "fig_s5/data/final_version_log_data/"

method_levels <- c("dmpp", "dmpp-dedup", "wgbstools", "ballc")
query_type_levels_all <- c("promoter", "genomic_window")
sample_order <- c("ENCFF790EEU", "ENCFF072EDU", "ENCFF980QLU")

method_cols <- c("dmpp" = "#1b9e77", "dmpp-dedup" = "#66c2a5",
                 "wgbstools" = "#d95f02", "ballc" = "#7570b3")

safe_max <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(0)
  max(x)
}

fmt_time <- function(x) {
  case_when(is.na(x) ~ "NA", x < 60 ~ sprintf("%.0fs", x),
            x < 3600 ~ sprintf("%.1fm", x / 60),
            TRUE ~ sprintf("%.1fh", x / 3600))
}

fmt_mem <- function(x) {
  sprintf("%.2fGB", x / 1024 / 1024)
}

summary_files <- list.files(
  root,
  pattern = "^summary.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

read_one <- function(f) {
  x <- read.table(
    f,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    na.strings = c("NA", "NaN", "")
  )
  
  if (!"sample" %in% colnames(x)) {
    sample_name <- basename(dirname(f))
    x$sample <- sample_name
  }
  
  if (!"input_route" %in% colnames(x)) {
    x$input_route <- NA_character_
  }
  
  x
}

data <- map_dfr(summary_files, read_one) %>%
  filter(status == "OK") %>%
  mutate(
    sample = as.character(sample),
    input_route = as.character(input_route),
    cpu_sec = user_sec + system_sec,
    time_sec = cpu_sec
  ) %>%
  filter(!is.na(time_sec))

# External dedup/index steps are shared by methods that use pre-deduplicated BAM.
shared_dedup <- data %>%
  filter(task %in% c("dedup_bam_sambamba", "index_dedup_bam")) %>%
  group_by(sample) %>%
  summarise(
    #shared_dedup_sec = sum(wall_sec, na.rm = TRUE),
    shared_dedup_sec = sum(time_sec, na.rm = TRUE),
    shared_dedup_max_rss_kb = safe_max(max_rss_kb),
    .groups = "drop"
  )

calling_map <- tribble(
  ~task,                                      ~method,
  "prepare_genome_ballc_cmeta",              "ballc",
  "call_methylation_allc_dedup",             "ballc",
  "call_methylation_ballc_from_allc_dedup",  "ballc",
  
  "prepare_genome_wgbstools",                "wgbstools",
  "call_methylation_wgbstools_dedup",        "wgbstools",
  
  "call_methylation_dmpp_bam",               "dmpp-dedup",
  "call_methylation_dmpp_dedup",             "dmpp"
)

calling_specific <- data %>%
  inner_join(calling_map, by = "task") %>%
  group_by(sample, method) %>%
  summarise(
    #specific_calling_sec = sum(wall_sec, na.rm = TRUE),
    specific_calling_sec = sum(time_sec, na.rm = TRUE),
    specific_calling_max_rss_kb = safe_max(max_rss_kb),
    .groups = "drop"
  )

methods_need_external_dedup <- c("ballc", "wgbstools", "dmpp")

calling_sum <- crossing(
  sample = unique(data$sample),
  method = method_levels
) %>%
  left_join(calling_specific, by = c("sample", "method")) %>%
  left_join(shared_dedup, by = "sample") %>%
  mutate(
    specific_calling_sec = replace_na(specific_calling_sec, 0),
    specific_calling_max_rss_kb = replace_na(specific_calling_max_rss_kb, 0),
    shared_dedup_sec = replace_na(shared_dedup_sec, 0),
    shared_dedup_max_rss_kb = replace_na(shared_dedup_max_rss_kb, 0),
    
    external_dedup_sec = if_else(method %in% methods_need_external_dedup,
                                 shared_dedup_sec, 0),
    external_dedup_max_rss_kb = if_else(method %in% methods_need_external_dedup,
                                        shared_dedup_max_rss_kb, 0),
    
    calling_sec = specific_calling_sec + external_dedup_sec,
    calling_hour = calling_sec / 3600,
    calling_max_rss_kb = pmax(specific_calling_max_rss_kb,
                              external_dedup_max_rss_kb)
  ) %>%
  filter(calling_sec > 0)

# =========================
# 3. Query benchmark
# =========================

query_actual <- data %>%
  filter(step %in% c("calculate_promoter_methylation",
                     "calculate_genomic_window_methylation")) %>%
  mutate(
    method = case_when(
      str_detect(task, "ballc") ~ "ballc",
      str_detect(task, "wgbstools") ~ "wgbstools",
      str_detect(task, "dmpp_bam") ~ "dmpp-dedup",
      str_detect(task, "dmpp_dedup") ~ "dmpp",
      TRUE ~ NA_character_
    ),
    query_type = case_when(
      str_detect(task, "promoter") ~ "promoter",
      str_detect(task, "genomic_window|window") ~ "genomic_window",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(method %in% method_levels, !is.na(query_type)) %>%
  group_by(sample, method, query_type) %>%
  summarise(
    #query_region_sec = mean(wall_sec, na.rm = TRUE),
    query_region_sec = mean(time_sec, na.rm = TRUE),
    query_region_max_rss_kb = safe_max(max_rss_kb),
    n_query_runs = n(),
    .groups = "drop"
  )

bed_prep_map <- tribble(
  ~task,                                  ~method,      ~query_type,
  "prepare_promoter_bed_wgbstools",       "wgbstools",  "promoter",
  "prepare_window_bed_wgbstools",         "wgbstools",  "genomic_window",
  "prepare_genomic_window_bed_wgbstools", "wgbstools",  "genomic_window"
)

query_bed <- data %>%
  inner_join(bed_prep_map, by = "task") %>%
  group_by(sample, method, query_type) %>%
  summarise(
    #bed_sec = mean(wall_sec, na.rm = TRUE),
    bed_sec = mean(time_sec, na.rm = TRUE),
    bed_max_rss_kb = safe_max(max_rss_kb),
    .groups = "drop"
  )

query_sum <- query_actual %>%
  full_join(query_bed, by = c("sample", "method", "query_type")) %>%
  mutate(
    bed_sec = replace_na(bed_sec, 0),
    query_region_sec = replace_na(query_region_sec, 0),
    bed_max_rss_kb = replace_na(bed_max_rss_kb, 0),
    query_region_max_rss_kb = replace_na(query_region_max_rss_kb, 0),
    
    query_total_sec = bed_sec + query_region_sec,
    query_total_hour = query_total_sec / 3600,
    query_total_max_rss_kb = pmax(bed_max_rss_kb, query_region_max_rss_kb)
  )

query_type_levels <- query_type_levels_all[
  query_type_levels_all %in% unique(query_sum$query_type)
]

e2e_sum <- calling_sum %>%
  crossing(query_type = query_type_levels) %>%
  left_join(query_sum, by = c("sample", "method", "query_type")) %>%
  mutate(
    query_total_sec = replace_na(query_total_sec, 0),
    query_total_hour = query_total_sec / 3600,
    query_total_max_rss_kb = replace_na(query_total_max_rss_kb, 0),
    total_sec = calling_sec + query_total_sec,
    total_hour = total_sec / 3600,
    total_max_rss_kb = pmax(calling_max_rss_kb, query_total_max_rss_kb),
    total_max_rss_gb = total_max_rss_kb / 1024 / 1024,
    
    method = factor(method, levels = method_levels),
    query_type = factor(query_type, levels = query_type_levels)
  ) %>%
  group_by(sample, query_type) %>%
  mutate(
    dmpp_total_sec = total_sec[method == "dmpp"][1],
    speedup_vs_dmpp = total_sec / dmpp_total_sec
  ) %>%
  ungroup()

e2e_sum$sample <- factor(e2e_sum$sample, levels = sample_order)

theme_e2e <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(
      fill = "#F5F5F2", color = "#C9C9C3", linewidth = 0.45),
    strip.text = element_text(
      face = "bold", size = 12, color = "#3B3B35",
      margin = margin(t = 5, r = 8, b = 5, l = 8)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right")

plot_query_type <- "promoter"

plot_base <- e2e_sum %>%
  filter(query_type == plot_query_type) %>%
  mutate(method = factor(method, levels = method_levels))

sample_levels <- plot_base %>%
  group_by(sample) %>%
  summarise(med_total = median(total_hour, na.rm = TRUE), 
            .groups = "drop") %>%
  arrange(med_total) %>% pull(sample)

plot_base <- plot_base %>%
  mutate(sample = factor(sample, levels = sample_levels))

pd <- position_dodge(width = 0.55)

time_df <- plot_base %>%
  group_by(sample) %>%
  mutate(
    max_calling = max(calling_hour, na.rm = TRUE),
    max_query = max(query_total_sec, na.rm = TRUE),
    max_total = max(total_hour, na.rm = TRUE)) %>%
  ungroup() %>%
  select(sample, method, calling_hour, query_total_sec, total_hour,
         max_calling, max_query, max_total) %>%
  pivot_longer(
    cols = c(calling_hour, query_total_sec, total_hour),
    names_to = "component", values_to = "value") %>%
  mutate(
    component = recode(
      component, 
      calling_hour = "Calling",
      query_total_sec = "Query",
      total_hour = "Total"),
    component = factor(component, levels = c("Calling", "Query", "Total")),
    relative_value = case_when(
      component == "Calling" ~ value / max_calling,
      component == "Query" ~ value / max_query,
      component == "Total" ~ value / max_total
    ),
    label = case_when(
      component == "Calling" ~ paste0("C: ", sprintf("%.1fh", value)),
      component == "Query" ~ paste0("Q: ", fmt_time(value)),
      component == "Total" ~ paste0("T: ", sprintf("%.1fh", value))
    ),
    y_lab = if_else(
      component == "Query",
      pmax(relative_value - 0.065, 0.03),
      pmin(relative_value + 0.055, 1.12)))

time_df$sample <- factor(time_df$sample, levels = sample_order)

p_time <- ggplot(time_df, aes(component, relative_value)) +
  geom_point(
    aes(color = method, shape = component),
    position = pd,
    size = 3.6,
    stroke = 0.35) +
  geom_text_repel(
    aes(label = label, color = method, y = y_lab),
    position = pd,
    size = 2.8,
    show.legend = FALSE) +
  geom_hline(
    yintercept = 0.2,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey55") +
  facet_grid(. ~ sample, scales = "free_x", space = "free_x") +
  scale_y_continuous(
    name = "Relative time within each sample",
    limits = c(0, 1.16),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
  scale_color_manual(name = "Tool", values = method_cols) +
  labs(x = NULL) +
  theme_e2e

pdf("fig_s5/fig_s5a.pdf", width = 8, height = 4)
p_time
dev.off()


library(tidyverse)
root_dir <- "fig_s5/data/final_version_log_data"

files <- list.files(
  root_dir,
  pattern = "\\.consistency_cor\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

cor_all <- map_dfr(files, function(f) {
  read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
    mutate(dataset = basename(dirname(f)))
})

method_order <- c("ballc", "wgbstools", "dmpp_dedup", "dmpp_bam")

pair_order <- c(
  "BAllCools vs wgbstools",
  "BAllCools vs dmpp-dedup",
  "BAllCools vs dmpp",
  "wgbstools vs dmpp-dedup",
  "wgbstools vs dmpp",
  "dmpp-dedup vs dmpp"
)

plot_df <- cor_all %>%
  mutate(x_id = match(x, method_order), y_id = match(y, method_order)) %>%
  filter(x_id < y_id) %>%
  mutate(x_name = recode(x, "ballc" = "BAllCools", 
                         "wgbstools" = "wgbstools",
                         "dmpp_dedup" = "dmpp-dedup", 
                         "dmpp_bam" = "dmpp"),
         y_name = recode(y, "ballc" = "BAllCools",
                         "wgbstools" = "wgbstools",
                         "dmpp_dedup" = "dmpp-dedup",
                         "dmpp_bam" = "dmpp"),
         pair = paste(x_name, y_name, sep = " vs "),
         pair = factor(pair, levels = rev(pair_order)),
         group = case_when(
           x_name == "BAllCools" | y_name == "BAllCools" ~ "BAllCools",
           x_name == "wgbstools" | y_name == "wgbstools" ~ "wgbstools",
           TRUE ~ "dmpp"),
         label = sprintf("%.3f", pearson_r))

dataset_order <- c("ENCFF790EEU", "ENCFF072EDU",
                   "ENCFF980QLU")

plot_df <- plot_df %>%
  mutate(dataset = factor(dataset, levels = dataset_order))


p <- ggplot(plot_df, aes(x = pearson_r, y = pair)) +
  geom_vline(xintercept = 0.90, color = "grey75", linewidth = 0.4)+
  geom_segment(aes(x = 0.90, xend = pearson_r, yend = pair),
               color = "grey78", linewidth = 0.6)+
  geom_point(aes(color = group), size = 2.8)+
  geom_text(aes(label = label), nudge_x = -0.003, nudge_y = 0.18,
            hjust = 1, vjust = 0, size = 3.0)+
  facet_wrap(~ dataset, nrow = 1) +
  scale_x_continuous(
    limits = c(0.90, 1.005),
    breaks = c(0.90, 0.95, 1.00),
    labels = c("0.90", "0.95", "1.00"),
    expand = expansion(mult = c(0.01, 0.03))) +
  scale_color_manual(
    values = c("BAllCools" = "#9bb7d4", "wgbstools" = "#4c78a8",
               "dmpp" = "#1f3b73")) +
  labs(x = "Pearson's r", y = NULL, color = NULL) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 10, color = "black",
                               angle = 15),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  )

pdf("fig_s5/fig_s5b.pdf", width = 8, height = 4)
p
dev.off()


tab <- read_excel("fig_s5/data/final_version_log_data/CpG_site_data.xlsx")

tab <- tab %>%
  mutate(total_mC = if_else(is.na(total_mC), C_mC + G_mC, total_mC),
         total_cov = if_else(is.na(total_cov), C_cov + G_cov, total_cov),
         ratio = total_mC / total_cov,
         C_site = if_else(is.na(C_mC), "–", paste0(C_mC, "/", C_cov)),
         G_site = if_else(is.na(G_mC), "–", paste0(G_mC, "/", G_cov)),
         CpG_total = paste0(total_mC, "/", total_cov),
         CpG_ratio = sprintf("%.2f", ratio),
         Method = factor(Method,
                         levels = rev(c("Reads", "dmpp", 
                                        "BAllCools", "wgbstools"))))

plot_tab <- tab %>%
  select(site_id, Method, C_site, G_site, CpG_total, CpG_ratio) %>%
  pivot_longer(cols = c(C_site, G_site, CpG_total, CpG_ratio),
               names_to = "Column", values_to = "Value") %>%
  mutate(Column = factor(Column, levels = c("C_site", "G_site", 
                                            "CpG_total", "CpG_ratio")),
         Fill = if_else(Method == "BAllCools", "BAllCools", "other"))

p_tab <- ggplot(plot_tab, aes(x = Column, y = Method)) +
  geom_tile(aes(fill = Fill), color = "grey82",
            linewidth = 0.6, width = 0.96, height = 0.92)+
  geom_text(aes(label = Value), size = 4.5, color = "black")+
  facet_wrap(~ site_id, ncol = 1) +
  scale_fill_manual(values = c("BAllCools" = "#eaf2ff", "other" = "white"))+
  scale_x_discrete(
    labels = c("C site\nmC/covC", "G site\nmC/covC", "CpG total\nmC/covC",
               "CpG\nratio"),position = "top")+
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10.5, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.ticks = element_blank(),
    legend.position = "none")

p_tab

pdf("fig_s5/fig_s5d.pdf", width = 4, height = 6)
p_tab
dev.off()
