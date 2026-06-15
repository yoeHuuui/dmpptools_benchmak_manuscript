library(ggplot2)
library(tidyverse)
library(ggsci)
library(patchwork)
library(ggbreak)
library(scales)

id_bench <-  read.table("fig_3/data/benchmark_id.tsv")

dm_list <- id_bench[id_bench$V1 %in% c("dm_list_outdir", "python_aggregate"), ]
merge_list <- id_bench[!id_bench$V1 %in% c("dm_list_outdir", "python_aggregate"), ]

dm_list %>% group_by(V4,V11) %>%
  summarise(wall=sum(V5), cpu=sum(V6+V7), mem=max(V9)) -> dm_list.summary

merge_list <- merge_list[,c("V1","V4","V11","V5","V6","V7","V9")]
merge_list$cpu <- merge_list$V6+merge_list$V7
merge_list <- merge_list[,c("V4","V11","V5","cpu","V9","V1")]
colnames(merge_list) <- c("V4", "V11", "wall", "cpu", "mem", "method")
dm_list.summary$method  <- "single dmpp"

plot.df <-  rbind(dm_list.summary, merge_list)
#plot.df <-  reshape2::melt(plot.df)

plot.df$V4  <- str_replace_all(plot.df$V4, "\\.rep[0-9]+","")
plot.df$V4  <- str_replace_all(plot.df$V4, "multi.n","")
plot.df[plot.df$V4 =="single","V4"] <- "1"
plot.df[plot.df$method !="single dmpp","method"] <- "mdmpp"
plot.df$V4 <- as.numeric(plot.df$V4)
plot.df$V11 <- str_replace_all(plot.df$V11, "\\.rep[0-9]", "")
plot.df$V11 <- str_replace_all(plot.df$V11, "size", "")
plot.df$V11 <- paste0("n=",plot.df$V11)
plot.df$V11 <- factor(plot.df$V11, levels=paste0("n=", c(16,32,64,128,256,512)))
plot.df$mem <- plot.df$mem/1024

p <- ggplot()+
  stat_summary(data = plot.df,
               aes(x = V4, y = wall, group = interaction(V11, method),
                   color = V11, linetype = method),
               fun = mean, geom = "line", linewidth = 0.7)+
  stat_summary(data = plot.df,
               aes(x = V4, y = wall, group = interaction(V11, method),
                   color = V11),
               fun = mean, geom = "point", size = 1.8)+
  stat_summary(data = plot.df,
               aes(x = V4, y = wall, group = interaction(V11, method),
                   color = V11),
               fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.08, linewidth = 0.4,
               show.legend = FALSE)+
  scale_x_continuous(trans = "log2")+
  scale_color_brewer(palette = "Blues")+
  theme_classic()+
  theme()+
  xlab("Number of queried samples")+ylab("Wall time (s)")+
  labs(color="Samples in mdmpp", linetype="Method")

pdf("fig_3/fig_3b.pdf",width = 6, height = 4)
p
dev.off()

###
df <- read.table("fig_3/data/benchmark_multi_group.tsv")
df %>% group_by(V1,V2,V9) %>%
  summarise(c_time=sum(V5+V6), rss=max(V7)) -> df.plot

df.plot$rss <- df.plot$rss/1024
df.plot$group <- str_split(df.plot$V9, "/", simplify = T)[,4]
df.plot$number <- str_count(df.plot$group, "group")
df.plot$group <- str_remove_all(df.plot$group, "group")


plot_df <- df.plot %>%
  mutate(K = number, method = V1) %>%
  separate_rows(group, sep = "_") %>%
  mutate(G = as.integer(group),
         K = as.integer(K))
plot_df <- subset(plot_df, method != "merge_onepass")

wgbs_tools.df <- read.table("fig_3/data/wgbstools_multigroup_summary.tsv",
                            header = T)

wgbs_tools.df2 <- data.frame(V1="wgbstools", V2=1, V9="NA",
                             c_time=wgbs_tools.df$user_sec+wgbs_tools.df$sys_sec,
                             rss=wgbs_tools.df$max_rss_kb/1024,
                             group=as.numeric(str_split(wgbs_tools.df$group_col,"_", simplify = T)[,2]),
                             number=1, K=1, method="wgbstools",
                             G=as.numeric(str_split(wgbs_tools.df$group_col,"_", simplify = T)[,2]))
plot_df <- rbind(as.data.frame(plot_df), wgbs_tools.df2)

anno <- c("dm_list_long" = "single dmpp + script", "merge_classic" = "mdmpp",
          "wgbstools" = "wgbstools")
plot_df$method <- unname(anno[plot_df$method])

label_df <- plot_df %>%
  filter(method %in% c("single dmpp + script", "mdmpp")) %>%
  group_by(method, K) %>%
  summarise(mean_time = mean(c_time), .groups = "drop") %>%
  group_by(method) %>%
  summarise(
    K1 = mean_time[K == 1],
    K6 = mean_time[K == 6],
    ratio = K6 / K1,
    label = sprintf("+%.1f%%", (ratio - 1) * 100),
    .groups = "drop"
  ) %>%
  mutate(
    x_bracket = 0.08,
    x_text = 1.08
  )


p <- ggplot(plot_df, aes(x = K, y = c_time, color = factor(G),
                         shape = method, group = interaction(G, method)))+
  stat_summary(fun = mean, geom = "line", linewidth = 0.85, alpha = 0.95)+
  stat_summary(fun = mean, geom = "point", size = 2.4, stroke = 0.8)+
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.05, linewidth = 0.35,
               alpha = 0.75)+
  geom_segment(data = label_df,  
               aes(x = 6+x_bracket, xend = 6+x_bracket, y = K1, yend = K6),
               inherit.aes = FALSE,linetype = "dashed", 
               linewidth = 0.55, color = "#2171b5")+
  geom_text(data = label_df,
            aes(x = 6-x_bracket, y = (K1+K6)/2, label = label),
            inherit.aes = FALSE, hjust = 1, vjust = 0, size = 3.8,
            fontface = "bold", color = "#2171b5") +
  scale_color_brewer(palette = "Blues", name = "Number of groups per scheme")+
  scale_x_continuous(
    breaks = sort(unique(plot_df$K)),
    limits = c(0.75, 6.15),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_y_break(c(400, 1175), ticklabels = c(1175, 1200, 1250),
                scales = 0.35, space = 0.10, symbol = "//")+
  labs(x = "Number of grouping schemes computed together",
       y = "Running Time (s)",
       shape="Method")+
  theme_bw()+
  theme(
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.title.x.top = element_blank(),
    axis.line.x.top = element_blank(),
    panel.grid.major = element_line(linewidth = 0.35, colour = "grey88"),
    panel.grid.minor = element_line(linewidth = 0.25, colour = "grey93"),
    legend.position = "right",
  )

pdf("fig_3/fig_3c.pdf", width = 8, height = 5)
p
dev.off()
