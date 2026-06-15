library(ggplot2)
library(tidyverse)
library(ggsci)
library(patchwork)
library(ggbreak)

df <- read.table("fig_s3/data/filesize_compare.tsv", header = T)
df <- df[,c(2,3,11,12)]
colnames(df)[3] <- "mdmpp"
colnames(df)[4] <- "single dmpp"
df <- reshape2::melt(df, c("n_samples", "rep"))
df$value <- df$value/1024
pd <- position_dodge(width = 0.75)

format_cols <- c("single dmpp" = "#3f88c5", "mdmpp" = "#b59ad9")

p <- ggplot()+
  stat_summary(data = df,
               aes(x = factor(n_samples), y = value, fill = variable, group = variable),
               fun = mean, geom = "bar", position = pd, width = 0.65)+
  stat_summary(data = df,
    aes(x = factor(n_samples), y = value, group = variable),
    fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar",
    position = pd, width = 0.15, linewidth = 0.4, show.legend = FALSE)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = format_cols)+
  theme_classic(base_size = 11)+
  theme(legend.position = "top")+
  xlab("Samples in mdmpp") +
  ylab("File size (GB)") +
  labs(fill = "")

pdf("fig_s3/fig_s3.pdf", width = 4, height = 3)
p
dev.off()