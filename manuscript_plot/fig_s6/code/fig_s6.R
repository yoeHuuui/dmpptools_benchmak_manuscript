library(ggplot2)
library(ggsci)

repeat_region <- read.table("fig_s6/data/file6.Col-PEK1.5_allrepeats.gff3")

processing_data <- function(x, y){
  cg <- read.table(paste0("fig_s6/data/",x,".allrepeat.CG.tsv"),
                   header = T)
  chg <- read.table(paste0("fig_s6/data/",x,".allrepeat.CHG.tsv"),
                    header = T)
  chh <- read.table(paste0("fig_s6/data/",x,".allrepeat.CHH.tsv"),
                    header = T)
  cg$group <- repeat_region$V3
  chg$group <- repeat_region$V3
  chh$group <- repeat_region$V3
  c.meth <- rbind(cg, chg, chh)
  c.meth$sample <- y
  return(c.meth)
}

SRX12151071 <- processing_data("SRX12151071", "WT")
SRX12151072 <- processing_data("SRX12151072", "ddcc")
SRX12151073 <- processing_data("SRX12151073", "met1-3")

count.df <- as.data.frame(table(repeat_region$V3))
index <- as.character(count.df[count.df$Freq >= 500, 1])

c.meth <- rbind(SRX12151071, SRX12151072, SRX12151073)


c.meth <- c.meth[c.meth$group %in% index, ]
c.meth <- subset(c.meth, group  != "Low_complexity")
c.meth <- subset(c.meth, group  != "SSR")
c.meth <- subset(c.meth, group  != "tandem_repeat")
c.meth <- na.omit(c.meth)

c.meth$sample <- factor(c.meth$sample, 
                        levels = c("WT", "ddcc", "met1-3"))

sample_cols <- c("WT"= "#4C78A8", "ddcc"= "#7B6FD0", "met1-3"="#B07AA1")

p <- ggplot(c.meth, aes(x=group, y=mean_meth, fill=sample))+
  geom_boxplot()+
  scale_fill_manual(values = sample_cols)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top")+
  xlab("")+ylab("DNA methylation level")+
  labs(fill="")+
  facet_wrap(.~context, ncol = 1)

pdf("fig_s6/fig_s6.pdf", width = 8, height = 6) 
p
dev.off()
