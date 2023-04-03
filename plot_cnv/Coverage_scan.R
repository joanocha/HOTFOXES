library(ggplot2)
library(qqman)
library(ggpubr)
library(ggrepel)

coverage_matrix <-read.csv("normalized_cov_old100kb.tsv", sep = '\t',header = TRUE, row.names=1)
sample_info<-read.delim2('Phd_genomic_dataset.txt', sep="\t", header=TRUE,row.names=1)
coverage_matrix2<-merge(sample_info,coverage_matrix, by='row.names', all=TRUE)
coverage_matrix3<-na.omit(coverage_matrix2)
colors <- c('#33a02c','#ff7f00', '#e31a1c','#6C1717', '#fb9a99')
coverage_matrix3$Group <- factor(coverage_matrix3$Group, levels = c("Vulpes vulpes (Eurasia)", "Vulpes vulpes (North Africa)", "Vulpes rueppellii (North Africa)", "Vulpes zerda (North Africa)", "Vulpes pallida (North Africa)"))
CRADDbar<-ggplot(data=coverage_matrix3, aes(x=NAME, y=scaffold_3.30500000.30600000, color=Group)) + 
  geom_bar(stat="identity", fill="white") + scale_color_manual(values = colors) + 
  ylab("Normalized coverage") + xlab("") + ylim(0,1) + 
  ggtitle("CRADD") + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("CRADDbar.pdf", CRADDbar, dpi=300, width=5.2,height=3, units="in")