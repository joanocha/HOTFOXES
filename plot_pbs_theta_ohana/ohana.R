setwd("~/Documents/work/PhD/Tasks/Task4_POPGEN_final/ohana/")
require(tidyr)
require(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ape)
library(phyclust)

K2_colors <- c('#e31a1c','#ff7f00') 
K3_colors <- c('#ff7f00','#33a02c','#e31a1c') 
K4_colors <- c('#1f78b4','#e31a1c','#ff7f00','#33a02c') 
K5_colors <- c('#fb9a99','#1f78b4','#ff7f00','#e31a1c','#33a02c')
K5_colors <- c('#33a02c','#1f78b4','#fb9a99','#ff7f00', '#e31a1c') 
K6_colors <- c('#ff7f00','#fb9a99','#33a02c','#b2df8a','#1f78b4','#e31a1c') 
K4_colors_allvulpes <- c('#fb9a99','#ff7f00','#e31a1c','#6C1717')        
# '#1f78b4' is blue
# '#e31a1c' is red
# #33a02c is green 
# #ff7f00 is orange 
#'#a6cee3' is light blue 
#''#b2df8a' is light green 
#''#fb9a99' is pink
#'#6a3d9a','#b15928','#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6','#ffff99','#4d6680', '#4d8051', '#804d57', '#604d80', '#806f4d')

K_plot<-function(filename, kcolors, mygrid,xlegend,ktitle) {
    print(filename)
    data <- read.table(filename, skip=1)
    key <- read.delim("ingroups_bams.txt", header=FALSE, sep = "\t", col.names=c('samples', 'pops', 'continent','species','old'))
    data$sample <- key$samples
    data2 <- merge(data, key, by.x="sample", by.y="samples")
    data2 <- data2[order(data2$species, data2$continent, data2$sample),]
    data_long <- gather(data2, -sample, -pops, -continent, -species, -old, key = "Component", value = "Proportion")
    ggplot() +
    geom_col(data=data_long, aes(x=sample, y=Proportion, fill = Component), position = "fill") + 
    scale_fill_manual(values = kcolors,
                      guide=FALSE) +
    facet_grid(mygrid,scales = "free",space="free", switch="x") +
    theme(
      panel.spacing = unit(0, "lines"), 
      strip.background = element_rect(colour="gray"),
      axis.ticks = element_blank(), 
      axis.text.x = element_blank(),
      axis.line.x = element_blank(), 
      strip.text.x = element_blank()) + 
      #strip.text.x = element_text(size = 4, angle=0)) + 
      xlab(xlegend) +
      ggtitle(ktitle) + theme(plot.title = element_text(size = 14)) 
}

K2_plot<-K_plot("subsets/ingroups_K_2_s_12_q", K2_colors, .~species, "", "K=2") 
K3_plot<-K_plot("subsets/ingroups_K_3_s_17_q", K3_colors, .~continent,"", "K=3")
K4_plot<-K_plot("subsets/ingroups_K_4_s_9_q", K4_colors, .~species,"", "K=4")
K5_plot<-K_plot("subsets/allvulpes_K_5_s_26_q", K5_colors, .~species,"", "K=5")
K6_plot<-K_plot("subsets/allvulpes_withArctic_k_6_1percent_subset_q", K6_colors, .~species,"", "K=6")

Figure2B<-ggarrange(K2_plot, K3_plot, K4_plot,  nrow=3, ncol=1, common.legend = TRUE, legend="bottom")
ggsave("Figure2B.tiff", Figure2B, width = 4, height = 4, "tiff", dpi=300, limitsize = TRUE, units = "in")
ggsave("Figure2B.svg", Figure2B, width = 4, height = 4, "svg", dpi=300, limitsize = TRUE, units = "in")

all_K<-ggarrange(K2_plot, K3_plot, K4_plot, K5_plot, K6_plot, nrow=5, common.legend = TRUE, legend="bottom")
print(all_K)
ggsave("all_K.tiff", all_K, width = 8, height = 9, "tiff", dpi=300, limitsize = TRUE, units = "in")

