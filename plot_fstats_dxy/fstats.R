setwd("~/Documents/work/PhD/Tasks/Task4_POPGEN_final/genomics_general")
library(ggplot2)
library(qqman)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(hrbrthemes)

dxy <-read.delim("vvNAvvEUvR_arctic_100kb20kb_haplo_fstats_genes.txt", sep = '\t', header = T)
dxy <-read.delim("vVvRvZvP_100kb20kb_haplo_fstats_genes.txt", sep = '\t', header = T)
#dxy<-dxy[dxy$nsites>5000,] 
scaffold_order<-read.delim("scaffolds_ordered.txt", header = FALSE)
scaffold_order$scaffold<-scaffold_order$V4
scaffold_order$ranking<-c(1:length(scaffold_order$V4))
test<-merge(scaffold_order, dxy, by="scaffold", all.y=TRUE)
test<-na.omit(test)
test2 <-test[
    with(test, order(test$V3, test$ranking, test$start)),
  ]
#test2<-subset(test2, V1!=39)
test2$Window<-c(1:length(test2$V2))



#cbPalette <- c( "#D55E00", "#CC79A7","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#999999", "#E69F00", "#56B4E9")
cbPalette <- c( "#D55E00", 'grey')

#test2$target <- "Genome-wide"
#test2$target[test2$V1 == 15 & test2$start >= 30360001 & test2$start <= 31730000] <- "tmp1" ### to remove CRADD
#test2<-subset(test2, target!="tmp1")

thresh_fd<-quantile(test2$fd, 0.999, na.rm = TRUE)
plotdxy<-function(small1, stat, poplegend, threshold1, cbPalette){
  n_of_chr<-length(table(small1$V3))
  ggplot(small1, aes(x=Window, y=stat, color=V3)) + ylab(poplegend) + xlab("") +
    geom_point(size=.4)+
    theme_classic(base_size=20)+ ylim(0,NA) + 
    theme(legend.position="none") +
    theme(axis.text.x=element_blank()) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_color_manual(values=rep(cbPalette,n_of_chr))  +
    geom_hline(yintercept =threshold1, linetype="dashed") 
}

# fstat

plot_fstat_vRvZ<-plotdxy(test2, test2$fd, "fd", thresh_fd, cbPalette) + theme(text = element_text(size=12)) +  ggtitle("") + xlab("Chromosome") + theme(plot.title = element_text(size = 12, face = "bold")) 
plot_fstat_vRvvNA<-plotdxy(test2, test2$fd, "fd", thresh_fd, cbPalette) + theme(text = element_text(size=12)) +  ggtitle("") + xlab("Chromosome") + theme(plot.title = element_text(size = 12, face = "bold")) 


outliers<-test2[ test2$fd>thresh_fd,] 
## Figure 3B top and bottom
ggsave("Figure3B_bottom_dingbat.pdf", plot_fstat_vRvZ, "pdf", dpi=300, width=5 , height=2, units="in", useDingbats = TRUE)
ggsave("Figure3B_bottom.tiff", plot_fstat_vRvZ, dpi=300, width=5 , height=2, units="in")

ggsave("Figure3B_top_dingbat.pdf", plot_fstat_vRvvNA, dpi=300, width=5 , height=2, units="in", useDingbats = TRUE)
ggsave("Figure3B_top.tiff", plot_fstat_vRvvNA, dpi=300, width=5 , height=2, units="in")


Figure3B_topANDbottom<-ggarrange(plot_fstat_vRvvNA, plot_fstat_vRvZ, ncol=1)
ggsave("Figure3B_top&bottom_dingbat.pdf", Figure3B_topANDbottom, "pdf", dpi=300, width=5 , height=4, units="in", useDingbats = TRUE)
ggsave("Figure3B_top&bottom_dingbat.tiff", Figure3B_topANDbottom, dpi=300, width=5 , height=4, units="in")

#write.csv(outliers, "~/Documents/work/PhD/Tasks/Task4_POPGEN_final/genomics_general/vRvZ_fstat_CRADD_genomics_general_99.9th_final.csv", row.names = FALSE)



