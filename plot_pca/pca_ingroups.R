setwd("~/Documents/work/PhD/Tasks/Task4_POPGEN_FINAL/pca/old_npy")
library(RcppCNPy)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(dplyr)

colors <- c('#33a02c','#1f78b4', '#ff7f00', '#e31a1c', '#b15928', "#D55E00", "#CC79A7","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#999999", "#E69F00", "#56B4E9")
plot_pca<-function(matrix,components, n1, n2, e1,e2){
    names<-read.table('ingroups_bams.txt', header=FALSE, stringsAsFactors = FALSE)
    dataset<-read.delim2('Phd_genomic_dataset.txt', sep="\t", header=TRUE)
    names$CODE <- gsub("/space/s1/joana/data/bam_final/", "", names$V1)
    names$CODE <- gsub(".sorted.bam", "", names$CODE)
    covmatrix<-npyLoad(matrix)
    e<-eigen(covmatrix)
    print(e$values[e1]/sum(e$values)*100)
    print(e$values[e2]/sum(e$values)*100)
    pca <- as.data.frame(e$vectors[,components])
    vars<-apply(pca, 2, var)
    vars/sum(vars)
    colnames(pca) <- c("x", "y")
    pca$CODE<-names$CODE
    ingroup_pca <- merge(dataset,pca,by="CODE")
    ingroup_pca$Species1<-factor(ingroup_pca$Species, levels = c("Red fox (Europe)", "Red fox (Asia)", "Red fox (North Africa)", "Rueppell's fox (North Africa)"))
    ggplot(ingroup_pca, aes(x=x, y=y, colour=Species), main=NULL) + geom_point(size=1.5) + 
    xlab(n1) + ylab(n2) + theme_test(base_size = 9) + scale_color_manual(values = colors) 
   # geom_label_repel( aes(label = Country),box.padding=0.15, point.padding = 0.01, label.size = 0.001, segment.size =0.1) 
  }

myplot1<-plot_pca('ingroups_qc_ind65cov2pval6_snps.beagle.gz.cov.npy',  1:2, "PC1", "PC2",1,2) + theme(legend.position = "none")
myplot2<-plot_pca('ingroups_qc_ind65cov2pval6_snps.beagle.gz.cov.npy', 3:4, "PC3", "PC4",3,4) + theme(legend.position = "none")

my_components<-ggarrange(myplot1, myplot2, nrow=1) 

#my_components<-ggarrange(myplot1)
print(my_components)

#+ theme(legend.title = element_text(color="black", size=10), legend.text=element_text(size = 10))

ggsave("pca_ingroups1.pdf", my_components, "pdf", dpi=400, width=4, height=2, units="in")

