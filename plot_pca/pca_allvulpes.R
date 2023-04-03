setwd("~/Documents/work/PhD/Tasks/Task4_POPGEN2/pca")
library(RcppCNPy)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(dplyr)



names<-read.table('all_vulpes_bams.txt', header=FALSE, stringsAsFactors = FALSE)
dataset<-read.table('Phd_genomic_dataset.txt', sep="\t", header=TRUE)
names$CODE <- gsub("/space/s1/joana/data/bam_final/", "", names$V1)
names$CODE <- gsub(".sorted.bam", " ", names$CODE)
covmatrix1<-npyLoad('allvulpes_qc_ind43cov1pval6_snps.cov.npy')
e1<-eigen(covmatrix1)
pca1 <- as.data.frame(e1$vectors[,1:2])
colnames(pca1) <- c("x", "y")
pca1$CODE<-names$CODE
all_pca <- merge(dataset,pca1,by="CODE")
all_pca$Group <- factor(all_pca$Group, levels = c("Vulpes vulpes (Iberia)", "Vulpes vulpes (Middle East)", "Vulpes vulpes (North Africa)", "Vulpes rueppellii (North Africa)", "Vulpes zerda (North Africa)", "Vulpes pallida (North Africa)"))
colors <- c('#1f78b4','#33a02c','#ff7f00', '#e31a1c','#b15928', "#D55E00", "#CC79A7","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#999999", "#E69F00", "#56B4E9")
myplot<-ggplot(all_pca, aes(x=x, y=y, colour=Group)) + geom_point() + xlab("PC1") + ylab("PC2") + theme_test(base_size = 10)
print(myplot + scale_color_manual(values = colors))
