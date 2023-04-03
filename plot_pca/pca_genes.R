setwd("~/Documents/work/PhD/Tasks/Task4_POPGEN/pca")
library(RcppCNPy)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

names<-read.table('ingroups_bams2.txt', header=FALSE, stringsAsFactors = FALSE)
#names<-read.table('vrvvvz_bams.txt', header=FALSE, stringsAsFactors = FALSE)


individual_genes<-function(input, mymatrix, myname){
  names<-read.table(input, header=FALSE, stringsAsFactors = FALSE)
  dataset<-read.table('Phd_genomic_dataset.txt', sep="\t", header=TRUE)
  names$CODE <- gsub("/space/s1/joana/data/bam_final/", "", names$V1)
  names$CODE <- gsub(".sorted.bam", "", names$CODE)
  e1<-eigen(mymatrix)
  e1$values[1]/sum(e1$values)*100
  e1$values[2]/sum(e1$values)*100
  pca1 <- as.data.frame(e1$vectors[,1:2])
  colnames(pca1) <- c("x", "y")
  pca1$CODE<-names$CODE
  vars<-apply(pca1, 2, var)
  vars/sum(vars)
  all_pca <- merge(dataset,pca1,by="CODE")
  #all_pca$Group <- factor(all_pca$Group, levels = c("Vulpes vulpes (Eurasia)", "Vulpes vulpes (North Africa)", "Vulpes rueppellii (North Africa)", "Vulpes zerda (North Africa)"))
  all_pca$Group <- factor(all_pca$Group, levels = c("Vulpes vulpes (Eurasia)", "Vulpes vulpes (North Africa)", "Vulpes rueppellii (North Africa)"))
  colors <- c('#33a02c','#ff7f00', '#e31a1c','#6C1717', '#fb9a99')
  ggplot(all_pca, aes(x=x, y=y, colour=Group)) + geom_point() + xlab("PC1") + ylab("PC2") + theme_test(base_size = 12)+ scale_color_manual(values = colors)+ ggtitle(myname)
  #myplot1<-ggplot(all_pca, aes(x=x, y=y, colour=Group)) + geom_point() + geom_label_repel(label=all_pca$CODE) + xlab("PC1") + ylab("PC2") + theme_test(base_size = 12)+ scale_color_manual(values = colors)+ ggtitle("")
}
covmatrix1<-as.matrix(read.table('SLC6A16_ingroups.covMat'))
covmatrix2<-as.matrix(read.table('MGAM_ingroups.covMat'))
covmatrix3<-as.matrix(read.table('HPS5_ingroups.covMat'))
covmatrix4<-as.matrix(read.table('vrvvvz_chr13_25Mb.covMat'))
covmatrix5<-as.matrix(read.table('CRADD_ingroups.covMat'))

SLC6A16<-individual_genes('ingroups_bams2.txt', covmatrix1, "SLC6A16")
MGAM<-individual_genes('ingroups_bams2.txt', covmatrix2, "MGAM")
HPS5<-individual_genes('ingroups_bams2.txt', covmatrix3, "HPS5")
CRADD<-individual_genes('ingroups_bams2.txt',  covmatrix5, "CRADD")
chr13<-individual_genes('vrvvvz_bams.txt', covmatrix4, "25 Mb region")
#myplot2<-ggplot(all_pca, aes(x=x, y=y, colour=Group)) + geom_point() + xlab("PC3") + ylab("PC4") + theme_test(base_size = 12)+ scale_color_manual(values = colors)
#myplot3<-ggplot(all_pca, aes(x=x, y=y, colour=Group)) + geom_point() + xlab("PC5") + ylab("PC6") + theme_test(base_size = 12)+ scale_color_manual(values = colors)
myplot1<-ggarrange(SLC6A16, MGAM, HPS5, common.legend=TRUE, legend = "bottom", nrow=1, ncol=3) 

ggsave("ED_fig4.pdf", myplot1,"pdf", dpi=300, width=9,height=3, units="in")
ggsave("ED_fig5C.pdf", chr13,"pdf", dpi=300, width=5,height=2, units="in")

ggsave("ED_fig3B.pdf", CRADD,"pdf", dpi=300, width=5.5,height=3, units="in")

names2<-read.table('ingroups_bams.txt', header=FALSE, stringsAsFactors = FALSE)
dataset<-read.table('Phd_genomic_dataset.txt', sep="\t", header=TRUE)
dataset2<-dataset[1:73,]
names2$CODE <- gsub("/space/s1/joana/data/bam_final/", "", names2$V1)
names2$CODE <- gsub(".sorted.bam", "", names2$CODE)
covmatrix2<-npyLoad('SLC6A16_ingroups.covMat')
e2<-eigen(covmatrix2)
pca2 <- as.data.frame(e2$vectors[,1:2])
vars<-apply(pca2, 2, var)
vars/sum(vars)
colnames(pca2) <- c("x", "y")
pca2$CODE<-names2$CODE
ingroup_pca <- merge(dataset2,pca2,by="CODE")
ingroup_pca$Species1<-factor(ingroup_pca$Group, levels = c("Red foxes (Iberia)", "Red foxes (Middle East)", "Red foxes (North Africa)", "Rueppell foxes (North Africa)"))






colors <- c('#1f78b4','#33a02c', '#ff7f00', '#e31a1c', '#b15928', "#D55E00", "#CC79A7","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#999999", "#E69F00", "#56B4E9")
myplot<-ggplot(ingroup_pca, aes(x=x, y=y, colour=Group), main=NULL) + geom_point() + xlab("PC1") + ylab("PC2") + theme_test(base_size = 16) + scale_color_manual(values = colors) + geom_label_repel(aes(label = CODE,box.padding=0.35, point.padding = 0.5))
print(myplot)
ggsave("pca_ingroups.tiff", myplot, "tiff", dpi=1000, width=6 , height=3, units="in")

dataset<-read.table('Phd_genomic_dataset.txt', sep="\t", header=TRUE)
names2<-read.table('vrueppellii_bams.txt', header=FALSE, stringsAsFactors = FALSE)
dataset2<-dataset[1:26,]
names2$CODE <- gsub("/space/s1/joana/data/bam_final/", "", names2$V1)
names2$CODE <- gsub(".sorted.bam", "", names2$CODE)
covmatrix2<-npyLoad('vrueppellii_qc_ind76cov2pval6_snps_mapDOG.cov.npy')
e2<-eigen(covmatrix2)
plot(e2$vectors)
pca2 <- as.data.frame(e2$vectors[,1:2])
vars<-apply(pca2, 2, var)
vars/sum(vars)
colnames(pca2) <- c("x", "y")
pca2$CODE<-names2$CODE
rueppellii_pca <- merge(dataset2,pca2,by="CODE")
#ingroup_pca$Species <-factor(ingroup_pca$Species, levels = c("Vulpes vulpes (Iberia)", "Vulpes vulpes (Middle East)", "Vulpes vulpes (North Africa)", "Vulpes rueppellii (North Africa)"))
colnames(pca2) <- c("x", "y")

colors <- c('#e31a1c','#b15928', "#D55E00", "#CC79A7","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#999999", "#E69F00", "#56B4E9")
myplot<-ggplot(rueppellii_pca, aes(x=x, y=y, colour=Country)) + geom_point(size=0.5) + 
  xlab("PC1") + ylab("PC2") + theme_test(base_size = 16) +
  geom_label_repel(aes(label =CODE,box.padding=0.35, point.padding = 0.5)) 
print(myplot+ scale_color_manual(values = colors))

dataset<-read.table('Phd_genomic_dataset.txt', sep="\t", header=TRUE)
names2<-read.table('vvulpes_bams.txt', header=FALSE, stringsAsFactors = FALSE)
dataset2<-dataset[27:75,]
names2$CODE <- gsub("/space/s1/joana/data/bam_final/", "", names2$V1)
names2$CODE <- gsub(".sorted.bam", "", names2$CODE)
covmatrix2<-npyLoad('vvulpes_qc_ind76cov2pval6_snps_mapDOG.cov.npy')
e2<-eigen(covmatrix2)
plot(e2$vectors)
pca2 <- as.data.frame(e2$vectors[,1:2])
vars<-apply(pca2, 2, var)
vars/sum(vars)
colnames(pca2) <- c("x", "y")
pca2$CODE<-names2$CODE
vvulpes_pca <- merge(dataset2,pca2,by="CODE")
#ingroup_pca$Species <-factor(ingroup_pca$Species, levels = c("Vulpes vulpes (Iberia)", "Vulpes vulpes (Middle East)", "Vulpes vulpes (North Africa)", "Vulpes rueppellii (North Africa)"))
colnames(pca2) <- c("x", "y")

colors <- c('#e31a1c','#b15928', "#D55E00", "#CC79A7","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#999999", "#E69F00", "#56B4E9")
myplot<-ggplot(vvulpes_pca, aes(x=x, y=y, colour=Country)) + geom_point(size=0.5) + 
  xlab("PC1") + ylab("PC2") + theme_test(base_size = 16) +
  geom_label_repel(aes(label =CODE,box.padding=0.35, point.padding = 0.5)) 
print(myplot+ scale_color_manual(values = colors))



dataset<-read.table('Phd_genomic_dataset.txt', sep="\t", header=TRUE)
names2<-read.table('vvEU_bams.txt', header=FALSE, stringsAsFactors = FALSE)
dataset2<-dataset[59:76,]
names2$CODE <- gsub("/space/s1/joana/data/bam_final/", "", names2$V1)
names2$CODE <- gsub(".sorted.bam", "", names2$CODE)
covmatrix2<-npyLoad('vvEU_qc_ind67cov2pval6_snps.cov.npy')
e2<-eigen(covmatrix2)
plot(e2$vectors)
pca2 <- as.data.frame(e2$vectors[,1:2])
vars<-apply(pca2, 2, var)
vars/sum(vars)
colnames(pca2) <- c("x", "y")
pca2$CODE<-names2$CODE
rueppellii_pca <- merge(dataset2,pca2,by="CODE")
#ingroup_pca$Species <-factor(ingroup_pca$Species, levels = c("Vulpes vulpes (Iberia)", "Vulpes vulpes (Middle East)", "Vulpes vulpes (North Africa)", "Vulpes rueppellii (North Africa)"))
colnames(pca2) <- c("x", "y")

colors <- c('#e31a1c','#b15928', "#D55E00", "#CC79A7","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#999999", "#E69F00", "#56B4E9")
myplot<-ggplot(rueppellii_pca, aes(x=x, y=y, colour=Country)) + geom_point(size=0.5) + 
  xlab("PC1") + ylab("PC2") + theme_test(base_size = 16) +
  geom_label_repel(aes(label =Province,box.padding=0.35, point.padding = 0.5)) 
print(myplot+ scale_color_manual(values = colors))








## heatmap / clustering / trees
covmatrix2<-npyLoad('ingroups_qc_snps.cov.npy')
info<-merge(names2, dataset2, by="CODE")
rownames(covmatrix2)<-info$wholegenome
colnames(covmatrix2)<-info$wholegenome
#heat map
heatmap(covmatrix2)
#neighbour joining
plot(ape::nj(covmatrix2))
plot(hclust(dist(covmatrix2), "ave"))