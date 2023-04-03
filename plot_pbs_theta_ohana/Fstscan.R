setwd("~/Documents/work/PhD/Tasks/Task4_POPGEN_final/pbs/")
library(ggplot2)
library(qqman)
library(ggpubr)
library(ggrepel)

#par(mfrow=c(1,1))
slidingFst <-read.csv("vr_vvEU_vvNA_qc_ind66cov1pval5_snps_50kb10kb.txt", sep = '\t',header = F)
colnames(slidingFst) <- c("chr", "pos_start", "pos_end", "fst01_w", "fst02_w", "fst12_w", "fst01_uw", "fst02_uw", "fst12_uw","PBS0_w", "PBS1_w", "PBS2_w", "PBS0_uw", "PBS1_uw", "PBS2_uw", "nsites", "genes")
#colnames(slidingFst) <- c("chr", "pos_start", "pos_end", "fst01_w", "fst02_w", "fst12_w", 
                        #  "PBS0_w", "PBS1_w", "PBS2_w", "nsites")
#slidingFst<-slidingFst[slidingFst$nsites>10,] 

scaffold_order<-read.delim("scaffolds_ordered.txt", header = FALSE)
scaffold_order$chrom<-scaffold_order$V4
scaffold_order$ranking<-c(1:length(scaffold_order$V4))


pbsdata<-data.frame(chrom=slidingFst$chr, start=slidingFst$pos_start, end=slidingFst$pos_end, 
                    PBS0_w=slidingFst$PBS0_w,
                    PBS1_w=slidingFst$PBS1_w,
                    PBS2_w=slidingFst$PBS2_w,
                    fst01_w=slidingFst$fst01_w,
                    fst02_w=slidingFst$fst02_w,
                    fst12_w=slidingFst$fst12_w,
                    nsites=slidingFst$nsites,
                    genes=slidingFst$genes)

test<-merge(scaffold_order, pbsdata, by="chrom", all.y=TRUE)
test<-na.omit(test)
test2 <-test[
  with(test, order(test$V3, test$ranking, test$start)),
  ]

test2<-subset(test2, V1!=39)



test2$Window<-c(1:length(test2$V3))


##### BARPLOT DIFFERENTIATION
test2$target <- "Genome-wide"
test2$target[test2$V4 == "scaffold_3" & test2$start >=30530000 & test2$start <= 31799999] <- "tmp1" ### to remove CRADD
test2<-subset(test2, target!="tmp1")
test2$target[test2$V4 == "scaffold_155"] <- "chr13 affected"
test2$target[test2$V4 == "scaffold_47"] <- "chr13 affected"
test2$target[test2$V4 == "scaffold_200"] <- "chr13 affected"
test2$target[test2$V4 == "scaffold_126"] <- "chr13 affected"
test2$target[test2$V4 == "scaffold_264"] <- "chr13 affected"
#test2<-subset(test2, Window>=103117)

#test2$target[test2$V4 == "scaffold_99" & test2$start >= 4738670  & test2$start <=4776939] <- "SLC6A16"
#test2$target[test2$V4 == "scaffold_39" & test2$start >= 10357249 & test2$start <=10402470] <- "HPS5"
#test2$target[test2$V4 == "scaffold_127" & test2$start >= 648938 & test2$start <=744024] <- "MGAM"
#test2$target[test2$V4 == "scaffold_36" & test2$start >= 7951199  & test2$start <=8021476] <- "NPR3"

#vp_vr_vv_ind66cov1pval5_allvulpes_50kb10kb.txt
fst01<-ggplot(test2, aes(x=target, y=fst01_w, color=target)) + xlab("") + ylab("Fst") +
  geom_boxplot()  + scale_color_manual(values=c("#E69F00", "#999999","#999999", "#999999","#999999","#999999")) +
  theme(legend.position="None", axis.title=element_text(size=12)) +
  scale_alpha_manual(values=c(1,0.1)) +scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme(legend.position = "none") + theme(text = element_text(size=12)) + ggtitle("Pale vs Rueppell's") + theme(plot.title = element_text(size = 8, face = "bold")) 
fst02<-ggplot(test2, aes(x=target, y=fst02_w, color=target)) + xlab("") + ylab("Fst") +
  geom_boxplot()  + scale_color_manual(values=c("#E69F00", "#999999","#999999","#999999","#999999","#999999")) +
  theme(legend.position="None", axis.title=element_text(size=12)) +
  scale_alpha_manual(values=c(1,0.1)) + scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme(legend.position = "none") + theme(text = element_text(size=12)) + ggtitle("Pale vs Red fox") + theme(plot.title = element_text(size = 8, face = "bold")) 
fst12<-ggplot(test2, aes(x=target, y=fst12_w, color=target)) + xlab("") + ylab("Fst") +
  geom_boxplot()  + scale_color_manual(values=c("#E69F00", "#999999","#999999","#999999","#999999","#999999")) +
  theme(legend.position="None", axis.title=element_text(size=12)) +
  scale_alpha_manual(values=c(1,0.1)) + scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme(legend.position = "none") + theme(text = element_text(size=12)) + ggtitle("Rueppell's vs Red fox") + theme(plot.title = element_text(size = 8, face = "bold")) 

gene<-ggarrange(fst01, fst02, fst12, ncol= 3) 
ggsave("chr13_affected_boxplot_vp_vr_vv.tiff", gene,"tiff", dpi=900, width=8,height=3, units="in")  





### MANHATTANS
thresh_pbs0<-quantile(test2$PBS0_w, 0.9995, na.rm = TRUE)
thresh_pbs1<-quantile(test2$PBS1_w, 0.9995, na.rm = TRUE)
thresh_pbs2<-quantile(test2$PBS2_w, 0.9995, na.rm = TRUE)
thresh_fst12<-quantile(test2$fst12_w, 0.9995, na.rm = TRUE)
thresh_fst01<-quantile(test2$fst01_w, 0.9995, na.rm = TRUE)
thresh_fst02<-quantile(test2$fst02_w, 0.9995, na.rm = TRUE)

cbPalette_gene<- c( "#d5d5d4")
cbPalette_vv<- c( "#d5d5d4","#e69f00")
cbPalette_vvNA<- c( "#ff7f00","#d5d5d4")
cbPalette_vvEU<- c( "#33a02c", "#d5d5d4")
cbPalette_vR<- c("#e31a1c", '#d5d5d4')
cbPalette_vZ<- c("#6C1717", '#d5d5d4')


#test2<-subset(test2, V1==13)


#cbPalette <- c("#D55E00", "#CC79A7","#009E73", "#F0E442", "#0072B2", "#999999", "#E69F00", "#56B4E9","#56B4E9" )

plotfst<-function(small1, stat, poplegend, threshold, cbPalette){
   n_of_chr<-length(table(small1$V3))
  ggplot(small1, aes(x=Window, y=stat, color=V3)) + ylab(poplegend) + xlab("") +
    geom_point(size=.1)+
    theme_classic(base_size=20)+ ylim(0,NA) + # xlim(7.5e+06,8.5e+06) +
    theme(legend.position="none") +
    theme(axis.text.x=element_blank()) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_color_manual(values=rep(cbPalette,n_of_chr))  +
    geom_hline(yintercept =threshold) 
    #geom_hline(yintercept = 1-threshold,linetype="dashed")
}



#test2$target[test2$V4 == "scaffold_47" & test2$start >= 3595567  & test2$start <= 3866575] <- "FRYL"
#pbs0<-plotfst(test2, test2$PBS0_w, "PBS", thresh_pbs0, cbPalette_vR) + theme(text = element_text(size=12)) +  ggtitle("Rueppell") + xlab("Chromosome") + theme(plot.title = element_text(size = 12, face = "bold")) 
pbs1<-plotfst(test2, test2$PBS1_w, "PBS", thresh_pbs1, cbPalette_vvEU) + theme(text = element_text(size=12)) +  ggtitle("Eurasian Red fox") + xlab("") + theme(plot.title = element_text(size = 12, face = "bold")) 
pbs1_nointro<-plotfst(test2, test2$PBS1_w, "PBS", thresh_pbs1, cbPalette_vvEU) + theme(text = element_text(size=12)) +  ggtitle("Eurasian Red fox") + xlab("") + theme(plot.title = element_text(size = 12, face = "bold")) 

pbs2<-plotfst(test2, test2$PBS2_w, "PBS", thresh_pbs2, cbPalette_vvNA) + theme(text = element_text(size=12)) + ggtitle("") + theme(plot.title = element_text(size = 12, face = "bold")) 


fst01<-plotfst(test2, test2$fst01_w, "Fst", thresh_fst01, cbPalette_gene) + theme(text = element_text(size=12)) + ggtitle("Fennec vs Rueppell's") + theme(plot.title = element_text(size = 12, face = "bold")) 
fst02<-plotfst(test2, test2$fst02_w, "Fst", thresh_fst02, cbPalette_gene)  + theme(text = element_text(size=12)) + ggtitle("Fennec vs Red fox") + theme(plot.title = element_text(size = 12, face = "bold")) 
fst12<-plotfst(test2, test2$fst12_w, "Fst", thresh_fst12, cbPalette_gene) + xlab("") + theme(text = element_text(size=12)) + ggtitle("Rueppell's vs Red fox") + theme(plot.title = element_text(size = 12, face = "bold")) 
gene<-ggarrange(fst01, fst02, fst12, nrow = 3)
ggsave("chr13_Fst_vz_vr_vvEU.tiff", gene,"tiff", dpi=900, width=4,height=7, units="in")  

pbs_plot<-ggarrange(pbs1, alpha_vvEU,  llratio_delta_EU, nrow = 3) 
ggsave("pbs_ohana_alpha_vvEU_plot.tiff", pbs_plot,"tiff", dpi=900, width=10,height=8, units="in")

Figure4BC<-ggarrange(pbs2,  alpha_vvNA,  ncol=1) 
ggsave("Figure4BC.pdf", Figure4BC,"pdf", dpi=300, width=10,height=5, units="in", useDingbats = TRUE)
ggsave("Figure4BC.svg", Figure4BC,"svg", dpi=300, width=10,height=5, units="in")
ggsave("Figure4BC.eps", Figure4BC,"eps", dpi=300, width=10,height=5, units="in")
ggsave("Figure4BC.tiff", Figure4BC,dpi=300, width=10,height=5, units="in")

outliers_pbs0<-test2[test2$PBS0_w>thresh_pbs0,] 
outliers_pbs1<-test2[test2$PBS1_w>thresh_pbs1,] 
outliers_pbs2<-test2[test2$PBS2_w>thresh_pbs2,] 
outliers_fst12<-test2[test2$fst12_w>thresh_fst12,] 
write.csv(outliers_pbs1, "~/Documents/work/PhD/Tasks/Task4_POPGEN_final/pbs/vvEU_99.95perc.csv", row.names = FALSE)



       












##### MANHATTAN STYLE
test2$midpos<-test2$start + ((test2$end-test2$start)/2)
manhattanqq<-data.frame(chr=test2$V1, scaffold=test2$V4, start=test2$start , end=test2$end , bp=test2$midpos, PBS0=test2$PBS0_w,  PBS1=test2$PBS1_w, PBS2=test2$PBS2_w, fst12=test2$fst12_w, snp=test2$nsites)
par(mfrow=c(1,1))
par(mar=c(3.9,4,2,2))
#PBSa<-manhattan(manhattanqq, chr="chr", bp= "bp", p= "PBS0", snp="snp", suggestiveline = F, genomewideline = F, logp=FALSE, annotateTop=TRUE, cex = 0.6, col=c( "#D55E00", 'grey'), chrlabs = c(1:38, "X"), ylab="PBS Rueppell") #, #ylim = c(-1, 0))
#Fstbc<-manhattan(manhattanqq, chr="chr", bp= "bp", p= "fst12", snp="snp", suggestiveline = F, genomewideline = F, logp=FALSE, annotateTop=TRUE, cex = 0.6, col=c( "#D55E00", 'grey'), chrlabs = c(1:38), ylab="Fst") #, #ylim = c(-15, 0))
PBSb<-manhattan(manhattanqq, chr="chr", bp= "bp", p= "PBS2", snp="snp", suggestiveline = F, genomewideline = F, logp=FALSE, annotateTop=TRUE, cex = 0.4, col=c( "#ff7f00","#d5d5d4"), chrlabs = c(1:38), ylab="PBS North African Red foxes", ylim = c(0, 1.2))
PBSc<-manhattan(manhattanqq, chr="chr", bp= "bp", p= "PBS1", snp="snp", suggestiveline = F, genomewideline = F, logp=FALSE, annotateTop=TRUE, cex = 0.4, col=c( "#33a02c", "#d5d5d4"), chrlabs = c(1:38), ylab="PBS Eurasian Red foxes", ylim = c(0, 1.2))


outliers_EU<-manhattanqq[manhattanqq$PBS1 > quantile(manhattanqq$PBS2, 0.999, na.rm=TRUE),] 
outliers_NA<-manhattanqq[manhattanqq$PBS2 > quantile(manhattanqq$PBS1, 0.999, na.rm=TRUE),] 
