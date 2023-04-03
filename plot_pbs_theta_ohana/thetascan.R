setwd("~/Documents/work/PhD/Tasks/Task4_POPGEN_final/theta/")
library(ggplot2)
library(qqman)
library(ggpubr)
library(ggrepel)


slidingAlpha <-read.delim("output_vvEU_noIntro.txt", header = TRUE)
colnames(slidingAlpha) <- c("chrom", "winstart", "winend", "wincenter","tP_x", "nsites_x", "relative_pi_x", "tP_y", "nsites_y", "relative_pi_y","tP_z", "nsites_z", "relative_pi_z", "relative_total_pi", "relative_alpha", "total_pi", "alpha", "genes")
#slidingAlpha<-slidingAlpha[slidingAlpha$total_pi>0,] 
scaffold_order<-read.delim("scaffolds_ordered.txt", header = FALSE)
scaffold_order$chrom<-scaffold_order$V4
scaffold_order$ranking<-c(1:length(scaffold_order$V4))
theta_data<-data.frame(chrom=slidingAlpha$chrom, start=slidingAlpha$winstart, end=slidingAlpha$winend, sites=slidingAlpha$nsites_x, pi=slidingAlpha$relative_pi_x, alpha=slidingAlpha$alpha, genes=slidingAlpha$genes)
test<-merge(scaffold_order, theta_data, by="chrom", all.y=TRUE)
test$logalpha<- -log10(test$alpha)
test<-na.omit(test)
test0 <-test[with(test, order(test$V3, test$ranking, test$start)),]
test0$Window<-c(1:length(test0$V3))
thresh_alpha0<-quantile(test0$logalpha, 0.9995, na.rm = TRUE)
thresh_alpha1<-quantile(test0$logalpha, 0.9995, na.rm = TRUE)
thresh_pi0<-quantile(test0$pi, 0.99, na.rm = TRUE)
thresh_pi1<-quantile(test0$pi, 0.05, na.rm = TRUE)
#test0<-subset(test0, V1==13)
#test0<-subset(test0, Window>=103255)
test0<-subset(test0, V1!=39)
cbPalette_vR<- c("#e31a1c", "#999999")
cbPalette_vvNA<- c( "#ff7f00","#d5d5d4")
cbPalette_vvEU<- c( "#33a02c","#d5d5d4")
cbPalette_vZ<- c("#6C1717", '#d5d5d4')
cbPalette_vP<- c('#fb9a99', '#d5d5d4')
#cbPalette <- c("#D55E00", "#CC79A7","#009E73", "#F0E442", "#0072B2", "#999999", "#E69F00", "#56B4E9","#56B4E9" )
plottheta<-function(small1, stat, poplegend, threshold1, cbPalette){
  n_of_chr<-length(table(small1$V3))
  ggplot(small1, aes(x=Window, y=stat, color=V3)) + ylab(poplegend) + xlab("") +
    geom_point(size=.1)+
    theme_classic(base_size=20)+ ylim(0.4,NA) +
    theme(legend.position="none") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_color_manual(values=rep(cbPalette,n_of_chr))  +
    geom_hline(yintercept =threshold1) #+ 
   # geom_hline(yintercept =threshold2, linetype="dotted")
}



alpha_vvEU<-plottheta(test0, test0$logalpha, "-log10(alpha)", thresh_alpha0, cbPalette_vvEU) + theme(text = element_text(size=12)) +  ggtitle("")+ xlab("") + theme(plot.title = element_text(size = 12, face = "bold")) 
alpha_vvNA<-plottheta(test0, test0$logalpha, "-log10(alpha)", thresh_alpha0, cbPalette_vvNA) + theme(text = element_text(size=12)) + ggtitle("") + theme(plot.title = element_text(size = 12, face = "bold")) 
#alpha_vR<-plottheta(test0, test0$logalpha, "-log10 alpha", thresh_alpha0, thresh_alpha1,cbPalette_vR) + theme(text = element_text(size=12)) + ggtitle("") + theme(plot.title = element_text(size = 12, face = "bold")) 
#alphaplot<-ggarrange(alpha_vvNA, alpha_vvEU, nrow = 2)
#outliers_alpha0<-test0[test0$logalpha>thresh_alpha0,] 
#write.csv(outliers_alpha0, "~/Documents/work/PhD/Tasks/Task4_POPGEN_final/theta/vvEUnoIntro_99.95.csv", row.names = FALSE)


cbPalette_vR<- c("#e31a1c", "#999999")
cbPalette_vvNA<- c( "#ff7f00","#d5d5d4")
cbPalette_vvEU<- c( "#33a02c","#d5d5d4")
cbPalette_vZ<- c("#6C1717", '#d5d5d4')
cbPalette_vP<- c('#fb9a99', '#d5d5d4')
localtheta <-function(myoutput, poplegend, cbPalette){
  slidingAlpha<-read.delim(myoutput, header = TRUE)
  colnames(slidingAlpha) <- c("chrom", "winstart", "winend", "wincenter","tP_x", "nsites_x", "relative_pi_x", "tP_y", "nsites_y", "relative_pi_y","tP_z", "nsites_z", "relative_pi_z", "relative_total_pi", "relative_alpha", "total_pi", "alpha", "genes")
  slidingAlpha<-slidingAlpha[slidingAlpha$total_pi>0,] 
  scaffold_order<-read.delim("scaffolds_ordered.txt", header = FALSE)
  scaffold_order$chrom<-scaffold_order$V4
  scaffold_order$ranking<-c(1:length(scaffold_order$V4))
  theta_data<-data.frame(chrom=slidingAlpha$chrom, start=slidingAlpha$winstart, end=slidingAlpha$winend, sites=slidingAlpha$nsites_x, pi=slidingAlpha$relative_pi_x, alpha=slidingAlpha$alpha, genes=slidingAlpha$genes)
  test<-merge(scaffold_order, theta_data, by="chrom", all.y=TRUE)
  test$logalpha<- -log10(test$alpha)
  test<-na.omit(test)
  test0 <-test[with(test, order(test$V3, test$ranking, test$start)),]
  test0$Window<-c(1:length(test0$V3))
  thresh_pi0<-quantile(test0$pi, 0.99, na.rm = TRUE)
  test0<-subset(test0, V1==13)
  test0<-subset(test0, Window>=103255)
  n_of_chr<-length(table(test0$V3))
  ggplot(test0, aes(x=Window, y=test0$pi, color=V3)) + ylab(poplegend) + xlab("") +
    geom_line()+
    theme_classic(base_size=20)+ ylim(0,NA) +
    theme(legend.position="none") +
    theme(axis.text.x=element_blank(),
    axis.ticks.x=element_blank()) +
    scale_color_manual(values=rep(cbPalette,n_of_chr))  +
    #geom_hline(yintercept =threshold1) #+ 
    geom_hline(yintercept =thresh_pi0, linetype="dotted")
}
theta_vvEU<-localtheta("output_vvEU.txt",  "π",  cbPalette_vvEU) + theme(text = element_text(size=12)) +  ggtitle("Red fox Eurasia")+ xlab("") + theme(plot.title = element_text(size = 12, face = "bold")) 
theta_vvNA<-localtheta("output_vvNA.txt",  "π", cbPalette_vvNA) + theme(text = element_text(size=12)) + ggtitle("Red fox North Africa") + theme(plot.title = element_text(size = 12, face = "bold")) 
theta_vR<-localtheta("output_vR.txt", "π", cbPalette_vR) + theme(text = element_text(size=12)) + ggtitle("Rueppell's fox") + theme(plot.title = element_text(size = 12, face = "bold")) 
theta_vZ<-localtheta("output_vz.txt", "π", cbPalette_vZ) + theme(text = element_text(size=12)) + ggtitle("Fennec") + theme(plot.title = element_text(size = 12, face = "bold")) 
theta_vP<-localtheta("output_vp.txt", "π", cbPalette_vP) + theme(text = element_text(size=12)) + ggtitle("Pale fox") + theme(plot.title = element_text(size = 12, face = "bold"))  
thetaplot<-ggarrange(theta_vP, theta_vZ,  theta_vvNA,theta_vvEU, theta_vR, nrow = 5)
ggsave("ED_Fig7C.pdf", thetaplot, dpi=300, width=6, height=10, units="in")

data_vZ<-theta_vZ$data
data_vP<-theta_vP$data
data_vR<-theta_vR$data
data_vvNA<-theta_vvNA$data
data_vvEU<-theta_vvEU$data
write.csv(data_vZ, "~/Documents/work/PhD/Tasks/Task4_POPGEN_final/theta/data_vZ.csv", row.names = FALSE)
write.csv(data_vP, "~/Documents/work/PhD/Tasks/Task4_POPGEN_final/theta/data_vP.csv", row.names = FALSE)
write.csv(data_vR, "~/Documents/work/PhD/Tasks/Task4_POPGEN_final/theta/data_vR.csv", row.names = FALSE)
write.csv(data_vvNA, "~/Documents/work/PhD/Tasks/Task4_POPGEN_final/theta/data_vvNA.csv", row.names = FALSE)

#alphaplot<-ggarrange(alpha_vR, alpha_vvNA, alpha_vvEU, nrow = 3)  
#alpha_pbs_vvNA<-ggarrange(alpha_vvNA, pbs2, nrow = 2) 
#alpha_pbs_vR<-ggarrange(alpha_vR, pbs0, nrow = 2) 
#ggsave("alpha_pbs_vvNA.tiff", alpha_pbs_vvNA, "tiff", dpi=800, width=8, height=6, units="in")
#ggsave("alpha_pbs_vR.tiff", alpha_pbs_vR, "tiff", dpi=800, width=8, height=6, units="in")
#ggsave("alpha.tiff", alphaplot,"tiff", dpi=900, width=10,height=7, units="in")
#ggsave("alpha.tiff", alphaplot,"tiff", dpi=300, width=10,height=5, units="in")



#outliers_alpha1<-test1[test1$logalpha>thresh_alpha1,] 
#outliers_alpha2<-test2[test2$logalpha>thresh_alpha2,] 


##### MANHATTAN STYLE
test0<-na.omit(test0)
test2<-na.omit(test2)
manhattanqq0<-data.frame(chr=test0$V1, scaffold=test0$V4, start=test0$start , end=test0$end , bp=test0$center, alpha=test0$alpha,  snp=test0$nsites)
manhattanqq2<-data.frame(chr=test2$V1, scaffold=test2$V4, start=test2$start , end=test2$end , bp=test2$center, alpha=test2$alpha,  snp=test2$nsites)
par(mfrow=c(3,1))
par(mar=c(3.9,4,2,2))
#PBSa<-manhattan(manhattanqq, chr="chr", bp= "bp", p= "PBS0", snp="snp", suggestiveline = F, genomewideline = F, logp=FALSE, annotateTop=TRUE, cex = 0.6, col=c( "#D55E00", 'grey'), chrlabs = c(1:38, "X"), ylab="PBS Rueppell") #, #ylim = c(-1, 0))
Fstbc<-manhattan(manhattanqq, chr="chr", bp= "bp", p= "fst12", snp="snp", suggestiveline = F, genomewideline = F, logp=FALSE, annotateTop=TRUE, cex = 0.6, col=c( "#D55E00", 'grey'), chrlabs = c(1:38, "X"), ylab="Fst") #, #ylim = c(-15, 0))
PBSb<-manhattan(manhattanqq, chr="chr", bp= "bp", p= "PBS1", snp="snp", suggestiveline = F, genomewideline = F, logp=FALSE, annotateTop=TRUE, cex = 0.6, col=c( "#D55E00", 'grey'), chrlabs = c(1:38, "X"), ylab="PBS1") #, #ylim = c(-15, 0))
PBSc<-manhattan(manhattanqq, chr="chr", bp= "bp", p= "PBS2", snp="snp", suggestiveline = F, genomewideline = F, logp=FALSE, annotateTop=TRUE, cex = 0.6, col=c( "#D55E00", 'grey'), chrlabs = c(1:38, "X"), ylab="PBS2")


outliers_EU<-manhattanqq[manhattanqq$PBS1 > quantile(manhattanqq$PBS2, 0.995, na.rm=TRUE),] 
outliers_NA<-manhattanqq[manhattanqq$PBS2 > quantile(manhattanqq$PBS1, 0.995, na.rm=TRUE),] 
