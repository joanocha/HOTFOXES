setwd("~/Documents/work/PhD/Tasks/Task4_POPGEN_final/ohana/")
library(ggplot2)
library(qqman)
library(ggpubr)
library(ggrepel)

par(mfrow=c(1,1))

#thresh_fst12<-quantile(test2$fst12_w, 0.995, na.rm = TRUE)
cbPalette_gene<- c( "#d5d5d4")
cbPalette_vv<- c( "#d5d5d4","#e69f00")
cbPalette_vvNA<- c( "#ff7f00","#d5d5d4")
cbPalette_vvEU<- c( "#33a02c", "#d5d5d4")
cbPalette_vR<- c("#e31a1c", '#d5d5d4')
cbPalette_vZ<- c("#6C1717", '#d5d5d4')


#cbPalette <- c("#D55E00", "#CC79A7","#009E73", "#F0E442", "#0072B2", "#999999", "#E69F00", "#56B4E9","#56B4E9" )
plotfst<-function(filename, poplegend, cbPalette){
  slidingFst <-read.csv(filename, sep = '\t',header = F)
  #colnames(slidingFst) <- c("chr", "pos_start", "pos_end", "llratio1",  "llratio2", "max", "delta", "genes")
  colnames(slidingFst) <- c("chr", "pos_start", "pos_end", "llratio1",  "genes")
  scaffold_order<-read.delim("scaffolds_ordered.txt", header = FALSE)
  scaffold_order$chrom<-scaffold_order$V4
  scaffold_order$ranking<-c(1:length(scaffold_order$V4))
  #pbsdata<-data.frame(chrom=slidingFst$chr, start=slidingFst$pos_start, end=slidingFst$pos_end, llratio1=slidingFst$llratio1, llratio2=slidingFst$llratio2, max=slidingFst$max, delta=slidingFst$delta, genes=slidingFst$genes)
  pbsdata<-data.frame(chrom=slidingFst$chr, start=slidingFst$pos_start, end=slidingFst$pos_end, delta=slidingFst$llratio1,genes=slidingFst$genes)
  test<-merge(scaffold_order, pbsdata, by="chrom", all.y=TRUE)
  test<-na.omit(test)
  test2 <-test[
    with(test, order(test$V3, test$ranking, test$start)),
  ]
  test2$Window<-c(1:length(test$V3))
  thresh_llratio<-quantile(test2$delta, 0.9995, na.rm = TRUE)
  n_of_chr<-length(table(test2$V3))
  test2<-subset(test2, V1!=39 )
  ggplot(test2, aes(x=Window, y=delta, color=V3)) + ylab(poplegend) + xlab("") +
    geom_point(size=.1)+
    theme_classic(base_size=20)+ 
    theme(legend.position="none") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + ylim(0, NA) + 
    scale_color_manual(values=rep(cbPalette,n_of_chr)) +
   geom_hline(yintercept =thresh_llratio) #linetype="dashed") 
}



llratio_delta_EU<-plotfst("vvEUcs4_vvNAcs1_window_output.txt", "Delta log10(llratio)", cbPalette_vvEU) + theme(text = element_text(size=12)) +  ggtitle("") + xlab("") + theme(plot.title = element_text(size = 12, face = "bold")) 
llratio_delta_NA<-plotfst("vvNAcs1_vvEUcs4_window_output.txt", "Delta log10(llratio)", cbPalette_vvNA) + theme(text = element_text(size=12)) +  ggtitle("") + xlab("") + theme(plot.title = element_text(size = 12, face = "bold")) 
llratio_delta_vR<-plotfst("vr_vv_window_output.top500.txt", "Delta log10(llratio)", cbPalette_vR) + theme(text = element_text(size=12)) +  ggtitle("") + xlab("") + theme(plot.title = element_text(size = 12, face = "bold")) 
llratio_vZ<-plotfst("vz_window_output_top500.txt", "Log10(Likelihood ratio)", cbPalette_vZ) + theme(text = element_text(size=12)) +  ggtitle("") + xlab("Chromosome") + theme(plot.title = element_text(size = 12, face = "bold")) 

test<-llratio_delta_vR$data

#pbs1<-plotfst(test2, test2$PBS1_w, "PBS", thresh_pbs1, cbPalette_vvEU) + theme(text = element_text(size=12)) +  ggtitle("Eurasian Red fox") + xlab("") + theme(plot.title = element_text(size = 12, face = "bold")) 
#pbs2<-plotfst(test2, test2$PBS2_w, "PBS", thresh_pbs2, cbPalette_vvNA) + theme(text = element_text(size=12)) + ggtitle("North African Red fox") + theme(plot.title = element_text(size = 12, face = "bold")) 

Figure5_AB<-ggarrange(llratio_delta_vR, llratio_vZ, nrow = 2) 
ggsave("Figure5AB_dingbats.pdf", Figure5_AB,"pdf", dpi=300, width=9,height=5, units="in", useDingbats = TRUE)
ggsave("Figure5AB.tiff", Figure5_AB, dpi=300, width=9,height=5, units="in")

thresh_llratio_vZ<-quantile(llratio_vZ$data$delta, 0.9995, na.rm = TRUE)
outliers_llratio<-llratio_vZ$data[llratio_vZ$data$delta>thresh_llratio_vZ,] 
write.csv(outliers_llratio, "~/Documents/work/PhD/Tasks/Task4_POPGEN_final/ohana/vZ_99.95_FINAL.csv", row.names = FALSE)

thresh_llratio_delta_vR<-quantile(llratio_delta_vR$data$delta, 0.9995, na.rm = TRUE)
outliers_llratio2<-llratio_delta_vR$data[llratio_delta_vR$data$delta>thresh_llratio_delta_vR,] 
write.csv(outliers_llratio2, "~/Documents/work/PhD/Tasks/Task4_POPGEN_final/ohana/vR_99.95_FINAL.csv", row.names = FALSE)



outliers_pbs1<-test2[test2$PBS1_w>thresh_pbs1,] 
outliers_pbs2<-test2[test2$PBS2_w>thresh_pbs2,] 

