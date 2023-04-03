setwd("~/Documents/work/PhD/Tasks/Task4_POPGEN_final/Dxyz")
library(ggplot2)
library(qqman)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(hrbrthemes)

#NC_006595.3:41200000-63241923 presumably introgressed 
#NC_006595.3:38236300-63241923  presumably introgressed 
#NC_006595.3:59784933-63241923  end part 

barplot_divergence<-function(myfile){
  dxy <-read.delim(myfile, sep = '\t', header = F)
  colnames(dxy) <- c("chrom", "winstart", "winend", "dxy", "nsites", "na", "genes")
  dxy$length<-dxy$winend - dxy$winstart
  dxy<-dxy[dxy$nsites>5000,] 
  dxy$dxy_relative <- dxy$dxy / dxy$nsites
  scaffold_order<-read.delim("scaffolds_ordered_DOG.txt", header = FALSE)
  scaffold_order$chrom<-scaffold_order$V2
  scaffold_order$ranking<-c(1:length(scaffold_order$V2))
  test<-merge(scaffold_order, dxy, by="chrom", all.y=TRUE)
  test<-na.omit(test)
  test2 <-test[
    with(test, order(test$V2, test$ranking, test$winstart)),
  ]
  test2<-subset(test2, V1!=39)
  test2$Window<-c(1:length(test2$V2))
  test2$target <- "Genome-wide"
  test2$target[test2$V1 == 15 & test2$winstart >33769447 & test2$winstart <= 34053174] <- "tmp1" ### to remove CRADD
  #test2$target[test2$V1 == 13 & test2$winstart >38236300 & test2$winstart < 59784933] <- "tmp2"
  test2$target[test2$V1 == 13 & test2$winstart >38236300] <- "chr13:affected"
 # test2$target[test2$V1 == 13 & test2$winstart >= 59784933] <- "chr13:affected end"
  test2<-subset(test2, target!="tmp1")
 # test2<-subset(test2, target!="tmp2")
#  test2<-subset(test2, target!="chr13:affected end")
#  test2$target[test2$V1 == 1 & test2$winstart >=  107149560 & test2$winend <=107184526] <- "SLC6A16"
#  test2$target[test2$V1 == 21 & test2$winstart >= 40726730 & test2$winend <=40764338 ] <- "HPS5"
#  test2$target[test2$V1 == 16 & test2$winstart >= 7053679 & test2$winend <=7269418] <- "MGAM"
#  test2$target[test2$V1 == 4& test2$winstart >=  74799962 & test2$winend <=74862770] <- "NPR3"
  ggplot(test2, aes(x=target, y=dxy_relative, color=target)) + xlab("") + ylab("dxy") +
    geom_boxplot() + scale_color_manual(values=c("#E69F00", "#999999", "#999999", "#999999", "#999999", "#56B4E9")) +
    theme(legend.position="None", axis.title=element_text(size=12)) +
    scale_alpha_manual(values=c(1,0.1)) +
    theme(legend.position = "none") 
}


vr_vvEU<-barplot_divergence("dxy_vr_vvEU_50kb10kb_sites_genes") + ggtitle("Rueppell's vs Red fox Eurasia") 
vr_vvNA<-barplot_divergence("dxy_vr_vvNA_50kb10kb_sites_genes") + ggtitle("Rueppell's and Red fox North Africa") 
vvNA_vvEU<-barplot_divergence("dxy_vvNA_vvEU_50kb10kb_sites_genes") + ggtitle("Red fox North Africa vs Eurasia") 
data_vr_vvNA<-vr_vvNA$data
#SLC6A16<-ggarrange(vr_vvEU, vr_vvNA,vvNA_vvEU, ncol=3)
ggarrange(vr_vvEU, vr_vvNA,vvNA_vvEU, ncol=3)

vr_vZ<-barplot_divergence("dxy_vr_vz_50kb10kb_sites_genes") + ggtitle("Rueppell's and Fennec fox") 
vr_vv<-barplot_divergence("dxy_vr_vv_50kb10kb_sites_genes") + ggtitle("Rueppell's and Red fox")
vZ_vv<-barplot_divergence("dxy_vv_vz_50kb10kb_sites_genes") + ggtitle("Fennec and Red fox")
A<-ggarrange(vr_vZ, vr_vv,vZ_vv, ncol=3)

vP_vR<-barplot_divergence("dxy_vr_vp_50kb10kb_sites_genes") + ggtitle("Pale and Rueppell's fox")
vP_vv<-barplot_divergence("dxy_vv_vp_50kb10kb_sites_genes") + ggtitle("Pale and Red fox")
vP_vZ<-barplot_divergence("dxy_vz_vp_50kb10kb_sites_genes") + ggtitle("Pale and Fenenc fox")
B<-ggarrange(vP_vR, vP_vv,vP_vZ, ncol=3)

data_vr_vv<-vr_vv$data
tmp<-ggarrange(A, B, nrow=2)
ggsave("ch13_dxy_end_sites.pdf", tmp, dpi=300, width=9 , height=7, units="in")

plotfst<-function(myfile){
  dxy <-read.delim(myfile, sep = '\t', header = F)
  colnames(dxy) <- c("chrom", "winstart", "winend", "dxy", "nsites", "na", "genes")
  dxy$length<-dxy$winend - dxy$winstart
  dxy<-dxy[dxy$nsites>1000,] 
  dxy$dxy_relative <- dxy$dxy / dxy$nsites
  scaffold_order<-read.delim("scaffolds_ordered_DOG.txt", header = FALSE)
  scaffold_order$chrom<-scaffold_order$V2
  scaffold_order$ranking<-c(1:length(scaffold_order$V2))
  test<-merge(scaffold_order, dxy, by="chrom", all.y=TRUE)
  test<-na.omit(test)
  test2 <-test[
    with(test, order(test$V2, test$ranking, test$winstart)),
  ]

  test2<-subset(test2, V1!=15) ### to remove CRADD
  test2<-subset(test2, V1!=39)
  test2<-subset(test2, V1!=13)
  threshold1<-quantile(test2$dxy_relative, 0.95, na.rm = TRUE)
  threshold2<-quantile(test2$dxy_relative, 0.05, na.rm = TRUE)
  #test2<-subset(test2, V1==13)
  test2<-subset(test2, V1 == 1 & test2$winstart >  107049560  & test2$winstart < 107284526 ) #SLC6A16
 # test2<-subset(test2, V1 == 16 & test2$winstart >   6974288   & test2$winstart < 7375902) #MGAM
  #test2<-subset(test2, V1 == 21 & test2$winstart >  40548003   & test2$winstart <  40890289) #HPS5
  #test2<-subset(test2, target=="tmp2")
  test2$Window<-c(1:length(test2$V2))
  #test2<-subset(test2, winstart>=59784933)
  #test2<-subset(test2, Window>=6145)
  cbPalette<- c( "#d5d5d4")
  n_of_chr<-length(table(test2$V2))
  ggplot(test2, aes(x=winstart, y=dxy_relative, color=V2)) + ylab("Dxy") + xlab("") +
    ggtitle("") +
    geom_point(size=.5)+
    theme_classic(base_size=20)  + #coord_cartesian(ylim = c(0,2))+
    theme(legend.position="None", axis.title=element_text(size=11)) +
    theme(plot.title = element_text(size = 11, face = "bold")) +
   # theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
    scale_color_manual(values=rep(cbPalette,n_of_chr))  +
    geom_hline(yintercept =threshold1,linetype="dashed") + 
    geom_hline(yintercept =threshold2,linetype="dashed")
}

vr_vvNA<-plotfst("dxy_vr_vvNA_25kb5kb_sites_genes") + ylab("dxy(Rueppell's,Red fox NA)") + xlab("")  #+
#  geom_label_repel(aes(label=ifelse(genes == "HPS5",as.character(genes),'')),box.padding= 0.35, point.padding = 0.5,size=3) + scale_x_continuous(breaks = seq(0, 3000, by = 100))
vr_vvEU<-plotfst("dxy_vr_vvEU_25kb5kb_sites_genes") + ylab("dxy(Rueppell's,Red fox EU)") #+
  #geom_label_repel(aes(label=ifelse(genes == "SLC6A16",as.character(genes),'')),box.paddi#ng= 0.35, point.padding = 0.5,size=3) 
# + scale_x_continuous(breaks = seq(0, 3000, by = 100))
vvNA_vvEU<-plotfst("dxy_vvNA_vvEU_25kb5kb_sites_genes") + ylab("dxy(Red fox NA,Red fox EU)") #+
  #geom_label_repel(aes(label=ifelse(genes == "MGAM",as.character(genes),'')),box.padding= 0.35, point.padding = 0.5,size=3) 
 #+ scale_x_continuous(breaks = seq(0, 3000, by = 100)) # + scale_x_continuous(breaks = seq(0, 3000, by = 100))
SLC6A16<-ggarrange(vr_vvEU,vr_vvNA, vvNA_vvEU, nrow = 3)
ggsave("SLC6A16.tiff", SLC6A16, "tiff", dpi=300, width=3 , height=7, units="in")
ggsave("SLC6A16.pdf", SLC6A16, "pdf", dpi=300, width=3 , height=7, units="in")
ggsave("SLC6A16_dingbats.pdf", SLC6A16, "pdf", dpi=300, width=3 , height=7, units="in", useDingbats=TRUE)
ggsave("SLC6A16.eps", SLC6A16, "eps", dpi=300, width=3 , height=7, units="in")
ggsave("SLC6A16.svg", SLC6A16, "svg", dpi=300, width=3 , height=7, units="in")


vr_vv<-plotfst("dxy_vr_vv_50kb10kb_sites_genes") + ylab("dxy(Rueppell's,Red fox)") # + scale_x_continuous(breaks = seq(0, 3000, by = 100))
vv_vZ<-plotfst("dxy_vv_vz_50kb10kb_sites_genes") + ylab("dxy(Fennec,Red fox)")
vr_vZ<-plotfst("dxy_vr_vz_50kb10kb_sites_genes") +  ylab("dxy(Rueppell's,Fennec)") #+ scale_x_continuous(breaks = seq(0, 3000, by = 100))
vP_vr<-plotfst("dxy_vr_vp_50kb10kb_sites_genes") + ylab("dxy(Pale,Rueppell'sfox)")
vP_vv<-plotfst("dxy_vv_vp_50kb10kb_sites_genes") + ylab("dxy(Pale,Red fox)")
vP_vZ<-plotfst("dxy_vz_vp_50kb10kb_sites_genes") + ylab("dxy(Pale,Fennec fox)")

#vr_vv<-plotfst("dxy_vr_vv_5000_1250_snps_genes") + ylab("dxy(Rueppell's,Red fox)") # + scale_x_continuous(breaks = seq(0, 3000, by = 100))
#vv_vZ<-plotfst("dxy_vz_vv_5000_1250_snps_genes") + ylab("dxy(Fennec,Red fox)")
#vr_vZ<-plotfst("dxy_vz_vr_5000_1250_snps_genes") +  ylab("dxy(Rueppell's,Fennec)") #+ scale_x_continuous(breaks = seq(0, 3000, by = 100))
#vP_vr<-plotfst("dxy_vp_vr_5000_1250_snps_genes") + ylab("dxy(Pale,Rueppell'sfox)")
#vP_vv<-plotfst("dxy_vp_vv_5000_1250_snps_genes") + ylab("dxy(Pale,Red fox)")
#vP_vZ<-plotfst("dxy_vp_vz_5000_1250_snps_genes") + ylab("dxy(Pale,Fennec fox)")
data_vr_vv<-vr_vv$data
data_vr_vZ<-vr_vZ$data
write.csv(data_vr_vZ, "~/Documents/work/PhD/Tasks/Task4_POPGEN_final/Dxyz/Dxy_chr13_end_vRvZ.csv", row.names = FALSE)

ch13_5000snps_1250<-ggarrange(vP_vZ, vP_vr, vr_vZ,  nrow = 3)
ch13_5000snps_1250v2<-ggarrange(vP_vv,vv_vZ,vr_vv, nrow = 3)
final<-ggarrange(ch13_5000snps_1250v2,ch13_5000snps_1250, ncol=2)
ggsave("chr13_all_dxy_final.tiff", final, "tiff", dpi=600, width=11 , height=7, units="in")


