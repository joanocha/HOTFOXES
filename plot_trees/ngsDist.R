setwd("~/Documents/work/PhD/Tasks/Task4_POPGEN_final/ngsDist/")

library(phangorn)
library(ggpubr)
library(Hmisc)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(plotly)
library(psych)
library(ape)
library(nlme)
library(GGally)
library(geiger)
library(ggtree)
library(phytools)
library(phylobase)
library(plyr)  ## or dplyr (transform -> mutate)
library(reshape2)

# WHOLE GENOME - JUST INGROUPS
key <- read.delim("allVulpes.txt", na.strings=c("", "NA"), header=FALSE, col.names=c('tip.label', 'pops', 'continent', 'species', 'old'), row.names=1)
tree_justVulpes<-read.tree("CRADD_bams.fasta.contree.tree")
tree_justVulpes <- midpoint(tree_justVulpes, node.labels = "support")
prun_list<-name.check(tree_justVulpes, key, data.names= key$tip.label)
key <-key[tree_justVulpes$tip.label,] #ensure key and tree are order the same way

tree_justVulpes$tip.label<-key$pops
tree_justVulpes<-drop.tip(tree_justVulpes, c("HW", "Grey_wolf4", "Culpeo_Fox", "Dingo", "Basenji_dog", "Grey_wolf3", "Island_fox",
                                             "TW1", "TW2", "TW3", "Grey_wolf5","Grey_wolf6",  "Grey_wolf2", "GJ2", "C_Mexico", "GS_dog", "BBJ",
                                             "EW", "GJ1", "AGW", "Dhole", "AWD", "C_Alaska", "NA26", "himalayan wolf", "Polar bear", "Gray fox", "Island fox", 
                                             "Red fox (Russian farm)"
                                            ))
tree_justVulpes<-drop.tip(tree_justVulpes, c("Dog (Tibetan Mastiff)", "Fenenc fox", "Pale fox", "NA"))
tree_justVulpes<-drop.tip(tree_justVulpes, c("dhole", "black_backed_jackal", "ethiopian_wolf", "dingo"))


plotTree(tree_justVulpes,ftype="i",fsize=1,lwd=3)
#nodelabels(raxml.tree$node.label, adj = c(1, 0), frame = "none")
wholegenome<-plot(tree_justVulpes, show.tip.label=FALSE, edge.color = TRUE, align.tip.label = FALSE)
ggtree(tree_justVulpes, layout="daylight", size=1) + geom_tiplab()


p<-ggtree(tree_justVulpes, color="dark grey",size=0.7) + geom_tiplab()
ggsave("CRADDonly.eps",p,height=5,width=3, dpi = 600)
