setwd("~/Documents/backup_10April2019/PhD/Tasks/Task2_Phenotyping/Analysis/phyenotypic_Analyses/R/myPhysiology/VULPES")
#useful packages
library(ggpubr)
library(Hmisc)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(psych)
library(car)
library(ggpubr)
## DATASETS AND SUBSETS
Foxbiochem<-read.delim('Fox_morpho_physiology_with_juveniles.txt')
Foxbiochem2<-read.delim('Fox_morpho_physiology.txt') # juveniles excluded
vulpes <-Foxbiochem2[24:34,]
zerda<-Foxbiochem2[1:10,]
ruep<-Foxbiochem2[11:23,]
cbPalette_vv<- c( "#d5d5d4","#e69f00")
cbPalette_vvNA<- c( "#ff7f00","#d5d5d4")
cbPalette_vvEU<- c( "#33a02c", "#d5d5d4")
cbPalette_vR<- c("#e31a1c", '#d5d5d4')
cbPalette_vZ<- c("#6C1717", '#d5d5d4')

#DESCRIPTIVE STATISTICS
#Bodymass and Biometrics
weight <-ggplot(Foxbiochem2, aes(x=Species1, y=body.mass.Kg)) + geom_boxplot(colour= c("orange", "#00AFBB", "#FC4E07")) + xlab("") + ylab("Weight (Kg)") 
tail<-ggplot(Foxbiochem2, aes(x=Species1, y=tail.length.mm)) + geom_boxplot(colour= c("orange", "#00AFBB", "#FC4E07")) + xlab("") + ylab("Tail lenght (mm)") 
body <-ggplot(Foxbiochem2, aes(x=Species1, y=body.length.mm)) +  geom_boxplot(colour= c("orange", "#00AFBB", "#FC4E07")) + xlab("") + ylab("Body length (mm)") 
hindfoot <-ggplot(Foxbiochem2, aes(x=Species1, y=hindfoot.mm)) + geom_boxplot(colour= c("orange", "#00AFBB", "#FC4E07")) + xlab("") + ylab("Hindfoot length (mm)")
ear <-ggplot(Foxbiochem2, aes(x=Species1, y= ear.length.mm)) + geom_boxplot(colour= c("orange", "#00AFBB", "#FC4E07")) + xlab("") + ylab("Ear length (mm)") 
bmi <-ggplot(Foxbiochem2, aes(x=Species1, y=BMI)) + geom_boxplot(colour= c("orange", "#00AFBB", "#FC4E07")) + xlab("") + ylab("Body Mass Index (Kg/m*m)") 

# CORRELATIONS: Between x variables 
bodymass.length<-ggscatter(Foxbiochem2, x = "body.length.mm", y = "body.mass.Kg", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Body length (mm)", ylab = "Weight (Kg)")

bodytail <-ggscatter(Foxbiochem2, x = "body.length.mm", y = "tail.length.mm", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Body length (mm)", ylab = "Tail lenght (mm)")

bodyhindfoot<-ggscatter(Foxbiochem2, x = "body.length.mm", y = "hindfoot.mm", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Body length (mm)", ylab = "Hindfoot length (mm)")

bodyear<-ggscatter(Foxbiochem2, x = "body.length.mm", y = "ear.length.mm", 
                        add = "reg.line", conf.int = TRUE, 
                        cor.coef = TRUE, cor.method = "spearman",
                        xlab = "Body length (mm)", ylab = "Ear length (mm)")

#ggarrange(weight,body,hindfoot,tail, labels = c("A", "B", "C", "D"), ncol=2, nrow=2)
FigureS_morphology<-ggarrange(bodymass.length, bodyhindfoot, bodytail, labels = c("A", "B", "C"), ncol=3, nrow=1)
ggsave("FigureS_morphology.pdf", FigureS_morphology,"pdf", dpi=400, width=10,height=5, units="in")

  

#Hormones
describeBy(Foxbiochem2$copeptin.pmol.L, Foxbiochem2$Species1)
describeBy(Foxbiochem2$aldosterone.pg.mL, Foxbiochem2$Species1)
describeBy(Foxbiochem2$T4.ug.dL, Foxbiochem2$Species1)

T4 <-ggplot(Foxbiochem2, aes(x=Species1, y=T4.ug.dL)) + geom_boxplot(colour= c("#ff7f00", "#e31a1c", "#6C1717"))  + scale_x_discrete(limits=c("Vulpes vulpes", "Vulpes rueppellii", "Vulpes zerda")) + xlab("") + ylab("Total T4 levels (ug/dL)") 
T4stat <-ggplot(Foxbiochem2, aes(x=Status, y=T4.ug.dL)) + geom_boxplot(colour= c("black", "black")) + scale_x_discrete(limits=c("captive", "wild")) + xlab("") + ylab("Total T4 levels (ug/dL)") 

copeptin <-ggplot(Foxbiochem2, aes(x=Species1, y=copeptin.pmol.L)) + geom_boxplot(colour= c("#ff7f00", "#e31a1c", "#6C1717")) + scale_x_discrete(limits=c("Vulpes vulpes", "Vulpes rueppellii", "Vulpes zerda")) + xlab("") + ylab("Copeptin (pmol/L)") 
Figure5_D <-ggplot(Foxbiochem2, aes(x=Species1, y=copeptin.pmol.L)) + geom_boxplot(colour= c("#ff7f00", "#e31a1c", "#6C1717")) + scale_x_discrete(limits=c("Vulpes vulpes", "Vulpes rueppellii", "Vulpes zerda")) + xlab("") + ylab("Copeptin (pmol/L)") + theme_classic()
ggsave("Figure5_D.pdf", Figure5_D,"pdf", dpi=400, width=4,height=3, units="in")

copstat <-ggplot(Foxbiochem2, aes(x=Status, y=copeptin.pmol.L)) + geom_boxplot(colour= c("black", "black")) + scale_x_discrete(limits=c("captive", "wild")) + xlab("") + ylab("Copeptin (pmol/L)") 
aldosterone <-ggplot(Foxbiochem2, aes(x=Species1, y=aldosterone.pg.mL)) +geom_boxplot(colour= c("#ff7f00", "#e31a1c", "#6C1717")) + scale_x_discrete(limits=c("Vulpes vulpes", "Vulpes rueppellii", "Vulpes zerda")) + xlab("") + ylab("Aldosterone (pg/mL)") 
aldstat <-ggplot(Foxbiochem2, aes(x=Status, y=aldosterone.pg.mL)) + geom_boxplot(colour= c("black", "black")) + scale_x_discrete(limits=c("captive", "wild")) + xlab("") + ylab("Aldosterone (pg/mL)") 

ggarrange(copeptin, aldosterone, T4, labels = c("A","B", "C" ),  ncol = 3)
#ggarrange(copstat, aldstat, T4stat, labels = c("A","B", "C" ),  ncol = 2, nrow = 2)

logcopeptin <-ggplot(Foxbiochem2, aes(x=Species1, y=log(copeptin.pmol.L))) + geom_boxplot(colour= c("#00AFBB", "orange", "#FC4E07")) + scale_x_discrete(limits=c("Vulpes vulpes", "Vulpes rueppellii", "Vulpes zerda")) + xlab("") + ylab("Ln Copeptin (pmol/L)")
logaldosterone <-ggplot(Foxbiochem2, aes(x=Species1, y=log(aldosterone.pg.mL))) +geom_boxplot(colour= c("#00AFBB", "orange", "#FC4E07")) + scale_x_discrete(limits=c("Vulpes vulpes", "Vulpes rueppellii", "Vulpes zerda")) + xlab("") + ylab("Ln Aldosterone (pg/mL)") 


#Metabolites, Proteins, Fat

## Summary stats
#all
describeBy(Foxbiochem$urea.mg.dL, Foxbiochem2$Species1)
describeBy(Foxbiochem$uric.acid.mg.dl, Foxbiochem2$Species1)
describeBy(Foxbiochem$creatinin.mg.dL, Foxbiochem2$Species1)
describeBy(Foxbiochem$albumin.g.dL, Foxbiochem2$Species1)
describeBy(Foxbiochem$cholesterol.mg.dL, Foxbiochem2$Species1)

#no juveniles
describeBy(Foxbiochem2$urea.mg.dL, Foxbiochem2$Species1)
describeBy(Foxbiochem2$uric.acid.mg.dl, Foxbiochem2$Species1)
describeBy(Foxbiochem2$creatinin.mg.dL, Foxbiochem2$Species1)
describeBy(Foxbiochem2$albumin.g.dL, Foxbiochem2$Species1)
describeBy(Foxbiochem2$cholesterol.mg.dL, Foxbiochem2$Species1)
urea <- ggplot(Foxbiochem2, aes(x=Species1, y=urea.mg.dL)) + geom_boxplot(colour= c("#ff7f00", "#e31a1c", "#6C1717")) + scale_x_discrete(limits=c("Vulpes vulpes", "Vulpes rueppellii", "Vulpes zerda")) + xlab("") + ylab("Urea (mg/dL)") 
uacid <- ggplot(Foxbiochem2, aes(x=Species1, y=uric.acid.mg.dl)) + geom_boxplot(colour= c("#ff7f00", "#e31a1c", "#6C1717")) + scale_x_discrete(limits=c("Vulpes vulpes", "Vulpes rueppellii", "Vulpes zerda"))+ xlab("") + ylab("Uric acid (mg/dL)") 
crea <- ggplot(Foxbiochem2, aes(x=Species1, y=creatinin.mg.dL)) + geom_boxplot(colour= c("#ff7f00", "#e31a1c", "#6C1717")) + scale_x_discrete(limits=c("Vulpes vulpes", "Vulpes rueppellii", "Vulpes zerda")) + xlab("") + ylab("Creatinine (mg/dL)") 
alb <- ggplot(Foxbiochem2, aes(x=Species1, y=albumin.g.dL)) + geom_boxplot(colour= c("#ff7f00", "#e31a1c", "#6C1717")) + scale_x_discrete(limits=c("Vulpes vulpes", "Vulpes rueppellii", "Vulpes zerda")) + xlab("") + ylab("Albumin (g/dL)") 
chol <-ggplot(Foxbiochem2, aes(x=Species1, y=cholesterol.mg.dL)) + geom_boxplot(colour= c("#ff7f00", "#e31a1c", "#6C1717")) + scale_x_discrete(limits=c("Vulpes vulpes", "Vulpes rueppellii", "Vulpes zerda")) + xlab("") + ylab("Cholesterol (mg/dL)")

#electrolytes and plasma osmolality - Just for VV and VR
Na<-ggplot(Foxbiochem, aes(x=Species1, y=sodium.mEq.L)) + geom_boxplot(colour= c("#ff7f00", "#e31a1c")) + scale_x_discrete(limits=c("Vulpes vulpes", "Vulpes rueppellii")) + xlab("") + ylab("Sodium (mEq/L)") 
Po<-ggplot(Foxbiochem, aes(x=Species1, y=potassium.mEq.L)) + geom_boxplot(colour= c("#ff7f00", "#e31a1c")) + scale_x_discrete(limits=c("Vulpes vulpes", "Vulpes rueppellii")) + xlab("") + ylab("Potassium (mEq/L)") 
Ch<-ggplot(Foxbiochem, aes(x=Species1, y=chloride.mEq.L)) + geom_boxplot(colour= c("#ff7f00", "#e31a1c")) + scale_x_discrete(limits=c("Vulpes vulpes", "Vulpes rueppellii")) + xlab("") + ylab("Chloride (mEq/L)") 
Osm<-ggplot(Foxbiochem, aes(x=Species1, y=Osmolality.mOsm.Kg)) + geom_boxplot(colour= c("#ff7f00", "#e31a1c")) + scale_x_discrete(limits=c("Vulpes vulpes", "Vulpes rueppellii")) + xlab("") + ylab("Osmolality (mOsm/Kg)") 

describeBy(Foxbiochem$sodium.mEq.L, Foxbiochem$Species1)
describeBy(Foxbiochem$potassium.mEq.L, Foxbiochem$Species1)
describeBy(Foxbiochem$chloride.mEq.L, Foxbiochem$Species1)
describeBy(Foxbiochem$Osmolality.mOsm.Kg, Foxbiochem$Species1)

param<-ggarrange(copeptin, aldosterone, T4, urea, uacid, crea, alb, chol,Osm, Na, Po, Ch, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"), ncol = 3, nrow = 4)
ggsave("param.pdf", param,"pdf", dpi=300, width=10,height=10, units="in")
ggsave("param.eps", dpi=300, width=11,height=11, units="in")
# CORRELATIONS: Between Hormones and Body mass
bmt4<-ggscatter(Foxbiochem2, x = "body.mass.Kg", y = "T4.ug.dL", title = "All species",add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
          xlab = "Body mass (Kg)", ylab = "T4 (ug/dL)")
vbmt4<-ggscatter(vulpes, x = "body.mass.Kg", y = "T4.ug.dL", title = "Red fox", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                 xlab = "Body mass (Kg)", ylab = "T4 (ug/dL)L")
rbmt4<-ggscatter(ruep, x = "body.mass.Kg", y = "T4.ug.dL", title = "Rueppell's fox",add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                  xlab = "Body mass (Kg)", ylab = "T4 (ug/dL)")
zbmt4<-ggscatter(zerda, x = "body.mass.Kg", y = "T4.ug.dL", title = "Fennec fox",add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                  xlab = "Body mass (Kg)", ylab = "T4 (ug/dL)")
bmcop<-ggscatter(Foxbiochem2, x = "body.mass.Kg", y = "copeptin.pmol.L", title = "All species",add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
          xlab = "Body mass (Kg)", ylab = "Copeptin (pmol/L)")
vbmcop<-ggscatter(vulpes, x = "body.mass.Kg", y = "copeptin.pmol.L", title = "Red fox",add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
          xlab = "Body mass (Kg)", ylab = "Copeptin (pmol/L)")
rbmcop<-ggscatter(ruep, x = "body.mass.Kg", y = "copeptin.pmol.L", title = "Rueppell's fox",add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                   xlab = "Body mass (Kg)", ylab = "Copeptin (pmol/L)")
zbmcop<-ggscatter(zerda, x = "body.mass.Kg", y = "copeptin.pmol.L", title = "Fennec fox",add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                   xlab = "Body mass (Kg)", ylab = "Copeptin (pmol/L)")
bmald<-ggscatter(Foxbiochem2, x = "body.mass.Kg", y = "aldosterone.pg.mL", title = "All species", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                    xlab = "Body mass (Kg)", ylab = "Aldosterone (pg/mL)")
vbmald<-ggscatter(vulpes, x = "body.mass.Kg", y = "aldosterone.pg.mL", title = "Red fox",add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                    xlab = "Body mass (Kg)", ylab = "Aldosterone (pg/mL)")
rbmald<-ggscatter(ruep, x = "body.mass.Kg", y = "aldosterone.pg.mL", title = "Rueppell's fox",add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                    xlab = "Body mass (Kg)", ylab = "Aldosterone (pg/mL)")
zbmald<-ggscatter(zerda, x = "body.mass.Kg", y = "aldosterone.pg.mL", title = "Fennec fox",add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                    xlab = "Body mass (Kg)", ylab = "Aldosterone (pg/mL)")

ggarrange(bmcop, vbmcop, rbmcop, zbmcop,bmald, vbmald, rbmald, zbmald, bmt4, vbmt4, rbmt4, zbmt4,
          labels = c("A", "d0", "d1", "d2", 
                     "B", "d0", "d1", "d2",
                     "C", "d0", "d1", "d1"))

# CORRELATIONS: Between Y variables
##### Hormones among each other
t4cop <-cor.test(Foxbiochem2$T4.ug.dL, Foxbiochem2$copeptin.pmol.L,  method = "spearman")
t4ald<-cor.test(Foxbiochem2$T4.ug.dL, Foxbiochem2$aldosterone.pg.mL,  method = "spearman")
copald<-cor.test(Foxbiochem2$copeptin.pmol.L, Foxbiochem2$aldosterone.pg.mL,  method = "spearman")

t4cop<-ggscatter(Foxbiochem2, x = "T4.ug.dL", y = "copeptin.pmol.L", title = "All species",add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
          xlab = "T4 (ug/dL)", ylab = "Copeptin (pmol/L)")
vt4cop<-ggscatter(vulpes, x = "T4.ug.dL", y = "copeptin.pmol.L", title = "Red fox",add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
          xlab = "T4 (ug/dL)", ylab = "Copeptin pmol/L")
rt4cop<-ggscatter(ruep, x = "T4.ug.dL", y = "copeptin.pmol.L", title = "Rueppell's fox",add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
          xlab = "T4 (ug/dL)", ylab = "Copeptin (pmol/L)")
zt4cop<-ggscatter(zerda, x = "T4.ug.dL", y = "copeptin.pmol.L", title = "Fennec fox",add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
          xlab = "T4 (ug/dL)", ylab = "Copeptin (pmol/L)")
t4ald<-ggscatter(Foxbiochem2, x = "T4.ug.dL", y = "aldosterone.pg.mL", title = "All species", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                 xlab = "T4 (ug/dL)", ylab = "Aldosterone (pg/mL)")
vt4ald<-ggscatter(vulpes, x = "T4.ug.dL", y = "aldosterone.pg.mL", title= "Red fox", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
          xlab = "T4 (ug/dL)", ylab = "Aldosterone (pg/mL)")
rt4ald<-ggscatter(ruep, x = "T4.ug.dL", y = "aldosterone.pg.mL", title= "Rueppell's fox", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
          xlab = "T4 (ug/dL)", ylab = "Aldosterone (pg/mL)")
zt4ald<-ggscatter(zerda, x = "T4.ug.dL", y = "aldosterone.pg.mL", title= "Fennec fox",add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
          xlab = "T4 (ug/dL)", ylab = "Aldosterone (pg/mL)")
copald<-ggscatter(Foxbiochem2, x = "copeptin.pmol.L", y = "aldosterone.pg.mL", title = "All species", add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Copeptin (pmol/L)", ylab = "Aldosterone (pg/mL)")

vcopald<-ggscatter(vulpes, x = "copeptin.pmol.L", y = "aldosterone.pg.mL", title =  "Red fox", add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman",
                  xlab = "Copeptin (pmol/L)", ylab = "Aldosterone (pg/mL)")
rcopald<-ggscatter(ruep, x = "copeptin.pmol.L", y = "aldosterone.pg.mL", title =  "Rueppell's fox", add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "spearman",
                   xlab = "Copeptin (pmol/L)", ylab = "Aldosterone (pg/mL)")
zcopald<-ggscatter(zerda, x = "copeptin.pmol.L", y = "aldosterone.pg.mL", title =  "Fennec fox", add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "spearman",
                   xlab = "Copeptin (pmol/L)", ylab = "Aldosterone (pg/mL)")

ggarrange(t4cop, vt4cop, rt4cop, zt4cop,t4ald, vt4ald, rt4ald, zt4ald, copald, vcopald, rcopald, zcopald, 
          labels = c("A", "d0", "d1", "d2",
                     "B", "d0", "d1", "d2",
                     "C", "d0", "d1", "d2"))

##### Hormones and Metabolites and Electrolytes
t4urea<-ggscatter(Foxbiochem2, x = "T4.ug.dL", y = "urea.mg.dL", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                  xlab = "T4 (ug/dL)", ylab = "Urea(mg/dL)")
t4uric<-ggscatter(Foxbiochem2, x = "T4.ug.dL", y = "uric.acid.mg.dl", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                  xlab = "T4 (ug/dL)", ylab = "Uric acid (mg/dL)")

t4crea<-ggscatter(Foxbiochem2, x = "T4.ug.dL", y = "creatinin.mg.dL", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
          xlab = "T4 (ug/dL)", ylab = "Creatinine (mg/dL)")

t4alb<-ggscatter(Foxbiochem2, x = "T4.ug.dL", y = "albumin.g.dL", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                 xlab = "T4 (ug/dL)", ylab = "Albumin (g/dL)")

t4chol<-ggscatter(Foxbiochem2, x = "T4.ug.dL", y = "cholesterol.mg.dL", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                  xlab = "T4 (ug/dL)", ylab = "Cholesterol (mg/dL)")

t4Na<-ggscatter(Foxbiochem2, x = "T4.ug.dL", y = "sodium.mEq.L", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                  xlab = "T4 (ug/dL)", ylab = "Sodium (mEq/L)")
t4K<-ggscatter(Foxbiochem2, x = "T4.ug.dL", y = "potassium.mEq.L", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                xlab = "T4 (ug/dL)", ylab = "Potassium (mEq/L)")

t4Cl<-ggscatter(Foxbiochem2, x = "T4.ug.dL", y = "chloride.mEq.L", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                xlab = "T4 (ug/dL)", ylab = "Chloride (mEq/L)")
t4osm<-ggscatter(Foxbiochem2, x = "T4.ug.dL", y = "Osmolality.mOsm.Kg", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                 xlab = "T4 (ug/dL)", ylab = "Osmolality.mOsm.Kg")
copurea<-ggscatter(Foxbiochem2, x = "copeptin.pmol.L", y = "urea.mg.dL", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                   xlab = "Copeptin (pmol/L)", ylab = "Urea(mg/dL)")
copuric<-ggscatter(Foxbiochem2, x = "copeptin.pmol.L", y = "uric.acid.mg.dl", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                  xlab = "Copeptin (pmol/L)", ylab = "Uric acid (mg/dL)")
copcrea<-ggscatter(Foxbiochem2, x = "copeptin.pmol.L", y = "creatinin.mg.dL", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                   xlab = "Copeptin (pmol/L)", ylab = "Creatinine (mg/dL)")
copalb<-ggscatter(Foxbiochem2, x = "copeptin.pmol.L", y = "albumin.g.dL", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                 xlab = "Copeptin (pmol/L)", ylab = "Albumin (g/dL)")
copchol<-ggscatter(Foxbiochem2, x = "copeptin.pmol.L", y = "cholesterol.mg.dL", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                  xlab = "Copeptin (pmol/L)", ylab = "Cholesterol (mg/dL)")

copNa<-ggscatter(Foxbiochem2, x = "copeptin.pmol.L", y = "sodium.mEq.L", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                xlab = "Copeptin (pmol/L)", ylab = "Sodium (mEq/L)")
copK<-ggscatter(Foxbiochem2, x = "copeptin.pmol.L", y = "potassium.mEq.L", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
               xlab = "Copeptin (pmol/L)", ylab = "Potassium (mEq/L)")

copCl<-ggscatter(Foxbiochem2, x = "copeptin.pmol.L", y = "chloride.mEq.L", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                xlab = "Copeptin (pmol/L)", ylab = "Chloride (mEq/L)")
coposm<-ggscatter(Foxbiochem2, x = "copeptin.pmol.L", y = "Osmolality.mOsm.Kg", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                 xlab = "Copeptin (pmol/L)", ylab = "Osmolality.mOsm.Kg")
aldurea<-ggscatter(Foxbiochem2, x = "aldosterone.pg.mL", y = "urea.mg.dL", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                   xlab = "Aldosterone (pg/mL)", ylab = "Urea(mg/dL)")
alduric<-ggscatter(Foxbiochem2, x = "aldosterone.pg.mL", y = "uric.acid.mg.dl", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                   xlab = "Aldosterone (pg/mL)", ylab = "Uric acid (mg/dL)")
aldcrea<-ggscatter(Foxbiochem2, x = "aldosterone.pg.mL", y = "creatinin.mg.dL", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                   xlab = "Aldosterone (pg/mL)", ylab = "Creatinine (mg/dL)")
aldalb<-ggscatter(Foxbiochem2, x = "aldosterone.pg.mL", y = "albumin.g.dL", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                  xlab = "aldosterone.pg.mL", ylab = "Albumin (g/dL)")
aldchol<-ggscatter(Foxbiochem2, x = "aldosterone.pg.mL", y = "cholesterol.mg.dL", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                   xlab = "Aldosterone (pg/mL)", ylab = "Cholesterol (mg/dL)")

aldNa<-ggscatter(Foxbiochem2, x = "aldosterone.pg.mL", y = "sodium.mEq.L", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                 xlab = "Aldosterone (pg/mL)", ylab = "Sodium (mEq/L)")
aldK<-ggscatter(Foxbiochem2, x = "aldosterone.pg.mL", y = "potassium.mEq.L", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                xlab = "Aldosterone (pg/mL)", ylab = "Potassium (mEq/L)")

aldCl<-ggscatter(Foxbiochem2, x = "aldosterone.pg.mL", y = "chloride.mEq.L", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                 xlab = "Aldosterone (pg/mL)", ylab = "Chloride (mEq/L)")
aldosm<-ggscatter(Foxbiochem2, x = "aldosterone.pg.mL", y = "Osmolality.mOsm.Kg", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                  xlab = "Aldosterone (pg/mL)", ylab = "Osmolality.mOsm.Kg")


# significant correlations
t4urea<-ggscatter(Foxbiochem2, x = "T4.ug.dL", y = "urea.mg.dL", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                  xlab = "T4 (ug/dL)", ylab = "Urea(mg/dL)")

aldcrea<-ggscatter(Foxbiochem2, x = "aldosterone.pg.mL", y = "creatinin.mg.dL", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman",
                   xlab = "Aldosterone (pg/mL)", ylab = "Creatinine (mg/dL)")

ggarrange(t4urea, aldcrea)
### STATISTICS - HORMONES 

### COPEPTIN
pairwise.t.test(Foxbiochem2$copeptin.pmol.L, Foxbiochem2$Species1,p.adjust.method = "bonferroni") # only significant between vv e vr
t.test(Foxbiochem2$copeptin.pmol.L~Foxbiochem2$Status) #non signficant for captive vs wild
t.test(Foxbiochem2$copeptin.pmol.L~Foxbiochem2$Sex) #non signficant for M vs F
t.test(vulpes$copeptin.pmol.L~vulpes$Status) #non significant within V. vulpes for captive and wild
t.test(vulpes$copeptin.pmol.L~vulpes$Country)

# LINEAR MODEL - COPEPTIN

fitcop <-lm(copeptin.pmol.L~ Species2 +  body.mass.Kg +  BMI  + Sex + Status, data = Foxbiochem2)
summary.lm(fitcop)
summary.aov(fitcop)
# Tested  for linear model assumptions
shapiro.test(as.numeric(residuals(object = fitcop))) #  normal
leveneTest(copeptin.pmol.L~Species1, data =Foxbiochem2) #not homogeneous
bartlett.test(copeptin.pmol.L~Species1, data =Foxbiochem2) # not homogenous

# log transformation
fitlogcop <-lm(log(copeptin.pmol.L)~ Species2 + body.mass.Kg + BMI + Sex + Status, data = Foxbiochem2)

plot(fitlogcop)
shapiro.test(as.numeric(residuals(object = fitlogcop))) # normal
leveneTest(log(copeptin.pmol.L)~Species2, data =Foxbiochem2) #homogeneous
bartlett.test(log(copeptin.pmol.L)~Species2, data =Foxbiochem2) #homogenous

summary.lm(fitlogcop)
summary.aov(fitlogcop)
TukeyHSD(aov(fitlogcop), which = 'Species2')

fitlogcop2 <-lm(log(copeptin.pmol.L)~ Species2 + T4.ug.dL + aldosterone.pg.mL + body.mass.Kg + BMI + Sex + Status, data = Foxbiochem2)
summary.lm(fitlogco2p)
summary.aov(fitlogcop2)
shapiro.test(as.numeric(residuals(object = fitlogcop2))) # normal
### ALDOSTERONE
pairwise.t.test(Foxbiochem2$aldosterone.pg.mL, Foxbiochem2$Species1,p.adjust.method = "bonferroni") 

t.test(Foxbiochem2$aldosterone.pg.mL~Foxbiochem$Status) #non signficant for captive vs wild
t.test(Foxbiochem2$aldosterone.pg.mL~Foxbiochem$Sex) #non signficant for M vs F
t.test(vulpes$aldosterone.pg.mL~vulpes$Status) #non significant within V. vulpes for captive and wild
t.test(vulpes$aldosterone.pg.mL~vulpes$Country) #nonn significant
# LINEAR MODEL - ALDOSTERONE
fitald <- lm(aldosterone.pg.mL~ Species2 + body.mass.Kg + BMI + Sex + Status, data = Foxbiochem2)
shapiro.test(as.numeric(residuals(object = fitald))) # not normal
leveneTest(aldosterone.pg.mL.L~Species2, data =Foxbiochem2) #not homogeneous
bartlett.test(aldosterone.pg.mL~Species2, data =Foxbiochem2) # not homogenous

fitlogald <- lm(log(aldosterone.pg.mL)~ Species2 + body.mass.Kg + BMI + Sex + Status, data = Foxbiochem2)
summary.lm(fitlogald)
summary.aov(fitlogald)

shapiro.test(as.numeric(residuals(object = fitlogald))) # normal
leveneTest(log(aldosterone.pg.mL)~Species2, data =Foxbiochem2) #homogeneous
bartlett.test(log(aldosterone.pg.mL)~Species2, data =Foxbiochem2) #homogenous

fitlogald2 <- lm(log(aldosterone.pg.mL)~ Species2 + T4.ug.dL + copeptin.pmol.L + body.mass.Kg + BMI + Sex + Status, data = Foxbiochem2)
summary.lm(fitlogald2)
summary.aov(fitlogald2)
shapiro.test(as.numeric(residuals(object = fitlogald2))) # normal


TukeyHSD(aov(fitlogald), which = 'Species2')

### T4
pairwise.t.test(Foxbiochem2$T4.ug.dL, Foxbiochem2$Species1,p.adjust.method = "bonferroni") # SIGNIFICANT: vv vs vz, vv vs vr 

t.test(Foxbiochem2$T4.ug.dL~Foxbiochem2$Status) #non signficant for captive vs wild
t.test(Foxbiochem2$T4.ug.dL~Foxbiochem2$Sex) #non signficant for M vs F
t.test(vulpes$T4.ug.dL~vulpes$Status) #non significant within V. vulpes for captive and wild
t.test(vulpes$T4.ug.dL~vulpes$Country) #non significant between Nth Africa and Iberian


#LINEAR MODEL - T4
# Tested T4 for linear model assumptions
fitT4 <- lm(T4.ug.dL~ Species2 + body.mass.Kg + BMI + Sex + Status, data = Foxbiochem2)


shapiro.test(as.numeric(residuals(object = fitT4))) # normal
bartlett.test(T4.ug.dL~Species2, data =Foxbiochem) # homogenous
plot(fitT4)

summary.lm(fitT4)
summary.aov(fitT4)
TukeyHSD(aov(fitT4), which = 'Species2') 

fitT42<-lm(T4.ug.dL~Species2 + aldosterone.pg.mL + copeptin.pmol.L + body.mass.Kg + BMI + Sex + Status, data = Foxbiochem2)
shapiro.test(as.numeric(residuals(object = fitT42))) # normal
summary.lm(fitT42)
summary.aov(fitT42)

### STATISTICS - METABOLITES
#UREA -> Not significant, not even logtransformed. BUT T4 was accused to be significant! 
pairwise.t.test(Foxbiochem2$urea.mg.dL, Foxbiochem2$Species2,p.adjust.method = "bonferroni") #NON SIGNIFICANT 
t.test(Foxbiochem2$urea.mg.dL~Foxbiochem2$Status) #non signficant for captive vs wild
t.test(Foxbiochem2$urea.mg.dL~Foxbioche2m$Age) #non signficant for adult vs juvenile
t.test(Foxbiochem2$urea.mg.dL~Foxbiochem2$Sex) #non signficant for M vs F
t.test(vulpes$urea.mg.dL~vulpes$Status) #non significant within V. vulpes for captive and wild

fiturea <-lm(urea.mg.dL~ Species2 + T4.ug.dL +body.mass.Kg + BMI + Sex + Status, data = Foxbiochem2)
summary.lm(fiturea)
summary.aov(fiturea)
shapiro.test(as.numeric(residuals(object = fiturea))) # not normal
leveneTest(urea.mg.dL~Species2, data =Foxbiochem2) # not homogeneous
bartlett.test(urea.mg.dL~Species2, data =Foxbiochem2) # not homogenous

fitlogurea <-lm(log(urea.mg.dL)~ Species2 + T4.ug.dL + body.mass.Kg + BMI + Sex + Status, data = Foxbiochem2) #NON SIGNIFICANT
summary.lm(fitlogurea)
summary.aov(fitlogurea)

shapiro.test(as.numeric(residuals(object = fitlogurea))) # normal
leveneTest(log(urea.mg.dL)~Species2, data =Foxbiochem2) # homogeneous
bartlett.test(log(urea.mg.dL)~Species2, data =Foxbiochem2) # homogenous

#URIC ACID -> Not significant (not even log transformed and accounting for other factors)
pairwise.t.test(Foxbiochem2$uric.acid.mg.dl, Foxbiochem2$Species2,p.adjust.method = "bonferroni") # NON-SIGNIFICANT
t.test(Foxbiochem2$uric.acid.mg.dl~Foxbiochem2$Status) #non signficant for captive vs wild
t.test(Foxbiochem2$uric.acid.mg.dl~Foxbiochem2$Sex) #non signficant for M vs F
t.test(vulpes$uric.acid.mg.dl~vulpes$Status) #non significant within V. vulpes for captive and wild

fituric<-lm(uric.acid.mg.dl~ Species2 + body.mass.Kg + BMI + Sex + Status, data = Foxbiochem2)
summary.lm(fituric)
summary.aov(fituric)
shapiro.test(as.numeric(residuals(object = fituric))) # not normal
leveneTest(uric.acid.mg.dl~Species2, data =Foxbiochem2) #  homogeneous
bartlett.test(uric.acid.mg.dl~Species2, data =Foxbiochem2) # homogenous

fitloguric<-lm(log(uric.acid.mg.dl)~ Species2 + body.mass.Kg + BMI + Sex + Status, data = Foxbiochem2) #NON SIGNIFICANT
summary.lm(fitloguric)
summary.aov(fitloguric)

shapiro.test(as.numeric(residuals(object = fitloguric))) # normal
leveneTest(log(uric.acid.mg.dl)~Species2, data =Foxbiochem2) # homogeneous
bartlett.test(log(uric.acid.mg.dl)~Species2, data =Foxbiochem2) # homogenous

#CREATININE - > Not significant (not even when log transformed and accounted for other factors)
pairwise.t.test(Foxbiochem2$creatinin.mg.dL, Foxbiochem2$Species2,p.adjust.method = "bonferroni") # NON-SIGNIFICANT
t.test(Foxbiochem2$creatinin.mg.dL~Foxbiochem2$Status) #non signficant for captive vs wild
t.test(Foxbiochem2$creatinin.mg.dL~Foxbiochem2$Sex) #non signficant for M vs F
t.test(vulpes$creatinin.mg.dL~vulpes$Status) #non significant within V. vulpes for captive and wild


fitcrea<-lm(creatinin.mg.dL~ Species2 + aldosterone.pg.mL + body.mass.Kg + BMI + Sex + Status, data = Foxbiochem2)
summary.lm(fitcrea)
summary.aov(fitcrea)
shapiro.test(as.numeric(residuals(object = fitcrea))) # not normal
leveneTest(creatinin.mg.dL~Species2, data =Foxbiochem2) #  homogeneous
bartlett.test(creatinin.mg.dL~Species2, data =Foxbiochem2) # homogenous

fitlogcrea<-lm(log(creatinin.mg.dL)~ Species2 + aldosterone.pg.mL + body.mass.Kg + BMI + Sex + Status, data = Foxbiochem2)
summary.lm(fitlogcrea)
summary.aov(fitlogcrea)
shapiro.test(as.numeric(residuals(object = fitlogcrea))) # normal
leveneTest(log(creatinin.mg.dL)~Species2, data =Foxbiochem2) #  homogeneous
bartlett.test(log(creatinin.mg.dL)~Species2, data =Foxbiochem2) # homogenous

### STATISTICS - PROTEINS AND FAT (CONTROLS, NOT EXPECTED TO CHANGE AMONG FOXES)

#ALBUMIN -> Not significant (no need to logtransform)
pairwise.t.test(Foxbiochem2$albumin.g.dL, Foxbiochem2$Species2,p.adjust.method = "bonferroni") # NOT SIGNIFICANT
t.test(Foxbiochem2$albumin.g.dL~Foxbiochem2$Status) #non signficant for captive vs wild
t.test(Foxbiochem2$albumin.g.dL~Foxbiochem2$Sex) #non signficant for M vs F
t.test(vulpes$albumin.g.dL~vulpes$Status) #non significant within V. vulpes for captive and wild

fitalb<-lm(albumin.g.dL~ Species2 + body.mass.Kg + BMI + Sex + Status, data = Foxbiochem2) #NOT SIGNIFICANT
summary.aov(fitalb)
plot(fitalb)
shapiro.test(as.numeric(residuals(object = fitalb))) # normal
leveneTest(albumin.g.dL~Species2, data =Foxbiochem2) #homogeneous
bartlett.test(albumin.g.dL~Species2, data =Foxbiochem2) # homogenous

#CHOLESTEROL -> Not significant (no need to logtransform)
pairwise.t.test(Foxbiochem2$cholesterol.mg.dL, Foxbiochem2$Species2,p.adjust.method = "bonferroni") # NOT SIGNIFICANT
t.test(Foxbiochem2$cholesterol.mg.dL~Foxbiochem2$Status) #non signficant for captive vs wild
t.test(Foxbiochem2$cholesterol.mg.dL~Foxbiochem2$Sex) #non signficant for M vs F
t.test(vulpes$cholesterol.mg.dL~vulpes$Status) #non significant within V. vulpes for captive and wild

fitchol<-lm(cholesterol.mg.dL~ Species2 + body.mass.Kg + BMI + Sex + Status, data = Foxbiochem2) #NOT SIGNIFICANT
summary.aov(fitchol)
shapiro.test(as.numeric(residuals(object = fitchol))) # normal
leveneTest(cholesterol.mg.dL~Species2, data =Foxbiochem2) #homogeneous
bartlett.test(cholesterol.mg.dL~Species2, data =Foxbiochem2) # homogenous

### STATISTICS - ELECTROLYTES AND PLASMA OSMOLALITY (NOT EXPECTED TO CHANGE AMONG FOXES)
t.test(Foxbiochem$sodium.mEq.L~Foxbiochem2$Species2) # NOT SIGNIFICANT
t.test(Foxbiochem$potassium.mEq.L~Foxbiochem2$Species2) # NOT SIGNIFICANT
t.test(Foxbiochem$chloride.mEq.L~Foxbiochem2$Species2) # NOT SIGNIFICANT
t.test(Foxbiochem$Osmolality.mOsm.Kg~Foxbiochem$Species2) # NOT SIGNIFICANT
