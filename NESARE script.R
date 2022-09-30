#NESARE project 


#Laura Rolon
#Last updated: 8/16/2022

#Load packages
library(ape)
library(vegan)
library(ggplot2)
library(phyloseq)
library(cowplot)
library(tidyr)
library(dplyr)
library(compositions)
library(zCompositions)
library(viridis)
library(readxl)
library(pairwiseAdonis)
library(SpadeR)
library(psych)
library(svglite)
library(agricolae)
library(decontam)


set.seed(336)

#### Control experiment 1 ####
#
setwd("C:/Users/lau_r/Google Drive/Penn State/Research/File for R/NESARE")

#Import data
exp1<-read_excel("Controls - Result sheet.xlsx", sheet=5, col_names = TRUE)

#Calculate statistics
exp1_logAPC<-describeBy(exp1$logAPC, group=exp1$Treatment, mat = TRUE)
exp1_logAPC$order<-c(2,3,1)
  
#ANOVA
#APC
exp1_anova_APC<-aov(logAPC ~ Treatment, data=exp1)
summary(exp1_anova_APC)

#tukey test
exp1_tukey_APC_all<-HSD.test(exp1_anova_APC, trt="Treatment") 

#Make dataframe with Tukey groups
tukey_APC_exp1_groups<-exp1_tukey_APC_all$groups
tukey_APC_exp1_groups$Sample<-as.character(rownames(tukey_APC_exp1_groups))
tukey_APC_exp1_groups_order<-tukey_APC_exp1_groups[with(tukey_APC_exp1_groups, order(Sample)),]
exp1_logAPC$TukeyGroups<-tukey_APC_exp1_groups_order$groups


#Plot
exp1_apc<- ggplot(exp1_logAPC, aes(x = reorder(group1,order), y = mean , fill = reorder(group1,order)))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.1,position=position_dodge(.4))+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  theme(legend.text=element_text(size=9), legend.title= element_blank(), legend.position = 'bottom') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow=1)) +
  ylab("log CFU/ml") + xlab("Treatment")+
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Control exp. 1 - Aerobic Plate count")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_y_continuous(limits=c(0,10) ,breaks= seq(0,10, 2))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='viridis', alpha = 1)
exp1_apc
ggsave("exp1_apc.png", plot=exp1_apc, device="png", width=4, height=5, units="in", dpi=600)
ggsave("exp1_apc.svg", plot=exp1_apc, device="svg", width=4, height=5, units="in", dpi=600)

#### Control experiment 2 ####

#Import data
exp2<-read_excel("Controls - Result sheet.xlsx", sheet=6, col_names = TRUE)

#Calculate statistics
exp2_logAPC<-describeBy(exp2$logAPC, group=exp2$Treatment, mat = TRUE)
exp2_logAPC$order<-c(4,1,2,3)

exp2_logMPN<-describeBy(exp2$logMPN, group=exp2$Treatment, mat = TRUE)
exp2_logMPN$order<-c(4,1,2,3)

#ANOVA
exp2_anova_APC<-aov(logAPC ~ Treatment, data=exp2)
summary(exp2_anova_APC)

exp2_anova_MPN<-aov(logMPN ~ Treatment, data=exp2)
summary(exp2_anova_MPN)

#tukey test
exp2_tukey_APC_all<-HSD.test(exp2_anova_APC, trt="Treatment") 
exp2_tukey_MPN_all<-HSD.test(exp2_anova_MPN, trt="Treatment") 

#Make dataframe with Tukey groups
tukey_APC_exp2_groups<-exp2_tukey_APC_all$groups
tukey_APC_exp2_groups$Sample<-as.character(rownames(tukey_APC_exp2_groups))
tukey_APC_exp2_groups_order<-tukey_APC_exp2_groups[with(tukey_APC_exp2_groups, order(Sample)),]
exp2_logAPC$TukeyGroups<-tukey_APC_exp2_groups_order$groups

tukey_MPN_exp2_groups<-exp2_tukey_MPN_all$groups
tukey_MPN_exp2_groups$Sample<-as.character(rownames(tukey_MPN_exp2_groups))
tukey_MPN_exp2_groups_order<-tukey_MPN_exp2_groups[with(tukey_MPN_exp2_groups, order(Sample)),]
exp2_logMPN$TukeyGroups<-tukey_MPN_exp2_groups_order$groups

#Plots
exp2_apc<- ggplot(exp2_logAPC, aes(x = reorder(group1,order), y = mean , fill = reorder(group1,order)))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.1,position=position_dodge(.4))+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  theme(legend.text=element_text(size=9), legend.title= element_blank(), legend.position = 'bottom') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow=1)) +
  ylab("log CFU/ml") + xlab("Treatment")+
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Control exp. 2 - Aerobic Plate count")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_y_continuous(limits=c(0,10) ,breaks= seq(0,10, 2))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='viridis', alpha = 1)
exp2_apc
ggsave("exp2_apc.png", plot=exp2_apc, device="png", width=4, height=5, units="in", dpi=600)
ggsave("exp2_apc.svg", plot=exp2_apc, device="svg", width=4, height=5, units="in", dpi=600)

exp2_mpn<- ggplot(exp2_logMPN, aes(x = reorder(group1,order), y = mean , fill = reorder(group1,order)))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.1,position=position_dodge(.4))+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  theme(legend.text=element_text(size=9), legend.title= element_blank(), legend.position = 'bottom') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow=1)) +
  ylab("log MPN/ml") + xlab("Treatment")+
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Control exp. 2 - MPN")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_y_continuous(limits=c(0,10) ,breaks= seq(0,10, 2))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='viridis', alpha = 1)
exp2_mpn
ggsave("exp2_mpn.png", plot=exp2_mpn, device="png", width=4, height=5, units="in", dpi=600)
ggsave("exp2_mpn.svg", plot=exp2_mpn, device="svg", width=4, height=5, units="in", dpi=600)


#### Obj2: Attachment of lactic acid bacteria in the presence of microbiomes ####

#Set working directory to where files are located
setwd("C:/Users/lau_r/Google Drive/Penn State/Research/File for R/NESARE/Obj2")

obj2<-read_excel('Obj2 - Result sheet.xlsx', sheet=1, col_names=TRUE)

#Calculate statistics
logAPC<-describeBy(obj2$logAPC, group=obj2$Sample, mat = TRUE)
logAPC$Facility<-c(rep("A", 5),rep("B",5),rep("C",5))
logAPC$Treatment<-rep(c("Initial","ED","LE","LL","NC"),3)
logAPC$order<-rep(c(0,3,4,2,1),3)

logLAC<-describeBy(obj2$logLAC, group=obj2$Sample, mat = TRUE) 
logLAC$Facility<-c(rep("A", 5),rep("B",5),rep("C",5))
logLAC$Treatment<-rep(c("Initial","ED","LE","LL","NC"),3)
logLAC$order<-rep(c(0,3,4,2,1),3)

#ANOVA
#APC
obj2_anova_APC_all<-aov(logAPC ~ Sample, data=obj2)
summary(obj2_anova_APC_all)

#tukey test
tukey_APC_all<-HSD.test(obj2_anova_APC_all, trt="Sample") 

#Make dataframe with Tukey groups
tukey_APC_all_groups<-tukey_APC_all$groups
tukey_APC_all_groups$Sample<-as.character(rownames(tukey_APC_all_groups))
tukey_APC_all_groups_order<-tukey_APC_all_groups[with(tukey_APC_all_groups, order(Sample)),]
logAPC$TukeyGroups<-tukey_APC_all_groups_order$groups

#LAC
obj2_anova_LAC_all<-aov(logLAC ~ Sample, data=obj2)
summary(obj2_anova_LAC_all)

#Tukey test
tukey_LAC_all<-HSD.test(obj2_anova_LAC_all, trt="Sample") 

tukey_LAC_all_groups<-tukey_LAC_all$groups
tukey_LAC_all_groups$Sample<-as.character(rownames(tukey_LAC_all_groups))
tukey_LAC_all_groups_order<-tukey_LAC_all_groups[with(tukey_LAC_all_groups, order(Sample)),]
logLAC$TukeyGroups<-tukey_LAC_all_groups_order$groups

#Plots
obj2_apc<- ggplot(logAPC, aes(x = reorder(Treatment,order), y = mean , fill = Facility))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) + 
  facet_grid(Facility~.)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.1,position=position_dodge(.4))+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  theme(legend.text=element_text(size=9), legend.title= element_blank(), legend.position = 'bottom') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow=1)) +
  ylab("log CFU/ml") + xlab("Treatment")+
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Obj2 - Aerobic Plate count")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_y_continuous(limits=c(0,10) ,breaks= seq(0,10, 2))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno', alpha = 0.8)
obj2_apc
ggsave("Obj2_apc.png", plot=obj2_apc, device="png", width=4, height=5, units="in", dpi=600)
ggsave("Obj2_apc.svg", plot=obj2_apc, device="svg", width=4, height=5, units="in", dpi=600)


obj2_lac<- ggplot(logLAC, aes(x = reorder(Treatment,order), y = mean , fill = Facility))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) + 
  facet_grid(Facility~.)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.1,position=position_dodge(.4))+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  theme(legend.text=element_text(size=9), legend.title= element_blank(), legend.position = 'bottom') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow=1)) +
  ylab("log CFU/ml") + xlab("Treatment")+
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Obj2 - Lactic acid bacteria count")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_y_continuous(limits=c(0,10) ,breaks= seq(0,10, 2))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno', alpha = 0.8)
obj2_lac
ggsave("Obj2_lac.png", plot=obj2_lac, device="png", width=4, height=5, units="in", dpi=600)
ggsave("Obj2_lac.svg", plot=obj2_lac, device="svg", width=4, height=5, units="in", dpi=600)



#ANOVA by facility
#Subset by facility
obj2_A<-subset(obj2, Facility=="A")
obj2_B<-subset(obj2, Facility=="B")
obj2_C<-subset(obj2, Facility=="C")


#APC-A
obj2_anova_APC_A<-aov(logAPC ~ Sample, data=obj2_A)
summary(obj2_anova_APC_A)

tukey_APC_A<-TukeyHSD(obj2_anova_APC_A) #Initial is different from NC, NC is different from LL, ED, LE
write.csv(as.data.frame(tukey_APC_A[1:1]), file='tukeyAPC_A.csv')

#APC-B
obj2_anova_APC_B<-aov(logAPC ~ Sample, data=obj2_B)
summary(obj2_anova_APC_B)

tukey_APC_B<-TukeyHSD(obj2_anova_APC_B)
write.csv(as.data.frame(tukey_APC_B[1:1]), file='tukeyAPC_B.csv')

#APC-C
obj2_anova_APC_C<-aov(logAPC ~ Sample, data=obj2_C)
summary(obj2_anova_APC_C)

tukey_APC_C<-TukeyHSD(obj2_anova_APC_C)
write.csv(as.data.frame(tukey_APC_C[1:1]), file='tukeyAPC_C.csv')

#LAC-A
obj2_anova_LAC_A<-aov(logLAC ~ Sample, data=obj2_A)
summary(obj2_anova_LAC_A)

tukey_LAC_A<-TukeyHSD(obj2_anova_LAC_A)
write.csv(as.data.frame(tukey_LAC_A[1:1]), file='tukeyLAC_A.csv')

#LAC-B
obj2_anova_LAC_B<-aov(logLAC ~ Sample, data=obj2_B)
summary(obj2_anova_LAC_B)

tukey_LAC_B<-TukeyHSD(obj2_anova_LAC_B)
write.csv(as.data.frame(tukey_LAC_B[1:1]), file='tukeyLAC_B.csv')

#LAC-C
obj2_anova_LAC_C<-aov(logLAC ~ Sample, data=obj2_C)
summary(obj2_anova_LAC_C)

tukey_LAC_C<-TukeyHSD(obj2_anova_LAC_C)
write.csv(as.data.frame(tukey_LAC_C[1:1]), file='tukeyLAC_C.csv')



#### Obj3: Inhibition of L. monocytogenes by lactic acid bacteria in the presence of microbiomes ####

#Set working directory to where files are located
setwd("C:/Users/lau_r/Google Drive/Penn State/Research/File for R/NESARE/Obj3")

obj3<-read_excel('Obj3 - Result sheet.xlsx', sheet=1, col_names=TRUE)

#Calculate statistics
logAPC_obj3<-describeBy(obj3$logAPC, group=obj3$Sample, mat = TRUE)
logAPC_obj3$Facility<-c(rep("A", 6),rep("B",6),rep("C",6))
logAPC_obj3$Treatment<-rep(c("Initial","ED","LE","LL","NC","PC"),3)
logAPC_obj3$order<-rep(c(0,4,5,3,1,2),3)

logLAC_obj3<-describeBy(obj3$logLAC, group=obj3$Sample, mat = TRUE) 
logLAC_obj3$Facility<-c(rep("A", 6),rep("B",6),rep("C",6))
logLAC_obj3$Treatment<-rep(c("Initial","ED","LE","LL","NC","PC"),3)
logLAC_obj3$order<-rep(c(0,4,5,3,1,2),3)

logLM_obj3<-describeBy(obj3$logMPN, group=obj3$Sample, mat = TRUE) 
logLM_obj3$Facility<-c(rep("A", 6),rep("B",6),rep("C",6))
logLM_obj3$Treatment<-rep(c("Initial","ED","LE","LL","NC","PC"),3)
logLM_obj3$order<-rep(c(0,4,5,3,1,2),3)

#ANOVA
#APC
obj3_anova_APC_all<-aov(logAPC~ Sample, data=obj3)
summary(obj3_anova_APC_all)

#tukey test
obj3_tukey_APC_all<-HSD.test(obj3_anova_APC_all, trt="Sample") 

#Make dataframe with Tukey groups
obj3_tukey_APC_all_groups<-obj3_tukey_APC_all$groups
obj3_tukey_APC_all_groups$Sample<-as.character(rownames(obj3_tukey_APC_all_groups))
obj3_tukey_APC_all_groups_order<-obj3_tukey_APC_all_groups[with(obj3_tukey_APC_all_groups, order(Sample)),]
logAPC_obj3$TukeyGroups<-obj3_tukey_APC_all_groups_order$groups

#LAC
obj3_anova_LAC_all<-aov(logLAC~ Sample, data=obj3)
summary(obj3_anova_LAC_all)

#tukey test
obj3_tukey_LAC_all<-HSD.test(obj3_anova_LAC_all, trt="Sample") 

#Make dataframe with Tukey groups
obj3_tukey_LAC_all_groups<-obj3_tukey_LAC_all$groups
obj3_tukey_LAC_all_groups$Sample<-as.character(rownames(obj3_tukey_LAC_all_groups))
obj3_tukey_LAC_all_groups_order<-obj3_tukey_LAC_all_groups[with(obj3_tukey_LAC_all_groups, order(Sample)),]
logLAC_obj3$TukeyGroups<-obj3_tukey_LAC_all_groups_order$groups

#LM
obj3_anova_LM_all<-aov(logMPN ~ Sample, data=obj3)
summary(obj3_anova_LM_all)

#tukey test
obj3_tukey_LM_all<-HSD.test(obj3_anova_LM_all, trt="Sample") 

#Make dataframe with Tukey groups
obj3_tukey_LM_all_groups<-obj3_tukey_LM_all$groups
obj3_tukey_LM_all_groups$Sample<-as.character(rownames(obj3_tukey_LM_all_groups))
obj3_tukey_LM_all_groups_order<-obj3_tukey_LM_all_groups[with(obj3_tukey_LM_all_groups, order(Sample)),]
logLM_obj3$TukeyGroups<-obj3_tukey_LM_all_groups_order$groups

#Plots
obj3_apc<- ggplot(logAPC_obj3, aes(x = reorder(Treatment,order), y = mean , fill = Facility))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) + 
  facet_grid(Facility~.)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.1,position=position_dodge(.4))+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  theme(legend.text=element_text(size=9), legend.title= element_blank(), legend.position = 'bottom') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow=1)) +
  ylab("log CFU/ml") + xlab("Treatment")+
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("obj3 - Aerobic Plate count")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_y_continuous(limits=c(0,10) ,breaks= seq(0,10, 2))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno', alpha = 0.8)
obj3_apc
ggsave("obj3_apc.png", plot=obj3_apc, device="png", width=4, height=5, units="in", dpi=600)
ggsave("obj3_apc.svg", plot=obj3_apc, device="svg", width=4, height=5, units="in", dpi=600)


obj3_lac<- ggplot(logLAC_obj3, aes(x = reorder(Treatment,order), y = mean , fill = Facility))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) + 
  facet_grid(Facility~.)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.1,position=position_dodge(.4))+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  theme(legend.text=element_text(size=9), legend.title= element_blank(), legend.position = 'bottom') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow=1)) +
  ylab("log CFU/ml") + xlab("Treatment")+
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("obj3 - Lactic acid bacteria count")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_y_continuous(limits=c(0,10) ,breaks= seq(0,10, 2))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno', alpha = 0.8)
obj3_lac
ggsave("obj3_lac.png", plot=obj3_lac, device="png", width=4, height=5, units="in", dpi=600)
ggsave("obj3_lac.svg", plot=obj3_lac, device="svg", width=4, height=5, units="in", dpi=600)

obj3_LM<- ggplot(logLM_obj3, aes(x = reorder(Treatment,order), y = mean , fill = Facility))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) + 
  facet_grid(Facility~.)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.1,position=position_dodge(.4))+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  geom_hline(yintercept=1.85, linetype=2, color="grey20")+
  theme(legend.text=element_text(size=9), legend.title= element_blank(), legend.position = 'bottom') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow=1)) +
  ylab("log MPN/sample") + xlab("Treatment")+
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("obj3 - Log MPN L. monocytogenes")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_y_continuous(limits=c(0,10) ,breaks= seq(0,10, 2))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno', alpha = 0.8)
obj3_LM
ggsave("obj3_LM.png", plot=obj3_LM, device="png", width=4, height=5, units="in", dpi=600)
ggsave("obj3_LM.svg", plot=obj3_LM, device="svg", width=4, height=5, units="in", dpi=600)


#Log reduction
obj3_red<-read_excel('Obj3 - Result sheet.xlsx', sheet=7, col_names=TRUE)

#Calculate mean by treatment and Facility
obj3_red_mean<-obj3_red %>%
  group_by(Facility, Treatment)%>%
  summarize(Mean=mean(LogRed), SD=sd(LogRed))%>%
  arrange(desc(Mean))

obj3_red_mean$SE<-obj3_red_mean$SD/sqrt(4)

#Anova
obj3_anova_red<-aov(Mean ~ Treatment+Facility, data=obj3_red_mean)
summary(obj3_anova_red)
#Note: No significant difference by treatment

tukey_red_fac<-HSD.test(obj3_anova_red, trt="Facility") #All facilities are sign. different
tukey_red_fac

#Plot Lm reductions
obj3_red<- ggplot(obj3_red_mean, aes(x = Treatment, y = Mean , fill = Facility))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) +
  facet_grid(Facility~.)+
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE),width=.1,position=position_dodge(.4))+
  #geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  geom_hline(yintercept = 1.82, type="dashed", color='grey90')+
  theme(legend.text=element_text(size=9), legend.title= element_blank(), legend.position = 'bottom') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow=1)) +
  ylab("log CFU/ml") + xlab("Treatment")+
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Obj3 - Lm Log reduction")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_y_continuous(limits=c(-4,6) ,breaks= seq(-4,6, 2))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno', alpha = 0.8)
obj3_red
ggsave("Obj3_red.png", plot=obj3_red, device="png", width=4, height=5, units="in", dpi=600)
ggsave("Obj3_red.svg", plot=obj3_red, device="svg", width=4, height=5, units="in", dpi=600)


#Compare log reductions from full stregnth and dilution experiment of Facility B
#Log reduction
dil_red<-read_excel('Obj3 - Result sheet.xlsx', sheet=2, col_names=TRUE)

#Calculate mean by treatment and Facility
dil_red_mean<-dil_red %>%
  group_by(Sample)%>%
  summarize(Mean=mean(Lm_red), SD=sd(Lm_red))%>%
  arrange(desc(Mean))

dil_red_mean$SE<-dil_red_mean$SD/sqrt(4)

#Anova
dil_anova_red<-aov(Lm_red ~ Sample, data=dil_red)
summary(dil_anova_red)
#Note: No significant difference by treatment


#Plot Lm reductions
dil_red<- ggplot(dil_red_mean, aes(x = Sample, y = Mean))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE),width=.1,position=position_dodge(.4))+
  #geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  geom_hline(yintercept = 0, color='black')+
  theme(legend.text=element_text(size=9), legend.title= element_blank(), legend.position = 'bottom') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow=1)) +
  ylab("log MPN/ml") + xlab("Treatment")+
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Dil - Lm Log reduction")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_y_continuous(limits=c(-4,6) ,breaks= seq(-4,6, 2))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno', alpha = 0.8)
dil_red
ggsave("Dilution experiment_reduction.png", plot=dil_red, device="png", width=4, height=5, units="in", dpi=600)
ggsave("Dilution experiment_reduction.svg", plot=obj3_red, device="svg", width=4, height=5, units="in", dpi=600)


