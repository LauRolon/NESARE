#Spot on lawn experiment
#LAB PS1155 and PS1156 against dairy-isolated L. monocytogenes

#By: MLR
#Last updated: 9/14/2022

setwd("C:/Users/lau_r/Google Drive/Penn State/Research/File for R/NESARE/Spot on lawn") #This is my working directory, change it for yours

#Input data (I removed the N/A so that thhe inhibition is read as numbers)
Data<-read.csv("Spot Data.csv", header = TRUE) #This is how you can attach a csv file. There are ways to import an Excel file, look into 'readxl' package

#Zones of inhibition were measured by Lm isolates, temperature and at two Lm lawn concentrations

#Attach libraries
library(ggplot2)
library(dplyr)
library(ggthemes)
library(tidyr)
library(agricolae)
library(psych)
library(svglite)



#Add column specifying combinations of temperature and concentration
Data$factorABC <- with(Data, interaction(PC, Temp, Concentration))

#Calculate statistics
Inhib<-describeBy(Data$`Inhibition`, group=list(  Data$Concentration, Data$Temp, Data$PC), mat = TRUE)


#ANOVA for the interaction effect (since I care about each independent variable)
anova_inhib<-aov(`Inhibition` ~ factorABC, data=Data)
summary(anova_inhib)


#tukey test
tukey_inhib<-HSD.test(anova_inhib, trt="factorABC") 
tukey_inhib



#Make dataframe with Tukey groups and order in the way as stat data
tukey_inhib_groups<-tukey_inhib$groups
tukey_inhib_groups$Sample<-as.character(rownames(tukey_inhib_groups))
tukey_inhib_groups<-tukey_inhib_groups[with(tukey_inhib_groups, order(Sample)),]
Inhib$TukeyGroups<-tukey_inhib_groups$groups

#Make plot
plot_summary<-ggplot(Inhib, aes(x=group2, y=mean, fill=group3))+
  geom_bar(colour="black",position="dodge", stat="identity", na.rm = TRUE)+  facet_grid(group1~group3,scales="free")+
  geom_errorbar(aes(ymin=mean-(sd/sqrt(n)), ymax=mean+(sd/sqrt(n))), width=0.25,position=position_dodge(5))+
  labs(x="Temperature (dC)",y="Average Inhibition(mm)")+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  theme(axis.text.x = element_text(angle=90, color='black', size=10))+
  labs(fill="LAB")+
  scale_y_continuous(limits=c(0,4))+theme(legend.text=element_text(size=10), legend.title= element_blank(), legend.position = 'bottom') +
  theme(axis.title=element_text(size=10),axis.text.y = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow=1)) +
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  ggtitle("Inhibition of Listeria monocytogenes by Protective Cultures")

ggsave("spot_summary.png", plot=plot_summary, device="png", width=5, height=4, units="in", dpi=600)
ggsave("spot_summary.svg", plot=plot_summary, device="svg", width=5, height=4, units="in", dpi=600)



#### Supplemental material ####
#Add column specifying combinations of temperature and concentration
Data$factorABCD <- with(Data, interaction(Isolate,PC, Temp, Concentration))

#Calculate statistics
Inhib_strain<-describeBy(Data$`Inhibition`, group=list(Data$Isolate,Data$Concentration, Data$Temp, Data$PC), mat = TRUE)


#ANOVA for the interaction effect (since I care about each independent variable)
anova_inhib_strain<-aov(`Inhibition` ~ factorABCD, data=Data)
summary(anova_inhib_strain)


#tukey test
tukey_inhib_strain<-HSD.test(anova_inhib_strain, trt="factorABCD") 
tukey_inhib_strain



#Make dataframe with Tukey groups and order in the way as stat data
tukey_inhib_groups_strain<-tukey_inhib_strain$groups
tukey_inhib_groups_strain$Sample<-as.character(rownames(tukey_inhib_groups_strain))
tukey_inhib_groups_strain<-tukey_inhib_groups_strain[with(tukey_inhib_groups_strain, order(Sample)),]

Inhib_strain<-Inhib_strain[order(Inhib_strain$group1,Inhib_strain$group4,Inhib_strain$group3,Inhib_strain$group2 ),]
Inhib_strain$TukeyGroups<-tukey_inhib_groups_strain$groups

#Make plot by temperature
plot_summary_strain<-ggplot(Inhib_strain, aes(x=group1, y=mean, fill=group3))+
  geom_bar(colour="black",position="dodge", stat="identity", na.rm = TRUE)+  facet_grid(group4+group2~group3,scales="free")+
  geom_errorbar(aes(ymin=mean-(sd/sqrt(n)), ymax=mean+(sd/sqrt(n))), width=0.1)+
  labs(x="Temperature (dC)",y="Average Inhibition(mm)")+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  theme(axis.text.x = element_text(angle=90, color='black', size=10))+
  labs(fill="L. monocytogenes strain")+
  scale_y_continuous(limits=c(0,4))+theme(legend.text=element_text(size=10), legend.title= element_blank(), legend.position = 'bottom') +
  theme(axis.title=element_text(size=10),axis.text.y = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow=1)) +
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=10),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_fill_brewer(palette='Dark2')+
  ggtitle("Inhibition of Listeria monocytogenes strains by LAB")
ggsave("spot_summary_strain.png", plot=plot_summary_strain, device="png", width=5, height=4, units="in", dpi=600)
ggsave("spot_summary_strain.svg", plot=plot_summary_strain, device="svg", width=5, height=4, units="in", dpi=600)


