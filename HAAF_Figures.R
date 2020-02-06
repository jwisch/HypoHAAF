library(ggplot2) # devtools::install_github("hadley/ggplot2")
library(ggalt)   # devtools::install_github("hrbrmstr/ggalt")
library(dplyr)   # for data_frame() & arrange()
library(tidyr)
library(gridExtra)

FILEPATH<-"C:/Users/julie.wisch/Documents/Hypo/"
MINBRAINREGIONINDEX<-5 #Update this if the columns change
MAXBRAINREGIONINDEX<-12 #Update this if the columns change


df <-read.csv(paste(FILEPATH, "MeanCBFFreeSurferReduced_20180919.csv", sep = ""), header = TRUE)
colnames(df)[3]<-"Glucose.level"
colnames(df)[1]<-"X"
df$Glucose.level<-factor(df$Glucose.level, ordered = TRUE, levels = c("90", "65", "55", "45"))
df$Day<-as.factor(df$Day)

#Drop all data more than 2 SD's from mean for the column
cutoff<-2

DropOutliers<-function(column, MEAN, cutoff, SD){
  Outliers<-ifelse((column <  MEAN - cutoff*SD) | (column > cutoff*SD + MEAN), 1, 0)
  List<-which(Outliers %in% 1) #gives indices of outliers
  for(i in 1:length(List)){
    column[List[i]]<-"NA"
  }
  column<-as.numeric(column)
  return(column)}


for(i in MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX){df[,i]<-DropOutliers(df[,i], mean(df[,i]), cutoff, sd(df[,i]))}

#if there are more than 2 na's for a single row, I drop the entire row
df$na_count <- apply(df, 1, function(x) sum(is.na(x)))
df<-df[(df$na_count < 3),]
df<-df[ , !names(df) %in% c("FD","na_count")] 
df$Group<-as.factor(df$Group)



############
df_long<-gather(df, Region, measurement, Frontal:InsulaAndAmygdala, factor_key=TRUE)
df_long$Cluster<-as.factor(paste(df_long$Group, df_long$Day))
levels(df_long$Cluster)<-c("Diabetic", "Control, Day 1", "Control, Day 2")
levels(df_long$Region)<-c("Frontal", "Orbital Frontal", "Motor", "Parietal",
                          "Temporal", "Occipital", "Basal Ganglia", "Thalamus", "Insula and Amygdala")

p1<-ggplot(df_long, aes(x=Glucose.level, y=measurement, fill=Cluster)) + 
  geom_boxplot() + facet_wrap(~Region, scales="free_y") + 
  xlab("Glucose Level") + ylab("Cerebral Blood Flow") +
  scale_fill_manual(values = c("#e66101", "#5e3c99", "#b2abd2")) +
  labs(fill = "") + 
  theme(legend.position = "bottom", strip.background = element_rect(fill="black"),
        strip.text = element_text(colour = 'white'),
        legend.background = element_rect(fill= "white"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey90"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey90")) 


se <- function(x) sqrt(var(x)/length(x))

agg = aggregate(df_long,
                by = list(df_long$Cluster, df_long$Glucose.level, df_long$Region),
                FUN = mean, na.rm = TRUE)
agg2 = aggregate(df_long,
                 by = list(df_long$Cluster, df_long$Glucose.level, df_long$Region),
                 FUN = se)

pointplot_df<-data.frame(cbind(agg[,c("Group.1", "Group.2", "Group.3", "measurement")], agg2[,"measurement"]))
rm(agg, agg2)
colnames(pointplot_df)<-c("Cluster", "Glucose Level", "Region", "Mean", "Standard Error")
pointplot_df$lower<-pointplot_df$Mean - 2*pointplot_df$`Standard Error`
pointplot_df$upper<-pointplot_df$Mean + 2*pointplot_df$`Standard Error`


p2<-ggplot(pointplot_df, aes(x=`Glucose Level`, y=Mean, colour=Cluster, group = Cluster)) + 
  geom_point() + geom_line() + facet_wrap(~Region) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) +
  xlab("Glucose Level") + ylab("Cerebral Blood Flow\n+/-2 Standard Error") +
  scale_colour_manual(values = c("#e66101", "#5e3c99", "#b2abd2")) +
  labs(colour = "") + 
  theme(legend.position = "bottom", strip.background = element_rect(fill="black"),
        strip.text = element_text(colour = 'white'),
        legend.background = element_rect(fill= "white"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey90"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey90")) 



hormones<-read.csv(paste(FILEPATH, "Hormones_Updated.csv", sep = ""))
hormones <- gather(hormones, ID, measurement, BOLD.01:BOLD.45, factor_key=TRUE)
df$X<-gsub("-", ".", df$X)

hormones$Day<-gsub("Day ", "", hormones$Day)

hormones$HYA<-factor(hormones$HYA, ordered = TRUE, levels = c("90", "65", "55", "45"))
hormones$Day<-as.factor(hormones$Day)

hormones<-merge(hormones, df, by.x = c("ID", "Day", "HYA"), by.y = c("X", "Day", "Glucose.level"))
hormones$Cluster<-as.factor(paste(hormones$Group, hormones$Day, sep = "-"))
levels(hormones$Cluster)<-c("Diabetic", "Control, Day 1", "Control, Day 2")

p3<-ggplot(hormones, aes(x=HYA, y=measurement, fill=Cluster)) + 
  geom_boxplot() + facet_wrap(~Hormone, scales="free_y") + 
  xlab("Glucose Level") + ylab("Cerebral Blood Flow") +
  scale_fill_manual(values = c("#e66101", "#5e3c99", "#b2abd2")) +
  labs(fill = "") + 
  theme(legend.position = "bottom", strip.background = element_rect(fill="black"),
        strip.text = element_text(colour = 'white'),
        legend.background = element_rect(fill= "white"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey90"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey90")) 

#Need to go long to wide for hormones
hormone_wide<-spread(hormones[,1:7], Hormone, measurement)
hormone_wide1<-hormone_wide[!is.na(hormone_wide[,"C-pep"]), c("ID", "Day", "HYA", "Group", "C-pep")]
hormone_wide2<-hormone_wide[!is.na(hormone_wide[,"cortisol"]), c("ID", "Day", "HYA", "Group", "cortisol")]
hormone_wide3<-hormone_wide[!is.na(hormone_wide[,"EPI"]), c("ID", "Day", "HYA", "Group", "EPI")]
hormone_wide4<-hormone_wide[!is.na(hormone_wide[,"ffa"]), c("ID", "Day", "HYA", "Group", "ffa")]
hormone_wide5<-hormone_wide[!is.na(hormone_wide[,"glucagon"]), c("ID", "Day", "HYA", "Group", "glucagon")]
hormone_wide6<-hormone_wide[!is.na(hormone_wide[,"glucose" ]), c("ID", "Day", "HYA", "Group", "glucose" )]
hormone_wide7<-hormone_wide[!is.na(hormone_wide[,"insulin"]), c("ID", "Day", "HYA", "Group", "insulin")]
hormone_wide8<-hormone_wide[!is.na(hormone_wide[,"lactate"]), c("ID", "Day", "HYA", "Group", "lactate")]
hormone_wide9<-hormone_wide[!is.na(hormone_wide[,"NG"]), c("ID", "Day", "HYA", "Group", "NG")]
hormone_wide10<-hormone_wide[!is.na(hormone_wide[,"NGP"]), c("ID", "Day", "HYA", "Group", "NGP")]
hormone_wide11<-hormone_wide[!is.na(hormone_wide[,"Norepi"]), c("ID", "Day", "HYA", "Group", "Norepi")]
hormone_wide12<-hormone_wide[!is.na(hormone_wide[,"TSS"]), c("ID", "Day", "HYA", "Group", "TSS")]

hormone_wide<-merge(hormone_wide1, hormone_wide2, by = c("ID", "Day", "HYA", "Group"))
hormone_wide<-merge(hormone_wide, hormone_wide3, by = c("ID", "Day", "HYA", "Group"))
hormone_wide<-merge(hormone_wide, hormone_wide4, by = c("ID", "Day", "HYA", "Group"))
hormone_wide<-merge(hormone_wide, hormone_wide5, by = c("ID", "Day", "HYA", "Group"))
hormone_wide<-merge(hormone_wide, hormone_wide6, by = c("ID", "Day", "HYA", "Group"))
hormone_wide<-merge(hormone_wide, hormone_wide7, by = c("ID", "Day", "HYA", "Group"))
hormone_wide<-merge(hormone_wide, hormone_wide8, by = c("ID", "Day", "HYA", "Group"))
hormone_wide<-merge(hormone_wide, hormone_wide9, by = c("ID", "Day", "HYA", "Group"))
hormone_wide<-merge(hormone_wide, hormone_wide10, by = c("ID", "Day", "HYA", "Group"))
hormone_wide<-merge(hormone_wide, hormone_wide11, by = c("ID", "Day", "HYA", "Group"))
hormone_wide<-merge(hormone_wide, hormone_wide12, by = c("ID", "Day", "HYA", "Group"))
rm(hormone_wide1, hormone_wide2,hormone_wide3,hormone_wide4,hormone_wide5,hormone_wide6,
   hormone_wide7,hormone_wide8,hormone_wide9,hormone_wide10,hormone_wide11,hormone_wide12)

hormone_wide$Cluster<-paste(hormone_wide$Group, hormone_wide$Day, sep = "-")

se <- function(x) sqrt(var(x)/length(x))

agg = aggregate(.~Cluster + Hormone + HYA,
                data = hormones[,c("measurement", "Cluster", "Hormone", "HYA")],
                FUN = mean, na.rm = TRUE)
agg2 = aggregate(.~Cluster + Hormone + HYA,
                data = hormones[,c("measurement", "Cluster", "Hormone", "HYA")],
                FUN = se)

agg<-merge(agg, agg2, by = c("Cluster", "Hormone", "HYA"))




p4<-ggplot(agg, aes(x=HYA, y= measurement.x, colour=Cluster, group = Cluster)) + 
  geom_point() + geom_line() + facet_wrap(~Hormone, scales="free_y") + 
  geom_errorbar(aes(ymin=measurement.x - 2*measurement.y, ymax=measurement.x + 2*measurement.y), width=.1) +
  xlab("Glucose Level") + ylab("Hormone level\n+/-2 Standard Error") +
  scale_colour_manual(values = c("#e66101", "#5e3c99", "#b2abd2")) +
  labs(colour = "") + 
  theme(legend.position = "bottom", strip.background = element_rect(fill="black"),
        strip.text = element_text(colour = 'white'),
        legend.background = element_rect(fill= "white"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey90"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey90")) 

pdf(paste(FILEPATH, "HypoHAAF_Figures.pdf", sep = ""), width = 8.5, height = 11)
p1
p2
p3
p4
dev.off()
