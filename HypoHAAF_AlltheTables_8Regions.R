library(DescTools)
library(stargazer)
library(lmerTest)

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



hormones<-read.csv(paste(FILEPATH, "Hormones_Updated.csv", sep = ""))
hormones <- gather(hormones, ID, measurement, BOLD.01:BOLD.45, factor_key=TRUE)
df$X<-gsub("-", ".", df$X)

hormones$Day<-gsub("Day ", "", hormones$Day)

hormones$HYA<-factor(hormones$HYA, ordered = TRUE, levels = c("90", "65", "55", "45"))
hormones$Day<-as.factor(hormones$Day)

hormones<-merge(hormones, df, by.x = c("ID", "Day", "HYA"), by.y = c("X", "Day", "Glucose.level"))
############
#Linear model day 1 control vs day 1 diabetic
#In these models, 90 is the reference glucose level
#Controls are the reference group
df$Group <- relevel(df$Group, ref="1")
hormones$Group <- relevel(hormones$Group, ref="1")

#Doing models for CBF
modelList<-list()
for(i in MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX){
  modelList[[i-(MINBRAINREGIONINDEX - 1)]] <- lm(df[df$Day == 1,i] ~ Group + Glucose.level + Group:Glucose.level, data = df[df$Day == 1,])}

stargazer(modelList[[1]], modelList[[2]], modelList[[3]], modelList[[4]],
          modelList[[5]], modelList[[6]], modelList[[7]], modelList[[8]],
          title = "Linear Models - Day 1 Controls vs. Diabetics",
          align= TRUE, type = "html", out = paste(FILEPATH, "LM_Day1ControlvsDiabetic.htm"),
                    omit.stat=c("LL","ser"),
          dep.var.labels = "Linear Regression Model",
          column.labels = c("Frontal", "Orbital Frontal", "Motor", "Parietal", "Temporal", 
                             "Occipital", "Basal Ganglia", "Thalamus"),
          covariate.labels=c("Diabetic", "Glucose Level - 65", "Glucose Level - 55", "Glucose Level - 45",
                             "Diabetic: Glucose Level - 65", "Diabetic: Glucose Level - 55",
                             "Diabetic: Glucose Level - 45", "Intercept"),
          ci =TRUE, star.cutoffs = c(0.05/8, 0.01/8, 0.001/8)) #Stars are corrected for multiple comparisons

modelList<-list()
for(i in MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX){
  modelList[[i-(MINBRAINREGIONINDEX - 1)]] <- lm(df[(df$Day == 1 & df$Group == 0) | df$Day == 2,i] ~ Group + Glucose.level + Group:Glucose.level, data = df[(df$Day == 1 & df$Group == 0) | df$Day == 2,])}

stargazer(modelList[[1]], modelList[[2]], modelList[[3]], modelList[[4]],
          modelList[[5]], modelList[[6]], modelList[[7]], modelList[[8]],
          title = "Linear Models - Day 2 Controls vs. Diabetics",
          align= TRUE, type = "html", out = paste(FILEPATH, "LM_Day2ControlvsDiabetic.htm"),
          omit.stat=c("LL","ser"),
          dep.var.labels = "Linear Regression Model",
          column.labels = c("Frontal", "Orbital Frontal", "Motor", "Parietal", "Temporal", 
                            "Occipital", "Basal Ganglia", "Thalamus"),
          covariate.labels=c("Diabetic", "Glucose Level - 65", "Glucose Level - 55", "Glucose Level - 45",
                             "Diabetic: Glucose Level - 65", "Diabetic: Glucose Level - 55",
                             "Diabetic: Glucose Level - 45", "Intercept"),
          ci =TRUE, star.cutoffs = c(0.05/8, 0.01/8, 0.001/8)) #Stars are corrected for multiple comparisons

#Linear model comparing controls on day 1 vs. day 2
modelList<-list()
for(i in MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX){
  modelList[[i-(MINBRAINREGIONINDEX - 1)]] <- lmerTest::lmer(df[df$Group == 1,i] ~ Day + Glucose.level + Day:Glucose.level + (1|X), data = df[df$Group == 1,])
  class(modelList[[i-(MINBRAINREGIONINDEX - 1)]]) <- "lmerMod"} #Have to change the class of the model to get it to work with stargazer

stargazer(modelList[[1]], modelList[[2]], modelList[[3]], modelList[[4]],
          modelList[[5]], modelList[[6]], modelList[[7]], modelList[[8]],
          title = "Linear Models - Day 1 Controls vs. Day 2 Controls",
          align= TRUE, type = "html", out = paste(FILEPATH, "LM_Day2ControlvsDay1Controls.htm"),
          omit.stat=c("LL","ser"),
          dep.var.labels = "Linear Mixed Effect Model",
          column.labels = c("Frontal", "Orbital Frontal", "Motor", "Parietal", "Temporal", 
                            "Occipital", "Basal Ganglia", "Thalamus"),
          covariate.labels=c("Day 2", "Glucose Level - 65", "Glucose Level - 55", "Glucose Level - 45",
                             "Day 2: Glucose Level - 65", "Day 2: Glucose Level - 55",
                             "Day 2: Glucose Level - 45", "Intercept"),
          ci =TRUE, star.cutoffs = c(0.05/8, 0.01/8, 0.001/8)) #Stars are corrected for multiple comparisons


###<------Doing models for hormones
modelList<-list()
HORMONES<-(unique(hormones$Hormone))
for(i in 1:length(HORMONES)){
  modelList[[i]] <- lm(hormones[hormones$Day == 1,HORMONES[i]] ~ Group + HYA + Group:HYA, data = hormones[hormones$Day == 1,])}

stargazer(modelList[[1]], modelList[[2]], modelList[[3]], modelList[[4]],
          modelList[[5]], modelList[[6]], modelList[[7]], modelList[[8]],
          modelList[[9]], modelList[[10]], modelList[[11]], modelList[[12]],
          title = "Linear Models - Day 1 Controls vs. Diabetics",
          align= TRUE, type = "html", out = paste(FILEPATH, "LM_Day1ControlvsDiabetic_Hormones.htm"),
          omit.stat=c("LL","ser"),
          dep.var.labels = "Linear Regression Model",
          column.labels = as.character(HORMONES),
          covariate.labels=c("Diabetic", "Glucose Level - 65", "Glucose Level - 55", "Glucose Level - 45",
                             "Diabetic: Glucose Level - 65", "Diabetic: Glucose Level - 55",
                             "Diabetic: Glucose Level - 45", "Intercept"),
          ci =TRUE, star.cutoffs = c(0.05/8, 0.01/8, 0.001/8)) #Stars are corrected for multiple comparisons

modelList<-list()
for(i in 1:length(HORMONES)){
  modelList[[i]] <- lm(hormones[hormones$Day == 2 | (hormones$Day == 1 & hormones$Group == 1),HORMONES[i]] ~ Day + HYA + Day:HYA, data = hormones[hormones$Day == 2 | (hormones$Day == 1 & hormones$Group == 1),])}

stargazer(modelList[[1]], modelList[[2]], modelList[[3]], modelList[[4]],
          modelList[[5]], modelList[[6]], modelList[[7]], modelList[[8]],
          modelList[[9]], modelList[[10]], modelList[[11]], modelList[[12]],
          title = "Linear Models - Day 1 Controls vs. Diabetics",
          align= TRUE, type = "html", out = paste(FILEPATH, "LM_Day2ControlvsDiabetic_Hormones.htm"),
          omit.stat=c("LL","ser"),
          dep.var.labels = "Linear Regression Model",
          column.labels = as.character(HORMONES),
          covariate.labels=c("Diabetic", "Glucose Level - 65", "Glucose Level - 55", "Glucose Level - 45",
                             "Diabetic: Glucose Level - 65", "Diabetic: Glucose Level - 55",
                             "Diabetic: Glucose Level - 45", "Intercept"),
          ci =TRUE, star.cutoffs = c(0.05/8, 0.01/8, 0.001/8)) #Stars are corrected for multiple comparisons


################################3
#Now generating all the contrast level tests
#Group 0 = Diabetic, Group 1 = Controls
result_Group<-list()
result_Glucose<-list()
result_interaction<-list()

df$Group<-ifelse(df$Group == 0, "Diabetic", "Control")
df$Cluster<-paste(df$Group, df$Day, sep = "-")
for(i in MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX){
mod <- aov(df[df$Day == 1,i] ~ Group + Glucose.level + Group:Glucose.level, data = df[df$Day == 1,])
summary(mod)
result_Group[[i-(MINBRAINREGIONINDEX - 1)]]<-data.frame(ScheffeTest(mod)$'Group')
result_Glucose[[i-(MINBRAINREGIONINDEX - 1)]]<-data.frame(ScheffeTest(mod)$'Glucose.level')
result_interaction[[i-(MINBRAINREGIONINDEX - 1)]]<-data.frame(ScheffeTest(mod)$'Group:Glucose.level')}

CreateTable<-function(RESULT, COMPARISON){
result_Group_df<-data.frame(rbind(RESULT[[1]], RESULT[[2]], RESULT[[3]], RESULT[[4]],
                                  RESULT[[5]], RESULT[[6]], RESULT[[7]], RESULT[[8]]))
result_Group_df$Regions<-c(rep("Frontal Lobe", length(RESULT[[1]]$diff)), 
                           rep("Orbital Frontal", length(RESULT[[1]]$diff)), 
                           rep("Motor", length(RESULT[[1]]$diff)), 
                           rep("Parietal", length(RESULT[[1]]$diff)), 
                           rep("Temporal", length(RESULT[[1]]$diff)), 
                           rep("Occipital", length(RESULT[[1]]$diff)),
                           rep("Basal Ganglia", length(RESULT[[1]]$diff)), 
                           rep("Thalamus", length(RESULT[[1]]$diff)))
result_Group_df$Comparison<-rep(COMPARISON, length(result_Group_df$Regions))
result_Group_df<-result_Group_df[,c("Comparison", "Regions", "diff", "lwr.ci", "upr.ci", "pval")]
colnames(result_Group_df)[length(result_Group_df)]<-"Uncorrectedpvalue"
return(result_Group_df)}

Group_Comparison<-CreateTable(result_Group, "Control, day 1 - Diabetic")
Glucose_Comparison<-CreateTable(result_Glucose, "Glucose Comparison - Control D1 v Diabetic")
Interaction_Comparison<-CreateTable(result_interaction, "Group:Glucose Comparison")

for(i in MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX){
  mod <- aov(df[df$Group == "Control",i] ~ Day + Glucose.level + Day:Glucose.level, data = df[df$Group == "Control",])
  summary(mod)
  result_Group[[i-(MINBRAINREGIONINDEX - 1)]]<-data.frame(ScheffeTest(mod)$'Day')
  result_Glucose[[i-(MINBRAINREGIONINDEX - 1)]]<-data.frame(ScheffeTest(mod)$'Glucose.level')
  result_interaction[[i-(MINBRAINREGIONINDEX - 1)]]<-data.frame(ScheffeTest(mod)$'Day:Glucose.level')}

Group_Comparison<-rbind(Group_Comparison, CreateTable(result_Group, "Control, day 1 - Control, day 2"))
Glucose_Comparison<-rbind(Glucose_Comparison, CreateTable(result_Glucose, "Glucose.level, Day 1 vs Day 2"))
Interaction_Comparison<-rbind(Interaction_Comparison, CreateTable(result_interaction, "Control Day: Glucose Comparison"))


for(i in MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX){
  mod <- aov(df[(df$Group == "Control" & df$Day == 2) | (df$Group == "Diabetic"),i] ~ Cluster + Glucose.level + Cluster:Glucose.level, data = df[(df$Group == "Control" & df$Day == 2) | (df$Group == "Diabetic"),])
  summary(mod)
  result_Group[[i-(MINBRAINREGIONINDEX - 1)]]<-data.frame(ScheffeTest(mod)$'Cluster')
  result_Glucose[[i-(MINBRAINREGIONINDEX - 1)]]<-data.frame(ScheffeTest(mod)$'Glucose.level')
  result_interaction[[i-(MINBRAINREGIONINDEX - 1)]]<-data.frame(ScheffeTest(mod)$'Cluster:Glucose.level')}

Group_Comparison<-rbind(Group_Comparison, CreateTable(result_Group, "Diabetic, day 1 - Control, day 2"))
Glucose_Comparison<-rbind(Glucose_Comparison, CreateTable(result_Glucose, "Glucose.level, Diabetic vs Control Day 2"))
Interaction_Comparison<-rbind(Interaction_Comparison, CreateTable(result_interaction, "Diabetic vs Day 2 Control: Glucose Comparison"))

#Methods: Anova with post hoc Tukey Test
write.csv(Group_Comparison, paste(FILEPATH, "GroupComparisons.csv", sep = ""), row.names = FALSE)
write.csv(Glucose_Comparison, paste(FILEPATH, "GlucoseComparisons.csv", sep = ""), row.names = FALSE)
write.csv(Interaction_Comparison, paste(FILEPATH, "InteractionComparisons.csv", sep = ""), row.names = FALSE)


#################################################################################################################
#################################################################################################################
#################################################################################################################
#Looking at effect sizes in our comparisons
df_Diabetic<-df[df$Group == 0,]
df_Day1<-df[df$Group == 1 & df$Day == 1,]
df_Day2<-df[df$Day == 2,]
CalculateEffectSize<-function(mean1, mean2, sd1, sd2, n1, n2){
  numerator<-abs(mean1-mean2)
  denom<-sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1+n2-2))
  effSize<-numerator/denom
  return(effSize)
}

df_Diabetic$Cluster<-as.factor(paste(df_Diabetic$Day, df_Diabetic$Group, df_Diabetic$Glucose.level, sep = ", "))
mean.Diabetic<-aggregate(df_Diabetic[, MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX], list(df_Diabetic$Cluster), mean, na.rm = TRUE)
sd.Diabetic<-aggregate(df_Diabetic[, MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX], list(df_Diabetic$Cluster), sd, na.rm = TRUE)

df_Day1$Cluster<-as.factor(paste(df_Day1$Day, df_Day1$Group, df_Day1$Glucose.level, sep = ", "))
mean.df_Day1<-aggregate(df_Day1[, MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX], list(df_Day1$Cluster), mean, na.rm = TRUE)
sd.df_Day1<-aggregate(df_Day1[, MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX], list(df_Day1$Cluster), sd, na.rm = TRUE)

df_Day2$Cluster<-as.factor(paste(df_Day2$Day, df_Day2$Group, df_Day2$Glucose.level, sep = ", "))
mean.df_Day2<-aggregate(df_Day2[, MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX], list(df_Day2$Cluster), mean, na.rm = TRUE)
sd.df_Day2<-aggregate(df_Day2[, MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX], list(df_Day2$Cluster), sd, na.rm = TRUE)

EffectSizesforComparison<-function(mean.df1, mean.df2, sd.df1, sd.df2){
  EffectSizeList<-list()
  for(i in 2:length(mean.df1)){
    EffectSizeList[i-1]<-CalculateEffectSize(mean.df1[,i], mean.df2[,i], sd.df1[,i], sd.df2[,i], length(mean.df1), length(mean.df2))
  }
  EffectSizeList<-as.data.frame(EffectSizeList)
  names(EffectSizeList)<-names(mean.df1[2:length(mean.df1)])
  return(EffectSizeList)}

#Right now this collapses everything down
#need to separate out by blood sugar level
EffectSizeDay1toDay2<-EffectSizesforComparison(mean.df_Day1, mean.df_Day2, sd.df_Day1, sd.df_Day2)
EffectSizeDay1toDay2$Comparison<-"Day1toDay2"
EffectSizeDay1toDiabetic<-EffectSizesforComparison(mean.df_Day1, mean.Diabetic, sd.df_Day1, sd.Diabetic)
EffectSizeDay1toDiabetic$Comparison<-"Day1toDiabetic"
EffectSizeDiabetictoDay2<-EffectSizesforComparison(mean.Diabetic, mean.df_Day2, sd.Diabetic, sd.df_Day2)
EffectSizeDiabetictoDay2$Comparison<-"DiabetictoDay2"
EffectSizesWide<-rbind(EffectSizeDay1toDay2, EffectSizeDay1toDiabetic,EffectSizeDiabetictoDay2)
EffectSizes<-gather(EffectSizesWide, Comparison, EffectSize, Frontal:Thalamus, factor_key=TRUE)
colnames(EffectSizes)[1]<-"Region"
EffectSizes$Comparison<-rep(c("Day1vsDay2", "Day1vsDiabetic", "Day2vsDiabetic"), 8)

write.csv(EffectSizes, paste(FILEPATH, "EffectSizes_GroupDifferences.csv", sep = ""))
#################################################################################################################
#Comparing effect of changing glucose level relative to baseline


mean.Diabetic.Glucose<-aggregate(df_Diabetic[, MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX], list(df_Diabetic$Glucose.level), mean, na.rm = TRUE)
sd.Diabetic.Glucose<-aggregate(df_Diabetic[, MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX], list(df_Diabetic$Glucose.level), sd, na.rm = TRUE)

mean.Day1.Glucose<-aggregate(df_Day1[, MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX], list(df_Day1$Glucose.level), mean, na.rm = TRUE)
sd.Day1.Glucose<-aggregate(df_Day1[, MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX], list(df_Day1$Glucose.level), sd, na.rm = TRUE)

mean.Day2.Glucose<-aggregate(df_Day2[, MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX], list(df_Day2$Glucose.level), mean, na.rm = TRUE)
sd.Day2.Glucose<-aggregate(df_Day2[, MINBRAINREGIONINDEX:MAXBRAINREGIONINDEX], list(df_Day2$Glucose.level), sd, na.rm = TRUE)


GluLevEffSize<-function(mean, sd, meangroup, sdgroup, glucoselevel){
  EffectSizesforComparison(subset(mean, meangroup == 90), 
                           subset(mean, meangroup == glucoselevel), 
                           subset(sd, sdgroup == 90), 
                           subset(sd, sdgroup == glucoselevel))
}

EffectSizeDiabetics65<-GluLevEffSize(mean.Diabetic.Glucose, sd.Diabetic.Glucose, mean.Diabetic.Glucose$Group.1, sd.Diabetic.Glucose$Group.1, 65)
EffectSizeDiabetics65$Comparison<-"Diabetic 90 to 65"
EffectSizeDiabetics55<-GluLevEffSize(mean.Diabetic.Glucose, sd.Diabetic.Glucose, mean.Diabetic.Glucose$Group.1, sd.Diabetic.Glucose$Group.1, 55)
EffectSizeDiabetics55$Comparison<-"Diabetic 90 to 55"
EffectSizeDiabetics45<-GluLevEffSize(mean.Diabetic.Glucose, sd.Diabetic.Glucose, mean.Diabetic.Glucose$Group.1, sd.Diabetic.Glucose$Group.1, 45)
EffectSizeDiabetics45$Comparison<-"Diabetic 90 to 45"

EffectSizeDay165<-GluLevEffSize(mean.Day1.Glucose, sd.Day1.Glucose, mean.Day1.Glucose$Group.1, sd.Day1.Glucose$Group.1, 65)
EffectSizeDay165$Comparison<-"Day1 90 to 65"
EffectSizeDay155<-GluLevEffSize(mean.Day1.Glucose, sd.Day1.Glucose, mean.Day1.Glucose$Group.1, sd.Day1.Glucose$Group.1, 55)
EffectSizeDay155$Comparison<-"Day1 90 to 55"
EffectSizeDay145<-GluLevEffSize(mean.Day1.Glucose, sd.Day1.Glucose, mean.Day1.Glucose$Group.1, sd.Day1.Glucose$Group.1, 45)
EffectSizeDay145$Comparison<-"Day1 90 to 45"

EffectSizeDay265<-GluLevEffSize(mean.Day2.Glucose, sd.Day2.Glucose, mean.Day2.Glucose$Group.1, sd.Day2.Glucose$Group.1, 65)
EffectSizeDay265$Comparison<-"Day2 90 to 65"
EffectSizeDay255<-GluLevEffSize(mean.Day2.Glucose, sd.Day2.Glucose, mean.Day2.Glucose$Group.1, sd.Day2.Glucose$Group.1, 55)
EffectSizeDay255$Comparison<-"Day2 90 to 55"
EffectSizeDay245<-GluLevEffSize(mean.Day2.Glucose, sd.Day2.Glucose, mean.Day2.Glucose$Group.1, sd.Day2.Glucose$Group.1, 45)
EffectSizeDay245$Comparison<-"Day2 90 to 45"



EffectSizesWide<-rbind(EffectSizeDiabetics65, EffectSizeDiabetics55,EffectSizeDiabetics45,
                       EffectSizeDay165, EffectSizeDay155, EffectSizeDay145,
                       EffectSizeDay265, EffectSizeDay255, EffectSizeDay245)
EffectSizes<-gather(EffectSizesWide, Region, EffectSize, Frontal:Thalamus, factor_key=TRUE)
EffectSizes$Comparison<-rep(c("Diabetic 90 to 65", "Diabetic 90 to 55", "Diabetic 90 to 45",
                              "Day 1 90 to 65", "Day 1 90 to 55", "Day 1 90 to 45",
                              "Day 2 90 to 65", "Day 2 90 to 55", "Day 2 90 to 45"), 8)

write.csv(EffectSizes, paste(FILEPATH, "EffectSizes_BaselineComparison.csv", sep = ""))
#################################################################################################################
