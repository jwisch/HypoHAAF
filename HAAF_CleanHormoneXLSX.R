library(readxl)    

FILEPATH<-"C:/Users/julie.wisch/Documents/Hypo/"

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

df <- read_excel_allsheets("HormonesforASL.xlsx")


CleanUpHormones<-function(EPI){
EPI<-EPI[!is.na(EPI$HYA),]
colnames(EPI)[2]<-"labeled"
index <- which(EPI$HYA=="T1DM")

df1<-data.frame(EPI[1:(index-1),])
df1<-df1[, -grep("X_", colnames(df1))]
df1$labeled<-as.numeric(df1$labeled)
df2<-data.frame(EPI[index+1:length(EPI),])
colnames(df2)<-EPI[index,]
colnames(df2)[2]<-"labeled"
df2$labeled<-as.numeric(df2$labeled)

df2$Day<-"Day 1"
df1$Day<-c(rep("Day 1", 12), rep("Day 2", 12))
EPI<-merge(df1, df2, by.x = c("HYA", "labeled", "Day"), by.y = c("T1DM", "labeled", "Day"), all = TRUE)
EPI<-EPI[with(EPI, order(labeled, Day)),]
EPI<-EPI[!is.na(EPI$labeled),]
EPI[4:length(EPI)] <- sapply(EPI[4:length(EPI)],as.numeric)


EPI<-EPI %>%
  group_by(HYA, Day) %>% 
  summarise_each(funs(mean(., na.rm = TRUE)))
EPI<-data.frame(EPI)
return(EPI)}

df[[1]]<-NULL 

HORMONES<-c("EPI", "Norepi", "ffa", "glucose", "C-pep",
            "glucagon", "cortisol", "lactate", "insulin",
            "TSS", "NGP", "NG")



hold<-CleanUpHormones(df[[1]])
hold$Hormone<-HORMONES[1]

for(i in 2:length(HORMONES)){
  hold2<-CleanUpHormones(df[[i]])
  hold2$Hormone<-HORMONES[i] 
  hold<-data.frame(bind_rows(hold, hold2))
}

#Dealing with the fact that the column names are all jacked up
CombineColumns<-function(COLUMN, othercolumn = NA, othercolumn2 = NA, othercolumn3 = NA){
hold[,COLUMN]<-paste(hold[,othercolumn], hold[,othercolumn2], hold[,othercolumn3], hold[,COLUMN])
hold[,COLUMN]<-gsub("NaN", "", hold[,COLUMN])
hold[,COLUMN]<-gsub("NA", "", hold[,COLUMN])
hold[,COLUMN]<-as.numeric(hold[,COLUMN])
return(hold[,COLUMN])}

hold$BOLD.59<-CombineColumns("BOLD.59", "BOLD.59.B", "BOLD.59.V1", "BOLD.59B")

CombineColumns<-function(COLUMN, othercolumn = NA){
  hold[,COLUMN]<-paste(hold[,COLUMN], hold[,othercolumn])
  hold[,COLUMN]<-gsub("NaN", "", hold[,COLUMN])
  hold[,COLUMN]<-gsub("NA", "", hold[,COLUMN])
  hold[,COLUMN]<-as.numeric(hold[,COLUMN])
  return(hold[,COLUMN])}

hold$BOLD.36<-CombineColumns("BOLD.36", "BOLD.36__1")
hold$BOLD.36<-CombineColumns("BOLD.36", "BOLF.36")
hold$BOLD.61<-CombineColumns("BOLD.61", "BOLD.61.1")
hold$BOLD.08<-CombineColumns("BOLD.08", "BOLD.08__1")
hold$BOLD.56<-CombineColumns("BOLD.56", "bold.56")
hold$BOLD.14<-CombineColumns("BOLD.14", "bold.14")
hold$BOLD.09<-CombineColumns("BOLD.09", "BOLD.9")
hold$BOLD.09<-CombineColumns("BOLD.09", "BOLD.09.1.5.16")

hold<-hold[ , !names(hold) %in% c("BOLD.59.B", "BOLD.59.V1", "BOLD.59B", "BOLD.36__1", "BOLF.36",
                                  "BOLD.61.1", "BOLD.08__1", "bold.56", "bold.14", "BOLD.9", "BOLD.09.1.5.16", "NA.")] 

hold<-hold[,c(1:3, 52, 4:51, 53:length(hold))]

groups<-read.csv(paste(FILEPATH, "HypoGlycemiaGroups.csv", sep = ""))
groups$X<-gsub("-", ".", groups$X)


write.csv(hold, paste(FILEPATH, "Hormones_Updated.csv", sep = ""), row.names = FALSE)
