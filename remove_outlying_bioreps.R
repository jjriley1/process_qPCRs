#Load libraries
library(dplyr)
library(data.table)
library(tidyr)
library(tibble)

#Take file input on command line to allow integration into pipelines (Rscript remove_outlying_techreps.R 'file.csv')
file_input = commandArgs(trailingOnly = T)
#Read as a data.table
file = fread(file_input[1])
file = file %>% separate(Name, c("gene", "isoform", "time", "rep"), " ")
#Group by name and get mean and sd
file_summary = file %>% group_by(time) %>% summarize( Mean_Ct = mean(Ct), SD_Ct = sd(Ct))
file = left_join(file, file_summary)
#Compared to the rest of the data, which sets of technical replicates have abnormally high SDs? 
SDs = file %>% ungroup %>% select(SD_Ct) %>% distinct() %>% unname() %>% unlist()
SDs_mean = mean(SDs)
SDs_SD = sd(SDs)
n = length(SDs)
error = qnorm(0.975)*SDs_SD/sqrt(n)
cutoff = SDs_mean+error

#failsafe in case error is very low, we don't want to exclude decent samples
if(error < 0.2){
  error = 0.2
}

#Add column with bool to distinguish groups containing outlier vs non-outliers
file = file %>% mutate(high_SD = SD_Ct>cutoff)

#To spot the outlier in the triplet we will compare the each data point to the other two, and omit the point which has a difference > cutoff in both comparisons. If all comparisons > cutoff then the entire sample will be omitted. 
test = file %>% filter(high_SD == TRUE)
test = test %>% add_column(drop=NA)
test = transform(test, drop=as.logical(drop))
for(x in unique(test$time)){
  values = test %>% filter(time==x) %>% ungroup() %>% select(Ct) %>% unlist() 
  print(values)
  comp1 = abs(values[1] - values[2])
  comp2 = abs(values[1] - values[3])
  comp3 = abs(values[2] - values[3])
  
  #outputs (F=keep, T=drop)
  #if 1-2 is < cutoff, then 3 caused the high sd
  if(comp1<cutoff){
    output = c(F,F,T)
  }
  #if 1-3 < cutoff, then 2 caused the high sd
  if(comp2<cutoff){
    output = c(F,T,F)
  }
  # if 2-3 < cutoff, then 1 caused the high sd
  if(comp3<cutoff){
    output = c(T,F,F) 
  }
  #if all greater than cutoff, then drop all
  if(comp1>=cutoff & comp2>=cutoff & comp3>=cutoff){
    output = c(T,T,T)
  }
  
  drop = test %>% filter(time==x) 
  drop[,10] = output
  test = full_join(test, drop)
  
}

#this duplicates them but with NA, remove duplicates
test = test %>% filter(!is.na(drop))


#Drop those where drop == TRUE 
#join test (high SD) and file (total)
file = left_join(file, test)

#set those NAs to FALSES in drop column
file[,10][is.na(file[,10])]=F

#only keep those where drop == FALSE
file = file %>% filter(drop==FALSE)

#Remove unnecessary columns to allow output
file = file[,1:6]
file = file %>% unite("Name", c(gene, isoform, time, rep), sep=" ")

#Output as data.table
fwrite(file, file_input[2])
print("Final output:")
print(file)

print("Dropped from final output: (where drop=TRUE)")
print(test)
