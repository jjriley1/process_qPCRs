#Load libraries
library(dplyr)
library(data.table)

#Take file input on command line to allow integration into pipelines (Rscript remove_outlying_techreps.R 'file.csv')
file_input = commandArgs(trailingOnly = T)

#Read as a data.table
file = fread(file_input)

#Group by name and get mean and sd
file = file %>% group_by(Name) %>% summarize(.groups="keep", Name = Name, Ct = Ct, Mean_Ct = mean(Ct), SD_Ct = sd(Ct))

#Compared to the rest of the data, which sets of technical replicates have abnormally high SDs? 
SDs = file %>% ungroup %>% select(SD_Ct) %>% unname() %>% unique() %>% unlist()
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
test[6]=""
colnames(test)[6]="drop"
test = transform(test, drop=as.logical(drop))
for(x in unique(test$Name)){
  values = test %>% filter(Name==x) %>% ungroup() %>% select(Ct) %>% unlist() 
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
  
  drop = test %>% filter(Name==x) 
  drop[6] = output
  colnames(drop)[6] = "drop"
  test = full_join(test, drop, by=c("Name", "Ct", "Mean_Ct", "SD_Ct", "high_SD", "drop"))
  
}

#this duplicates them but with NA, remove duplicates
test = test %>% filter(!is.na(drop))


#Drop those where drop == TRUE 
#join test (high SD) and file (total)
file = left_join(file, test, by=c("Name", "Ct", "Mean_Ct", "SD_Ct", "high_SD"))

#set those NAs to FALSES in drop column
file[6][is.na(file[6])]=F

#only keep those where drop == FALSE
file = file %>% filter(drop==FALSE)

#Remove unnecessary columns to allow output
file = file[,1:2]

#Output as data.table
dir.create("outliers_removed")
fwrite(file, paste0("outliers_removed/no_outliers_",file_input))
print("Final output:")
print(file)

#Output log of those removed due to being outlier (also print the rest of the group so it can be verified)
dir.create("outliers_removed/dropped")
fwrite(test, paste0("outliers_removed/dropped/dropped_", file_input))
print("Dropped from final output: (where drop=TRUE)")
print(test)