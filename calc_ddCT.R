#Load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

#Load files (spliced and retained are arg[1:2] output is arg[3])
file_input = commandArgs(trailingOnly = T)

#Function to determine which is spliced vs retained
if(str_detect(file_input[1], "ss")==TRUE & str_detect(file_input[2], "iu")==TRUE){
  spliced = fread(file_input[1])
  retained = fread(file_input[2])
}

if(str_detect(file_input[1], "iu")==TRUE & str_detect(file_input[2], "ss")==TRUE){
  retained = fread(file_input[1])
  spliced = fread(file_input[2])
}

#Split Name column into the ID data
spliced = spliced %>% separate(Name, c("gene", "isoform", "time", "rep"), " ")
retained = retained %>% separate(Name, c("gene", "isoform", "time", "rep"), " ")

#Get average of spliced for each timepoint so can do deltaCT
mean_spliced = spliced %>% group_by(time) %>% summarize(meanCT_spliced = mean(Ct))

#Do delta CT against spliced
spliced = left_join(spliced, mean_spliced, by="time") %>% mutate(dCT = Ct-meanCT_spliced)
retained = left_join(retained, mean_spliced, by="time") %>% mutate(dCT = Ct-meanCT_spliced)

#Now need to do ddCT for each isoform against timepoint 0 of its condition
dCT0_spliced = spliced %>% 
  group_by(time) %>% 
  summarize(average_dCT = mean(dCT)) %>% 
  filter(time=="t0") %>% 
  select(average_dCT) %>% 
  unlist() %>% unname()

dCT0_retained = retained %>% 
  group_by(time) %>% 
  summarize(average_dCT = mean(dCT)) %>% 
  filter(time=="t0") %>% 
  select(average_dCT) %>% 
  unlist() %>% unname()

spliced = spliced %>% add_column(dCT0_spliced=dCT0_spliced)
retained = retained %>% add_column(dCT0_retained=dCT0_retained)

#ddCT calculation
spliced = spliced %>% mutate(ddCT = dCT - dCT0_spliced) %>% select(-c(meanCT_spliced, dCT0_spliced))
retained = retained %>% mutate(ddCT = dCT - dCT0_retained) %>% select(-c(meanCT_spliced, dCT0_retained))

#2^-ddCT calculation
spliced = spliced %>% mutate(`2^-ddCT`=2^-ddCT)
retained = retained %>% mutate(`2^-ddCT`=2^-ddCT)

#Output
output = rbind(spliced, retained)

#Write out
fwrite(output, file_input[3])