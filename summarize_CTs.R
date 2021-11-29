
library(data.table)
library(dplyr)

#Take file input on command line to allow integration into pipelines (Rscript summarize_CTs.R 'file.csv' 'output_file.csv')
file_input = commandArgs(trailingOnly = T)

#Read as a data.table
file = fread(file_input[1])

#Summarize, get mean, sd and N, and use this to make sem
summary = file %>% 
  group_by(Name) %>% 
  summarize(Ct_mean=mean(Ct), sd = sd(Ct), n=n()) %>% 
  mutate(Ct_sem = sd/sqrt(n)) %>%
  select(-c(sd, n))

#Change col names back to Ct and SEM
colnames(summary)=c("Name", "Ct", "SEM")

#Write out
fwrite(summary, file_input[2])


