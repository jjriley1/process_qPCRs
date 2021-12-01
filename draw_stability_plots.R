#Load libraries
library(dplyr)
library(ggplot2)
library(data.table)
library(stringr)

#Take file input on command line to allow integration into pipelines (Rscript draw_graphs.R 'file.csv' 'outputs.png')
file_input = commandArgs(trailingOnly = T)

#Read as a data.table
file = fread(file_input[1])

#Force time to be numeric 
file = file %>% mutate(time=as.numeric(str_remove_all(time,"t")))

#Summarize and get mean +/- SEM
file_summary = file %>% group_by(isoform, time) %>% summarize(mean=mean(`2^-ddCT`), sd = sd(`2^-ddCT`), n=n()) %>% mutate(sem=sd/sqrt(n))

#Get gene name
gene = file[1,1] %>% unlist() %>% unname() %>% toupper()

#Plot with ggplot2
file_summary %>% ggplot(aes(x=time, y=mean, group=isoform, color=isoform)) + geom_pointrange(aes(ymin=mean-sem, ymax=mean+sem)) + geom_line() + 
  labs(title=paste0(gene, " stability"),
       x = "Time (hours)",
       y = expression("2"^("-"*Delta*Delta*"CT")))

#Save as GENENAME_stability_plot.png
ggsave(file_input[2], width=8, height=6)