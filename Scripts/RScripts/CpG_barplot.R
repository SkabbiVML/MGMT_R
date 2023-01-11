setwd("~/R/cas9R")
setwd("~/MGMT/Methylation_call_files")

Patient <- "T21-216"

Methylation_summary <- read.delim("Summarized_methylation_T21-216.tsv", header = T)

# How many replicates you want of each row. How many CpG sites are in each motif
duptimes <- Methylation_summary$num_motifs_in_group

# Create an index of the rows you want with duplications
idx <- rep(1:nrow(Methylation_summary), duptimes)

#Use that index to genderate your new data frame
MetSum <- Methylation_summary[idx,]


library(dplyr)

CpG_window <- MetSum %>% filter(start > 131264926 & start < 131265795)

CpG_window$CpG <- c(1:100)

Big_window <- MetSum %>% filter(start > 131264926 & start < 131271868)

Big_window$CpG <- c(1:nrow(Big_window))



library(ggplot2)
library(scales)

ggplot(CpG_window, aes(x=CpG,y=methylated_frequency)) + 
  geom_bar(stat = "identity", fill = "red", color = "blue", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format())+
  geom_vline(xintercept = c(75.5,79.5), linetype = "dashed" , colour = "red", size = 1) +
  ggtitle(Patient)+
  theme_bw()


ggplot(Big_window, aes(x=CpG,y=methylated_frequency)) + 
  geom_bar(stat = "identity", fill = "blue", color = "grey") +
  scale_y_continuous(labels = scales::percent_format())+
  geom_vline(xintercept = c(75.5,83.5), linetype = "dashed" , colour = "red", size = 1) +
  ggtitle(Patient)+
  theme_bw()
