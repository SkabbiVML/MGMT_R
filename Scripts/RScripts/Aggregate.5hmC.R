library(dplyr)
library(tidyr)
library(readr)

# Import the EPIC sites
Five.mC.ROI.Radium2 <- read.table("MGMT/Analysis/Results/5hmC/Radium2/5mC.ROI.Radium2.bed", sep = "\t")

Five.mC.ROI.Radium2 <- separate(Five.mC.ROI.Radium2, "V1", into = c("Sample_ID", "Chromosome"), sep = ", ")

Five.mC.ROI.Radium2 <- separate(Five.mC.ROI.Radium2, "Sample_ID", into = c("Sample_ID", "Series"), sep = "/")

Five.mC.ROI.Radium1$Series <- "Radium"

Five.mC.ROI.Radium2 <- Five.mC.ROI.Radium2 %>% select(1,2,4,6,7,12:16)

names(Five.mC.ROI.Radium2) <- c("Sample_ID", "Series", "Pos", "Type", "Quality", "Depth", "Percent_Methylated", "UnMethylated", "Methylated", "No_call")

ALL.samples.5mc <- rbind(Five.mC.ROI.Areeba, Five.mC.ROI.Radium1,Five.mC.ROI.Radium2)

ALL.samples.5hmc <- rbind(Five.hmC.ROI.Areeba,Five.hmC.ROI.Radium1,Five.hmC.ROI.Radium2)

write.csv(ALL.samples.5hmc, "Documents/GitHub/MGMT_R/Data/5hmC.full.csv", row.names = F)
write.csv(ALL.samples.5mc, "Documents/GitHub/MGMT_R/Data/5mC.full.csv", row.names = F)

Metsym.ROI.5hmC <- ALL.samples.5hmc %>% select(Sample_ID,Pos,Percent_Methylated) %>% spread(Pos,Percent_Methylated) %>% na.omit()
Metsym.ROI.5mC <- ALL.samples.5mc %>% select(Sample_ID,Pos,Percent_Methylated) %>% spread(Pos,Percent_Methylated) %>% na.omit()

write.csv(Metsym.ROI.5hmC, "Documents/GitHub/MGMT_R/Data/Metsum.ROI.5hmC.csv", row.names = F)
write.csv(Metsym.ROI.5mC, "Documents/GitHub/MGMT_R/Data/Metsum.ROI.5mC.csv", row.names = F)
