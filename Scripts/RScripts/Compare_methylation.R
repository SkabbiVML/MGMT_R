MGMT_RunSum <- read.csv("/media/vml/Ghibli/MGMT/SampleMetadata/MGMT_RunSum.csv", sep = ",", header = T)
all_results_MGMT <- read.csv("/media/vml/Ghibli/MGMT/Results/Aggregated_MGMT_08102022.csv", row.names = 1)
Radium_four <- read.csv("/media/vml/Ghibli/MGMT/SampleMetadata/Radium_FOUR_CpGs.csv", sep = ",", header = T)

library(dplyr)
library(tidyr)

FOUR_CpGs_MGMT <- all_results_MGMT %>% filter(Pos > 131265517) %>% filter(Pos < 131265537)

FOUR_CpGs_MGMT$Pos <- gsub("131265518", "CpG1", FOUR_CpGs_MGMT$Pos)
FOUR_CpGs_MGMT$Pos <- gsub("131265521", "CpG2", FOUR_CpGs_MGMT$Pos)
FOUR_CpGs_MGMT$Pos <- gsub("131265525", "CpG3", FOUR_CpGs_MGMT$Pos)
FOUR_CpGs_MGMT$Pos <- gsub("131265535", "CpG4", FOUR_CpGs_MGMT$Pos)
colnames(FOUR_CpGs_MGMT) <- gsub("Pos", "Nano", colnames(FOUR_CpGs_MGMT))
colnames(FOUR_CpGs_MGMT) <- gsub("Sample", "Sample_ID", colnames(FOUR_CpGs_MGMT))

FOUR_CpGs_MGMT <- FOUR_CpGs_MGMT %>% select(c(DataSet, Sample_ID, Nano, Percent_Methylated)) %>% 
  spread(Nano, Percent_Methylated)

colnames(FOUR_CpGs_MGMT) <- gsub("CpG", "Nano_CpG", colnames(FOUR_CpGs_MGMT))

Compare_CpGs <- inner_join(FOUR_CpGs_MGMT,Radium_four, By = "Sample_ID")

Compare_CpGs <- Compare_CpGs[!duplicated(Compare_CpGs), ]

Compare_CpGs <- Compare_CpGs[complete.cases(Compare_CpGs), ]

Compare_CpGs <- Compare_CpGs[,c(1,2,7,3:6,8:11)]

Compare_CpGs <- Compare_CpGs %>% gather(key = Pos, value = "Meth", -c(1:3)) %>% 
  separate(Pos, into = c("set", "Pos"), sep = "_") %>% 
  spread(set,Meth)

library(ggplot2)

# Scatter plot of all samples
ggplot(Compare_CpGs, aes(x=Nano, y=Pyro, color=Pos))+
  geom_point(size=3, alpha = 0.6)+
  theme_bw()+
  scale_color_viridis_d()

#######################################
Medium_coverage_samples <- MGMT_RunSum %>% filter(On_target_seqs > 9) %>% select(Sample_ID)

Medium_coverage_samples <- unique(Medium_coverage_samples$Sample_ID) 

High_coverage_samples <-  MGMT_RunSum %>% filter(On_target_seqs > 29) %>% select(Sample_ID)

High_coverage_samples <- unique(High_coverage_samples$Sample_ID)
###########################

Compare_CpGs_HIGH <- Compare_CpGs %>% filter(Sample_ID %in% High_coverage_samples)

Compare_CpGs_Medium <- Compare_CpGs %>% filter(Sample_ID %in% Medium_coverage_samples)


# Scatter plot of only high coverage
ggplot(Compare_CpGs, aes(x=Pyro, y=Nano, color = Pos))+
  geom_smooth(method = lm, color="black")+
  stat_regline_equation(label.y = 2, label.x = 6,aes(label = ..rr.label..),color = "black")+
  geom_point(aes(size=On_target_seqs), alpha = 0.6)+
  scale_size_continuous(name="Depth",breaks = c(5,15,25,100,150))+
  geom_hline(yintercept = 10)+
  geom_vline(xintercept = 10)+
  theme_bw()+
    scale_color_brewer(palette = "Paired")+
  scale_y_sqrt()+
  scale_x_sqrt()+
  facet_wrap(.~Pos)

#   OR
library(ggpubr)
P1 <- ggplot(Compare_CpGs_HIGH, aes(x=Pyro, y=Nano, color = Pos))
P2 <- P1 + geom_smooth(method = "lm", color="black")
P2 + geom_point(size=3, alpha = 0.6)+
  stat_regline_equation(label.y = 25, label.x = 75,aes(label = ..rr.label..),color = "black")+
  theme_bw()+
  scale_color_brewer(palette = "Paired")+
  ylim(0,100)+
  xlim(0,100)



######################################################
######### CpG group average ########################3
################################################3
Denstem <- MGMT_RunSum %>% filter(Series == "DenStem") %>% select(c("Series", "Sample_ID", "Pyro_Methylation_Status", "Pyro_Methylation_percent"))

Denstem <- Denstem[!duplicated(Denstem),]

#######
Compare_CpGs_GROUP <- FOUR_CpGs_MGMT %>% 
  gather(key = "CpG", value = "Nano_Meth", -c(1,2)) %>%
  group_by(Sample_ID,DataSet) %>%
  summarise("Nano_4_CpG_average" = mean(Nano_Meth)) %>%
  filter(DataSet != "Areeba")

Compare_CpGs_GROUP$DataSet <- gsub("VML", "DenStem", Compare_CpGs_GROUP$DataSet)
Compare_CpGs_GROUP$DataSet <- gsub("Radium1", "Radium", Compare_CpGs_GROUP$DataSet)
Compare_CpGs_GROUP$DataSet <- gsub("Radium2", "Radium", Compare_CpGs_GROUP$DataSet)


Radium_Pyro_GROUP <- Radium_four %>% 
  gather(key = "CpG", value = "Pyro_Meth", -c(1,2)) %>%
  group_by(Sample_ID) %>%
  summarise("Pyro_4_CpG_average" = mean(Pyro_Meth))

Radium_Pyro_GROUP$DataSet <- "Radium"

Denstem_Pyro_GROUP <- MGMT_RunSum %>% filter(Series == "DenStem") %>%
  select(c("Sample_ID", "Series", "Pyro_Methylation_percent"))

colnames(Denstem_Pyro_GROUP) <- c("Sample_ID", "DataSet", "Pyro_4_CpG_average")

Pyro_GROUP <- rbind(Denstem_Pyro_GROUP,Radium_Pyro_GROUP)
Pyro_GROUP <- Pyro_GROUP[!duplicated(Pyro_GROUP),]

Compare_CpGs_GROUP <- left_join(Compare_CpGs_GROUP, Pyro_GROUP, by = c("Sample_ID","DataSet") )

Compare_CpGs_GROUP <- Compare_CpGs_GROUP[complete.cases(Compare_CpGs_GROUP), ]

###### Select High and medium coverage samples
Compare_CpGs_GROUP_HIGH <- Compare_CpGs_GROUP %>% filter(Sample_ID %in% High_coverage_samples)

Compare_CpGs_GROUP_Medium <- Compare_CpGs_GROUP %>% filter(Sample_ID %in% Medium_coverage_samples)
#################

############# PLOT

ggplot(Compare_CpGs_GROUP, aes(x=Pyro_4_CpG_average, y=Nano_4_CpG_average,color=DataSet))+
  geom_point(alpha = 0.6, aes( size=Coverage))+
    geom_smooth(method = lm)+
  stat_regline_equation(label.y = 2, label.x = 6,aes(label = ..rr.label..),color = "black")+
  scale_color_brewer(palette = "Paired")+
  xlab("Pyro")+
  ylab("Nano")+
  geom_hline(yintercept = 9)+
  geom_vline(xintercept = 9)+
  scale_y_sqrt()+
  scale_x_sqrt()+
  #facet_wrap(.~DataSet)+
  theme_bw()

### OR

P3 <- ggplot(Compare_CpGs_GROUP, aes(x=Pyro_4_CpG_average, y=Nano_4_CpG_average, color = DataSet))
P4 <- P3 + geom_smooth(method = "lm", color="black")
P4 + geom_point(aes(size=Coverage), alpha = 0.6)+
  stat_regline_equation(label.y = 2, label.x = 6,aes(label = ..rr.label..),color = "black", size=8)+
  scale_color_brewer(palette = "Paired")+
  scale_size_continuous(name="Depth",breaks = c(10,100,850), range = c(1,8) )+
  xlab("Pyro")+
  ylab("Nano")+
  geom_hline(yintercept = 9)+
  geom_vline(xintercept = 9)+
  scale_y_sqrt()+
  scale_x_sqrt()+
  facet_wrap(.~DataSet)+
  theme_bw()


