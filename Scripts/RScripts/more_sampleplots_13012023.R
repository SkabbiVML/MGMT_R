library(tidyr)
library(dplyr)
library(ggplot2)

FOUR_CpGs <- read.csv("Data/FOUR_CpGs_MGMT_AllDataSets_04012023.csv")

Four_CpGs_group <- FOUR_CpGs %>% select(c(Sample_ID,Calls_num,Percent_Methylated)) %>%
  group_by(Sample_ID) %>% summarise(Average_Methylation = mean(Percent_Methylated))

Four_CpGs_group <- na.omit(Four_CpGs_group)

p <- left_join(MGMT_RunSum,Four_CpGs_group)

ggplot(p, aes(x=Pyro_Methylation_Status, y=Average_Methylation, color = Series))+
  geom_hline(yintercept = 9)+
  geom_point(aes(size = On_target_seqs), alpha = 0.6, position = "jitter")+
    scale_size_continuous(name="Depth",breaks = c(5,10,50,100), range = c(1,10) )+
  scale_color_brewer(palette = "Set1", guide="none")+
  xlab("Classification by pyrosequencing")+
  ylab("Average methylation per nanopore")+
  facet_grid(.~Series)+
  theme_bw()

# Adding the cumulative methylation numbers (MGMT_Histogram.R) to p
p <- left_join(p,t3)

#### can filter to only show GBMs

p <- p %>% filter(Diagnosis == "Glioblastoma")

# Average methylation of 4 CpGs (nanopore) versus cumulative methylation

ggplot(p, aes(x=CumMeth, y=Average_Methylation, color=Pyro_Methylation_Status))+
  geom_point(aes(size = On_target_seqs), alpha = 0.6, position = "jitter")+
  scale_size_continuous(name="Depth",breaks = c(5,10,50,100), range = c(1,10) )+
  scale_color_brewer(palette = "Set1", guide="none")+
  xlab("Cumulative methylated CpGs (max 98)")+
  ylab("Average methylation % of 4CpGs (nanopore)")+
  theme_bw()+
  facet_grid(.~Series)


# Average methylation of 4 CpGs (pyrosequencing) versus cumulative methylation

ggplot(p, aes(x=CumMeth, y=Pyro_Methylation_percent, color=Pyro_Methylation_Status))+
  geom_point(aes(size = On_target_seqs), alpha = 0.6, position = "jitter")+
  scale_size_continuous(name="Depth",breaks = c(5,10,50,100), range = c(1,10) )+
  scale_color_brewer(palette = "Set1", guide="none")+
  xlab("Cumulative methylated CpGs (max 98)")+
  ylab("Average methylation % of 4CpGs (pyro)")+
  theme_bw()+
  facet_grid(.~Series)

# Beeswarm plots

# Average Nano 4CpGs
ggplot(p, aes(x=Pyro_Methylation_Status, y=Average_Methylation, color = Series))+
  geom_hline(yintercept = 9)+
  geom_beeswarm(aes(size = On_target_seqs), alpha = 0.6, cex=3)+
  scale_size_continuous(name="Depth",breaks = c(5,10,50,100), range = c(1,10) )+
  scale_color_brewer(palette = "Set1", guide="none")+
  xlab("Known classification")+
  ylab("4 CpGs, average methylation by nanopore")+
  facet_grid(.~Series)+
  theme_bw()

# Cumulative number of methylated CpGs
ggplot(p, aes(x=Pyro_Methylation_Status, y=CumMeth, color = Series))+
  #geom_hline(yintercept = 9)+
  geom_beeswarm(aes(size = On_target_seqs), alpha = 0.6, cex=3)+
  scale_size_continuous(name="Depth",breaks = c(5,10,50,100), range = c(1,10) )+
  scale_color_brewer(palette = "Set1", guide="none")+
  xlab("Classification by pyrosequencing")+
  ylab("Cumulative methylated CpGs (max 98)")+
  facet_grid(.~Series)+
  theme_bw()