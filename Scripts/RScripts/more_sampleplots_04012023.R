library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(ggbeeswarm)


FOUR_CpGs <- read.csv("Data/FOUR_CpGs_MGMT_AllDataSets_04012023.csv")

Four_CpGs_group <- FOUR_CpGs %>% select(c(Sample_ID,Calls_num,Percent_Methylated)) %>%
  group_by(Sample_ID) %>% summarise(Average_Methylation = mean(Percent_Methylated))

Four_CpGs_group <- na.omit(Four_CpGs_group)

p <- left_join(MGMT_RunSum,Four_CpGs_group)

p <- p %>% select(Sample_ID,Series,On_target_seqs, Pyro_Methylation_Status,Average_Methylation)
p$Nano_Methylation_Status <- if_else(p$Average_Methylation > 9, "Methylated", "UnMethylated")
p$Nano_Pyro_Concordance <- if_else(p$Pyro_Methylation_Status == p$Nano_Methylation_Status, "Concordant", "Discordant")

p2 <- p %>% select(Sample_ID, Series, Nano_Pyro_Concordance) %>%distinct()  %>% group_by(Series, Nano_Pyro_Concordance) %>% tally()

knitr::kable(p2)

ggplot(p, aes(x=Pyro_Methylation_Status, y=Average_Methylation, color = Series))+
  geom_hline(yintercept = 9)+
  geom_beeswarm(aes(size = On_target_seqs), alpha = 0.6, cex=3)+
  scale_size_continuous(name="Depth",breaks = c(5,10,50,100), range = c(1,10) )+
  scale_color_brewer(palette = "Set1", guide="none")+
  xlab("Classification by pyrosequencing")+
  ylab("Average methylation per nanopore")+
  facet_grid(.~Series)+
  theme_bw()


