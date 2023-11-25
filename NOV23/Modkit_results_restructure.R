library(readr)
library(dplyr)
library(tidyr)

All_MGMT <- read.csv("Documents/GitHub/MGMT_R/NOV23/SAMPLES.cpg.methyl.bed", sep ="\t", header = F)
Samples <- read.csv("Documents/GitHub/MGMT_R/NOV23/Samples.csv", sep =",", header = T)
Pyro <- read.csv("Documents/GitHub/MGMT_R/NOV23/Pyro_data.csv", header = T)

All_MGMT <- All_MGMT %>% select(2,10:19) %>% relocate(11)

names(All_MGMT) <- c("SampleID", "Pos", "Valid_cov", "Methylation_percent", "N_mod", "N_canon", "N_other_mod", "N_delete", "N_fail", "N_diff", "N_nocall")

All_MGMT_Island <- All_MGMT %>% filter(between(Pos, 129466683, 129467448))

cov <- All_MGMT_Island %>% group_by(SampleID) %>% summarise(AVG_cov = mean(Valid_cov), MED_cov = median(Valid_cov))

MGMT_long <- All_MGMT_Island %>% 
  select(c("SampleID", "Pos", "Methylation_percent")) %>%
  pivot_wider(names_from = Pos, values_from=Methylation_percent)

MGMT_RunSum <- inner_join(Samples,cov)

All_MGMT_Pyro <- All_MGMT %>% filter(between(Pos, 129467253, 129467272)) %>% select(1:4)

FOUR_CpGs <- inner_join(Pyro, All_MGMT_Pyro)

FOUR_CpGs$Pos <- as.factor(FOUR_CpGs$Pos)


