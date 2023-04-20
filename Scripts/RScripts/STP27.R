all_results_MGMT <- read.csv("R-markdown/Data/Aggregated_MGMT_24102022.csv")
MGMT_RunSum <- read.csv("R-markdown/Data/newRunsum.csv", row.names = 1)

# Fetch the positions for the two STP27 CpGs
STP27_CpGs <- all_results_MGMT %>% filter(Pos == 131265208 | Pos == 131265574 )

# make a column with the official epic array CpG IDs
STP27_CpGs$CpG_ID <- ifelse(STP27_CpGs$Pos == 131265208, "cg12434587", "cg12981137")

write.csv(STP27_CpGs, "R-markdown/Data/STP27_All_samples.csv")

library(dplyr)
library(tidyr)

STP27_CpGs <- read.csv("R-markdown/Data/STP27_All_samples.csv")

STP27_CpGs <- STP27_CpGs %>% select(c(DataSet,Sample, Percent_Methylated, CpG_ID)) %>%
                            spread(key=CpG_ID, value = Percent_Methylated)

names(STP27_CpGs) <- c("Series", "SampleID", "cg12434587","cg12981137")

STP27_CpGs <- left_join(STP27_CpGs, MGMT_RunSum)

STP27_CpGs <- na.omit(STP27_CpGs)

STP27_CpGs <- STP27_CpGs %>% filter(Series != "DenStem")

library(ggplot2)

# dotplot of the two CpGs

ggplot(STP27_CpGs, aes(x=cg12434587, y=cg12981137, color=Known_status)) + 
  geom_point(aes(size=Depth), alpha = 0.6)+
  scale_size_continuous(name="Depth",breaks = c(5,15,25,100,150), range = c(1,8))+
  scale_color_brewer(name="Known Status", palette = "Set1")+
 scale_y_sqrt(breaks=c(10,25,50,75,100))+
 scale_x_sqrt(breaks=c(10,25,50,75,100))+
  theme_bw()
  #facet_wrap(.~Series)


ggplot(STP27_CpGs, aes(x=cg12434587, y=cg12981137, color = Known_status)) + 
  geom_density_2d()+
  geom_point(aes(size=3), alpha = 0.6)+
  scale_size_continuous(name="Depth",breaks = c(5,15,25,100,150), range = c(1,8), guide="none")+
  scale_color_brewer(name="Known Status", palette = "Set1")+
  scale_y_sqrt(breaks=c(10,25,50,75,100))+
  scale_x_sqrt(breaks=c(10,25,50,75,100))+
  theme_bw(base_size = 18)
