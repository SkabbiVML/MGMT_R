MGMT_RunSum <- read.csv("/home/vml/MGMT/MGMT_RunSum.csv", sep = ",", header = T)
MGMT_RunSum <- read.csv("/media/vml/Ghibli/MGMT/SampleMetadata/MGMT_RunSum.csv", sep = ",", header = T)

library(ggplot2)
library(ggbreak) 
library(patchwork)
library(dplyr)


MGMT_RunSum <- MGMT_RunSum %>% 
  mutate(Hits_ratio = 100*On_target_seqs/Read_number) %>% 
  filter(Comment != "Fail" & On_target_seqs > 1) %>% 
  group_by(Series) %>% 
  mutate(Series_median=median(On_target_seqs)) %>% 
  group_by(Comment) %>% 
  mutate(Comment_median=median(On_target_seqs)) %>%
  group_by(Series,Comment) %>%
  mutate(Series_Comment_median=median(On_target_seqs)) %>%
  group_by(Series,Comment) %>% 
  mutate(Series_Comment_mean=mean(On_target_seqs)) %>% 
  group_by(Pyro_Methylation_Status) %>% 
  mutate(Reads_per_status_median = median(On_target_seqs)) %>%
  group_by(Series,Pyro_Methylation_Status)%>%
  mutate(Reads_per_status_series_median = median(On_target_seqs))

MGMT_RunSum$Series_Comment_mean <- as.numeric(round(MGMT_RunSum$Series_Comment_mean,1))
MGMT_RunSum$Reads_per_status_series_median <- as.numeric(round(MGMT_RunSum$Reads_per_status_series_median,1))


# Distribution of hits in RIO by Methylation status and sequencing method
ggplot(MGMT_RunSum, aes(x=Comment, y=On_target_seqs))+
  geom_point(size = 5, position = position_jitterdodge(0.2,0), alpha = 0.6, aes(color=Pyro_Methylation_Status))+
  scale_color_viridis_d(name="Mehtylation Status")+
  scale_fill_discrete(guide="none")+
  stat_summary(fun=median, geom="point", shape=18, size=6, color="red",position = position_dodge(0.8), aes(fill=Pyro_Methylation_Status),name = "")+
  #ylim(0,125)+
  #scale_y_break(c(150, 250)) +
  scale_y_log10()+
  xlab ("")+
  ylab("Sequence depth") +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust=1))

# Distribution of hits in RIO by Series and sequencing method
ggplot(MGMT_RunSum, aes(x=Series, y=On_target_seqs))+
  geom_point(size = 5, position = position_jitterdodge(0.1,0), alpha = 0.6, aes(color=Comment)) +
  stat_summary(fun=median, geom="point", shape=18, size=6, color="red",position = position_dodge(0.8), aes(fill=Comment),name = "")+
  scale_color_viridis_d(name="Method")+
  scale_fill_discrete(guide="none")+
  labs(fill = "", color = "Method")+
  scale_color_viridis_d()+
  scale_y_log10()+
  xlab ("")+
  ylab("Sequence depth") +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust=1))

# mean hits on RIO per Series by method
ggplot(MGMT_RunSum, aes(x=Series, y=Series_Comment_median, fill=Comment))+
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_text(aes(label = Series_Comment_median), color="darkgrey", position = position_dodge(.9), vjust = -0.2)+
   labs(fill = "", color = "", y="Median on target")+
    scale_fill_viridis_d()+
  # ylim(0,100)+
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position="top")

# median hits on RIO per Series by Status
ggplot(MGMT_RunSum, aes(x=Series, y=Reads_per_status_series_median, fill=Pyro_Methylation_Status))+
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_text(aes(label = Reads_per_status_series_median), color="darkgrey", position = position_dodge(.9), vjust = -0.2)+
  labs(fill = "", color = "", y="Median on target")+
   scale_fill_viridis_d()+
 #  ylim(0,100)+
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position="top")

# Total reads versus reads in ROI
ggplot(MGMT_RunSum, aes(x=log(Read_number), y=log(On_target_seqs), color=Series))+
  geom_point(size = 5, alpha = 0.8)+
  #geom_smooth(method = lm)+
  scale_color_viridis_d()+
  theme_bw(base_size = 14) + theme(legend.position="top")

ggplot(MGMT_RunSum, aes(x=log(Read_number), y=log(On_target_seqs), color=Comment))+
  geom_point(size = 5, alpha = 0.8)+
 # geom_smooth(method = lm)+
  scale_color_viridis_d()+
  theme_bw(base_size = 14) + theme(legend.position="top")

# density plot of sequnce depth
ggplot(MGMT_RunSum, aes(x=On_target_seqs))+
  geom_density(aes(color=Comment),alpha = 0.5,position="identity")+
  scale_x_sqrt()+
  scale_color_viridis_d()+
  #xlim(0,200)+
  theme_bw()
