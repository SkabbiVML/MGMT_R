pheatmap(Metsum_ROI_long,
         cluster_rows = T, 
         cluster_cols = F,
         clustering_method = "ward.D",
         cutree_rows = 2,
         treeheight_row = 20,
         border_color = NA,
         scale = "none",
         drop_levels = T,
         color = rev(brewer.pal(n = 10, name = "Spectral")),
         fontsize_col = 4,
         fontsize_row = 4,
         fontsize = 6,
         legend = T,
         show_rownames = F,
         legend_breaks = seq(0,100,10),
         annotation_col = mat_col,
         annotation_row = row_anno,
         annotation_colors = ann_colors,
         main = "All samples,  (n=145)",
         legend_labels = c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"))

#########################



ggplot(FOUR_CpGs, aes(x=Percent_Meth_Pyro, y=Percent_Meth_Nano, color = Pos))+
  geom_smooth(method = lm, color="black")+
  stat_regline_equation(label.y = 20, label.x = 50,aes(label = ..rr.label..),color = "black")+
  geom_point(aes(size=Depth), alpha = 0.6)+
  scale_size_continuous(name="Depth",breaks = c(5,15,25,100,150), range = c(1,8))+
  # geom_hline(yintercept = 10)+
  # geom_vline(xintercept = 10)+
  xlab("Pyro sequencing")+
  ylab("Nanopore sequencing")+
  theme_bw()+
  theme( aspect.ratio = 1)+
  scale_color_brewer(palette = "Set2", guide="none")+
  # scale_y_sqrt(breaks=c(10,25,50,75,100))+
  # scale_x_sqrt(breaks=c(10,25,50,75,100))+
  facet_wrap(.~Pos)

########################

GROUP <- FOUR_CpGs %>% 
          select(1:3,6:9) %>%
          group_by(Depth,Series, SampleID) %>%
          summarise(Average_Meth_Nano = mean(Percent_Meth_Nano),Average_Meth_Pyro = mean(Percent_Meth_Pyro)) 

ggplot(GROUP, aes(x=Average_Meth_Pyro, y=Average_Meth_Nano)) +
  geom_smooth(method = "lm", color="black") +
  geom_point(aes(size=Depth), alpha = 0.6, color="CornflowerBlue")+
  stat_regline_equation(label.y = 20, label.x = 60,aes(label = ..rr.label..),color = "black", size=4)+
  #scale_color_brewer(palette = "Set1", guide="none")+
  scale_size_continuous(name="Depth",breaks = c(5,10,50,100,150), range = c(1,16) )+
  xlab("Pyro sequencing")+
  ylab("Nanopore sequencing")+
  geom_hline(yintercept = 9)+
  geom_vline(xintercept = 9)+
  # scale_y_sqrt(breaks=c(10,25,50,75,100))+
  # scale_x_sqrt(breaks=c(10,25,50,75,100))+
  #facet_wrap(.~DataSet)+
  theme_bw()+
  theme(aspect.ratio = 1)

Compare_CpGs_Radium <- Compare_CpGs_GROUP %>% filter(DataSet == "Radium")

ggplot(Compare_CpGs_Radium, aes(x=Pyro_4_CpG_average, y=Nano_4_CpG_average)) +
  geom_smooth(method = "lm", color="black") +
  geom_point(aes(size=Coverage), alpha = 0.8, color="CornflowerBlue")+
  stat_regline_equation(label.y = 20, label.x = 60,aes(label = ..rr.label..),color = "black", size=4)+
  #scale_color_brewer(palette = "Set1", guide="none")+
  scale_size_continuous(name="Depth",breaks = c(5,10,50,100,150), range = c(1,10) )+
  xlab("Pyro sequencing")+
  ylab("Nanopore sequencing")+
  geom_hline(yintercept = 9)+
  geom_vline(xintercept = 9)+
  theme(aspect.ratio = 1)+
  theme_bw()

###################

Four_CpGs_group <- FOUR_CpGs %>% select(c(SampleID,Calls_num,Percent_Meth_Nano)) %>%
  group_by(SampleID) %>% summarise(Average_Meth_Nano = mean(Percent_Meth_Nano))

Four_CpGs_group <- na.omit(Four_CpGs_group)

p <- left_join(MGMT_RunSum,Four_CpGs_group)

p <- p %>% select(SampleID,Series,Depth, Known_status,Average_Meth_Nano)
p$Nano_Methylation_Status <- if_else(p$Average_Meth_Nano > 9, "Methylated", "UnMethylated")
p$Nano_Pyro_Concordance <- if_else(p$Known_status == p$Nano_Methylation_Status, "Concordant", "Discordant")

p2 <- p %>% select(SampleID, Series, Nano_Pyro_Concordance) %>%distinct()  %>% group_by(Series, Nano_Pyro_Concordance) %>% tally()

p <- p %>% filter(Series == "Radium")
########################
# confusion matrix
Nanopore <- factor(c("Methylated", "Methylated", "UnMethylated", "UnMethylated"))
PSQ <- factor(c("Methylated", "UnMethylated", "Methylated", "UnMethylated"))
Y      <- c(30, 7, 0, 31)
df <- data.frame(Nanopore, PSQ, Y)

library(ggplot2)
ggplot(data =  df, mapping = aes(x = PSQ, y = Nanopore)) +
  geom_tile(aes(fill = Y), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Y)), vjust = 1, size = 18) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_bw(base_size = 18) + theme(legend.position = "none")
################
#knitr::kable(p2)

ggplot(p, aes(x=Known_status, y=Average_Meth_Nano))+
  geom_hline(yintercept = 9)+
  geom_beeswarm(aes(size = 4, color = Known_status), alpha = 0.7, cex=3.5)+
  scale_size_continuous(guide="none")+
  scale_color_brewer(palette = "Set1", guide="none")+
  xlab("Known methylation status (PSQ)")+
  ylab("Nanopore methylation (%)")+
  #facet_grid(.~Series)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

```
```{r Concordance, echo=FALSE}
knitr::kable(p2)
