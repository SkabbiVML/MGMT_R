library(readr)
library(dplyr)
library(tidyr)

# Load the relevant data necessary to complete the analysis. Should be packaged into Rdata

# Contains all methylation sites in all samples
All_MGMT <- read.csv("GitHub/MGMT_R/NOV23/SAMPLES.cpg.methyl.bed", sep ="\t", header = F)
# Sample info
Samples <- read.csv("GitHub/MGMT_R/NOV23/Samples.csv", sep =",", header = T)
# Pyrosequencing results for the Radium samples
Pyro <- read.csv("GitHub/MGMT_R/NOV23/Pyro_data.csv", header = T)

# select the relevant columns
All_MGMT <- All_MGMT %>% select(2,10:19) %>% relocate(11)

# name the columns
names(All_MGMT) <- c("SampleID", "Pos", "Valid_cov", "Methylation_percent", "N_mod", "N_canon", "N_other_mod", "N_delete", "N_fail", "N_diff", "N_nocall")

# Extract the CpGs in the CpG island (98)
All_MGMT_Island <- All_MGMT %>% filter(between(Pos, 129466683, 129467448))

#or the extended CpG island (115)
All_MGMT_ext_Island <- All_MGMT %>% filter(between(Pos, 129466250, 129467831))

# Calculate mean and median coverage 
cov <- All_MGMT_Island %>% group_by(SampleID) %>% summarise(AVG_cov = mean(Valid_cov), MED_cov = median(Valid_cov))

# Re-structure the data to a table where each sample is a row and every position is a column
MGMT_long_island <- All_MGMT_Island %>% 
  select(c("SampleID", "Pos", "Methylation_percent")) %>%
  pivot_wider(names_from = Pos, values_from=Methylation_percent) %>%
  column_to_rownames("SampleID")

MGMT_long_ext_island <- All_MGMT_ext_Island %>% 
  select(c("SampleID", "Pos", "Methylation_percent")) %>%
  pivot_wider(names_from = Pos, values_from=Methylation_percent) %>%
  column_to_rownames("SampleID")


# Add the CpG methylation coverage to the sample info table
MGMT_RunSum <- inner_join(Samples,cov)

# Extract the CpG sites that correspond to the MGMT Pyro kit sites
All_MGMT_Pyro <- All_MGMT %>% filter(between(Pos, 129467253, 129467272)) %>% select(1:4)

# Add the Pyro kit results to the table
FOUR_CpGs <- inner_join(Pyro, All_MGMT_Pyro)
FOUR_CpGs$Pos <- as.factor(FOUR_CpGs$Pos)

FOUR_CpGs_Group <- FOUR_CpGs %>% 
  group_by(SampleID) %>% 
  summarise(Average_Nano=mean(Methylation_percent), 
            Average_Pyro=mean(Percent_Meth_Pyro),
            Average_Cov=mean(Valid_cov)) %>%
  inner_join(Samples) %>%
  filter(Average_Cov>2) 

# Depth represented as dot size
ggplot(FOUR_CpGs_Group, aes(x=Known_status, y=Average_Nano, color = Known_status))+
  geom_hline(yintercept = 10)+
  geom_beeswarm(aes(size = Average_Cov), alpha = 0.6, cex=3)+
  scale_size_continuous(name="Depth",breaks = c(5,10,50,100), range = c(1,10) )+
  scale_color_brewer(palette = "Set1", guide="none")+
  xlab("Known methylation status")+
  ylab("Nanopore methylation (%)")+
  facet_grid(.~Series)+
  theme_bw()+
  theme(aspect.ratio = 1.5)

# Universal dot size
ggplot(FOUR_CpGs_Group, aes(x=Known_status, y=Average_Nano, color = Known_status))+
  geom_hline(yintercept = 10)+
  geom_beeswarm(aes(size = Average_Cov), alpha = 0.8,size = 5, cex=4)+
  #scale_size_continuous(guide="none")+
  scale_color_brewer(palette = "Set1", guide="none")+
  xlab("Known methylation status (PSQ)")+
  ylab("Nanopore methylation (%)")+
  # facet_grid(.~Series)+
  theme_bw(base_size = 16)+
  coord_fixed(0.025)

# Make STP27 table

STP27 <- All_MGMT %>% 
  filter(Pos == 129466944 | Pos == 129467310) %>% 
  select(1:2,4) %>% 
  pivot_wider(names_from = Pos, values_from=Methylation_percent) %>%
  inner_join(Samples)

colnames(STP27) <- c("SampleID", "cg12434587", "cg12981137", "Series", "Method", "Known_status", "Diagnosis" )

ggplot(STP27, aes(y=cg12434587, x=cg12981137, color = Known_status)) +
  #  stat_function(fun = fun.1, color = "grey", size = 2) +
  geom_point(size = 5, alpha = 0.6)+
#  scale_size_continuous(name="Depth",breaks = c(5,15,30), range = c(1,8))+
  scale_color_brewer(name="Known Status", palette = "Set1")+
  # scale_x_sqrt(limits = c(0,100),breaks=c(1,10,25,50,75,100))+
  # scale_y_sqrt(limits = c(0,100),breaks=c(1,10,25,50,75,100))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw(base_size = 18)+
 #   facet_wrap(.~Series)+
  theme(aspect.ratio = 1)

############ Correlation plot

rownames(MGMT_long) <- MGMT_long$SampleID
column_to_rownames()
MGMT_long <- na.omit(MGMT_long)


corrplot(cor(MGMT_long),
         method = "color",
         order = "original",
         #hclust.method = "ward.D",
         tl.cex = 0.4,
         tl.col = 'black')

############# Headmap of all Samples

## Make annotation table
load("GitHub/MGMT_R/R-markdown/Data/Annotations.Rdata")

library(stringi)

# 
row_anno <- MGMT_RunSum %>%
  dplyr::select(SampleID,Known_status,Diagnosis) %>%
  distinct() %>% column_to_rownames("SampleID") 
%>%
  filter(Diagnosis == "Glioblastoma")

row_anno <- na.omit(row_anno[rownames(MGMT_long_ext_island),])
MGMT_long_ext_island <- MGMT_long_ext_island[rownames(row_anno),]

row_anno$Diagnosis <- ifelse(
  row_anno$Diagnosis == "Astrocytoma" | 
    row_anno$Diagnosis == "Oligodendroglioma" | 
    row_anno$Diagnosis == "Astrocytoma HG", "IDH_glioma", 
  ifelse(row_anno$Diagnosis == "Glioblastoma", "Glioblastoma", "Other"))

colnames(row_anno) <- c("Known status", "Diagnosis")

ann_colors <- list("Known status"=c("UnMethylated"="#377EB8", 
                                    "Methylated"="#E41A1C"), 
                   "Diagnosis"=c("Glioblastoma"="#31688EFF", 
                                 "IDH_glioma"="#440154FF", 
                                 "Other"="#35B779FF"))

library(pheatmap)
library(grid)
library(RColorBrewer)
library(tibble)

pheatmap(MGMT_long_ext_island,
                cluster_rows = T,
                cluster_cols = F,
                 clustering_method = "ward.D",
                 cutree_rows = 3,
                 treeheight_row = 30,
                 border_color = NA,
                 scale = "none",
                 drop_levels = T,
                 color = rev(brewer.pal(n = 10, name = "Spectral")),
                 fontsize_col = 8,
                 fontsize_row = 4,
                 fontsize = 14,
                 legend = T,
                 show_rownames = F,
                 legend_breaks = seq(0,100,10),
     # #           annotation_col = col_anno,
                 annotation_row = row_anno,
                 annotation_colors = ann_colors,
     #            # main = "All samples,  (n=148)",
                 legend_labels = c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")
)


