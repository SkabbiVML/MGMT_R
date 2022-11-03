all_results_MGMT <- read.csv("F:/MGMT_stuff/Aggregated_MGMT_05102022.csv", row.names = 1)
all_results_MGMT_full <- read.csv("D:/MGMT/Results/Aggregated_MGMT_08102022.csv", row.names = 1)
all_results_MGMT_full <- read.csv("Data/Aggregated_MGMT_24102022.csv", row.names = 1)

library(tidyr)
library(dplyr)
library(readr)
library(corrplot)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(viridis)

#
#all_results_MGMT <- all_results_MGMT %>% filter(Calls_num > 9)

# Trim off Dataset,call numbers and features. add these back for annotation of the plots
all_results_MGMT <- all_results_MGMT_full[,c(2:3,6)]


Metsum_ROI_long <- all_results_MGMT %>% pivot_wider(names_from = c(Pos), values_from = Percent_Methylated)
Metsum_ROI_long <- column_to_rownames(Metsum_ROI_long,var="Sample")

#### Some of the samples have missing values, mostly due to low coverage
#### There are a few that have decent coverage and only one missing value that 
#### means they have to be dropped from downstream analysis.
#### I've chosen to keep samples that only have a single missing value and
#### impute the missing value by the median methylation of the entire sample.
#### This can be discussed

# I've also run this with only complete cases which excludes a sample if there is a single NA
# somewhere in the data. I felt this was too stringent as is threw out some good samples
# Metsum_ROI_long[complete.cases(Metsum_ROI_long),]

library(Hmisc)

to_keep <- rownames(Metsum_ROI_long[rowSums(is.na(Metsum_ROI_long))< 2,])

Metsum_ROI_long <- Metsum_ROI_long[to_keep,]

Metsum_t <- as.data.frame(t(Metsum_ROI_long))

Metsum_t_impute <- impute(Metsum_t, median)

Metsum_ROI_long <- as.data.frame(t(Metsum_t_impute))

#### renaming chromosome position to CpG number

colnames(Metsum_ROI_long) <- c(-7:-1, 1:98,c("+1","+2", "+3","+4", "+5", "+6", "+7", "+8", "+9", "+10"))

########## Filtering MEDIUM and HIGH coverage samples
#MGMT_RunSum <- read.csv("D:/MGMT/SampleMetadata/MGMT_RunSum.csv", sep = ",", header = T) 

MGMT_RunSum <- read.csv("/media/vml/Ghibli/MGMT/SampleMetadata/MGMT_RunSum.csv", sep = ",", header = T) # Ghibli# Ghibli

### for simplification
MGMT_RunSum$Series <- gsub("Radium1", "Radium", MGMT_RunSum$Series)
MGMT_RunSum$Series <- gsub("Radium2", "Radium", MGMT_RunSum$Series)


Medium_coverage_samples <- MGMT_RunSum %>% filter(On_target_seqs > 9) %>% select(Sample_ID)

Medium_coverage_samples <- unique(Medium_coverage_samples$Sample_ID) 

High_coverage_samples <-  MGMT_RunSum %>% filter(On_target_seqs > 29) %>% select(Sample_ID)

High_coverage_samples <- unique(High_coverage_samples$Sample_ID)

Metsum_ROI_long_MEDIUM <- Metsum_ROI_long[Medium_coverage_samples,]
Metsum_ROI_long_HIGH <- Metsum_ROI_long[High_coverage_samples,]

################## Plotting

# correlation of CpGs
corrplot(cor(Metsum_ROI_long), 
         method = "color",
         order = "original",
         #hclust.method = "ward.D",
         tl.cex = 0.4,
         tl.col = 'black')

# correlation of samples
corrplot(cor(t(Metsum_ROI_long)), 
         method = "color",
         order = "hclust",
         hclust.method = "ward.D",
         addrect = 3,
         tl.cex = 0.3,
         tl.col = 'black')

######## Heatmaps

## Annotation tables

#feature annotation
mat_col <- all_results_MGMT_full[,c(3,7)] 
mat_col$Pos <- as.factor(mat_col$Pos)
mat_col <- mat_col %>%  group_by(Pos,feature) %>% summarise()
mat_col <- data.frame(Feature =mat_col$feature)
rownames(mat_col) <- colnames(Metsum_ROI_long)

# sample annotation
row_anno <- unique(MGMT_RunSum[,c("Sample_ID","Series","Pyro_Methylation_Status")])
rownames(row_anno) <- NULL
row_anno <- column_to_rownames(row_anno, var="Sample_ID")
colnames(row_anno) <- c("Series", "Known_status")
row_anno$Series <- gsub("-","_",row_anno$Series)

# annotation colors

ann_colors = list(
 #Series = c(DenStem = "#33A02C", Radium1 = "#1F78B4", Radium2 = "#A6CEE3", Rapid_CNS = "#B2DF8A"),
  Known_status = c(UnMethylated = "grey", Methylated = "black"),
  Series = c(DenStem = "#377EB8",  Radium = "#E41A1C", Rapid_CNS = "#440154"),
  Feature = c(Promoter = "#7FC97F", Exon1 = "#377EB8", Intron1 = "#FDC086")
)

# ann_colors = list(
#   Series = c(DenStem = "#33A02C", Radium = "#1F78B4", Rapid_CNS = "#B2DF8A"),
#   Meth_Status = c(UnMethylated = "grey", Methylated = "black", Unknown ="white"),
#   Feature = c(Promoter = "#7FC97F", Exon1 = "#377EB8", Intron1 = "#FDC086")
# )


library(pheatmap)
library(grid)

pheatmap(Metsum_ROI_long,
         cluster_rows = T, 
         cluster_cols = F,
         clustering_method = "ward.D",
         cutree_rows = 4,
         scale = "none",
         border_color = "grey",
         color = rev(brewer.pal(n = 10, name = "Spectral")),
         fontsize = 6,
         fontsize_col = 4,
         fontsize_row = 4,
        # kmeans_k = 4,
        show_rownames = F,
         treeheight_row = 30,
         legend = T,
         legend_breaks = seq(0,100,10),
          annotation_col = mat_col,
          annotation_row = row_anno,
         annotation_colors = ann_colors,
         #main = "Sequence depth > 29,  (n=45)",
         legend_labels = c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")
         )


Metsum_ROI_long_ISLAND <- Metsum_ROI_long[c(2:nrow(Metsum_ROI_long)),c(8:105)]

pheatmap(Metsum_ROI_long_ISLAND,
         cluster_rows = T, 
         cluster_cols = F,
         clustering_method = "ward.D",
         cutree_rows = 2,
         scale = "none",
         border_color = NA,
         color = rev(brewer.pal(n = 10, name = "Spectral")),
         fontsize = 6,
         fontsize_col = 4,
         fontsize_row = 4,
         # kmeans_k = 4,
         show_rownames = T,
         treeheight_row = 30,
         legend = T,
         legend_breaks = seq(0,100,10),
         annotation_col = mat_col,
         annotation_row = row_anno,
         annotation_colors = ann_colors,
         #main = "Sequence depth > 29,  (n=45)",
         legend_labels = c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")
)

GBMs <- MGMT_RunSum %>% select(Sample_ID,Diagnosis) %>% filter(Diagnosis == "Glioblastoma") %>% distinct()

Metsum_ROI_long_GBM <- na.omit(Metsum_ROI_long_ISLAND[GBMs$Sample_ID,])

pheatmap(Metsum_ROI_long_GBM,
         cluster_rows = T, 
         cluster_cols = F,
         clustering_method = "ward.D",
         cutree_rows = 6,
         scale = "none",
         border_color = NA,
         color = rev(brewer.pal(n = 10, name = "Spectral")),
         fontsize = 6,
         fontsize_col = 4,
         fontsize_row = 4,
         # kmeans_k = 4,
         #annotation = "CpG",
         show_rownames = T,
         treeheight_row = 30,
         legend = T,
         legend_breaks = seq(0,100,10),
         annotation_col = mat_col,
         annotation_row = row_anno,
         annotation_colors = ann_colors,
         main = "Glioblastoma samples ONLY",
         legend_labels = c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")
)
