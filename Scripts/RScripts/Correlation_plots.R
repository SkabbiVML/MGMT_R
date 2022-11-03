setwd("~/MGMT/MGMT_radium/Shorty_calls/") # VML desktop
setwd("~/MGMT/MGMT_VML/Shorty_calls/") # VML desktop
setwd("~/MGMT_stuff") # Home desktop
setwd("D:/MGMT_stuff") # USB stick

library(dplyr)
library(data.table)

#read file paths

all_paths <- list.files(path ="~/MGMT/MGMT_radium/Shorty_calls/", pattern ="*.tsv",full.names = T )
all_paths <- list.files(path ="~/MGMT/MGMT_VML/Shorty_calls/", pattern ="*.tsv",full.names = T )

all_paths <- list.files(path ="~/MGMT_stuff/MGMT_VML/Shorty_calls/", pattern ="*.tsv",full.names = T )
all_paths <- list.files(path ="~/MGMT_stuff/MGMT_radium/Shorty_calls/", pattern ="*.tsv",full.names = T )

#USB stick
all_paths <- list.files(path ="MGMT_radium/Shorty_calls/", pattern ="*.tsv",full.names = T )
all_paths <- list.files(path ="MGMT_VML/Shorty_calls/", pattern ="*.tsv",full.names = T )



#read file content

all_content <- all_paths %>% lapply(fread, header =T, select = c(2:8))

#read file name

all_filenames <- all_paths %>% basename() %>% as.list()

#combine file content list and filename

all_lists <- mapply(c, all_content, all_filenames, SIMPLIFY = F)

# unlist and change column name

all_results <- rbindlist(all_lists, fill = T)

# change column name

names(all_results)[8] <- "Sample"

all_results$Sample <- gsub("Summarized_methylation_","",all_results$Sample)
all_results$Sample <- gsub("-shorty.tsv","",all_results$Sample)

Methylation_summary <- all_results

# Filter results to only contain CpGs within the Cas9 region 
Metsum_ROI <- Methylation_summary %>% filter(start > 131264784 & start < 131266094)

# Metsum_ROI$start <- as.factor(Metsum_ROI$start)
# Metsum_ROI$end <- as.factor(Metsum_ROI$end)
Metsum_ROI$num_motifs_in_group <- as.factor(Metsum_ROI$num_motifs_in_group)
Metsum_ROI$called_sites <- as.factor(Metsum_ROI$called_sites)
Metsum_ROI$called_sites_methylated <- as.factor(Metsum_ROI$called_sites_methylated)
Metsum_ROI$group_sequence <- as.factor(Metsum_ROI$group_sequence)


# # How many replicates you want of each row. How many CpG sites are in each motif
# duptimes <- Methylation_summary$num_motifs_in_group
# 
# # Create an index of the rows you want with duplications
# idx <- rep(1:nrow(Methylation_summary), duptimes)
# 
# #Use that index to genderate your new data frame
# MetSum <- Methylation_summary[idx,]





library(tidyr)

Metsum_ROI_trim <- Metsum_ROI[,c(1:3,6:8)]

Metsum_ROI_trim$feature <- ifelse(Metsum_ROI_trim$start < 131265505, "Promoter", ifelse(Metsum_ROI_trim$start > 131265560,"Intron1","Exon1"))

Metsum_ROI_trim$start <- as.factor(Metsum_ROI_trim$start)
Metsum_ROI_trim$end <- as.factor(Metsum_ROI_trim$end)

Metsum_ROI_long <- spread(Metsum_ROI_trim, Sample, methylated_frequency)

# How many replicates you want of each row. How many CpG sites are in each motif
duptimes <- Metsum_ROI_long$num_motifs_in_group

# Create an index of the rows you want with duplications
idx <- rep(1:nrow(Metsum_ROI_long), duptimes)

#Use that index to generate your new data frame
Metsum_ROI_long <- Metsum_ROI_long[idx,]

Metsum_t <- t(Metsum_ROI_long)

colnames(Metsum_t) <- c(1:113)
colnames(Metsum_t) <- c(-6:-1, 1:98,c("+1","+2", "+3","+4", "+5", "+6", "+7", "+8", "+9"))

#write.csv(Metsum_t, "Radium_MGMT_Methylation_summary.csv")
#write.csv(Metsum_t, "VML_MGMT_Methylation_summary.csv")

Radium <- read.csv("MGMT_Radium/Radium_MGMT_Methylation_summary.csv", header = T, row.names = 1)
colnames(Radium) <- c(-6:-1, 1:98,c("+1","+2", "+3","+4", "+5", "+6", "+7", "+8", "+9"))

VML <- read.csv("MGMT_VML/VML_MGMT_Methylation_summary.csv", header = T, row.names = 1)
colnames(VML) <- c(-6:-1, 1:98,c("+1","+2", "+3","+4", "+5", "+6", "+7", "+8", "+9"))

# ### Combining VML and Radium Samples 
Metsum_t2 <- rbind(VML[c(6:nrow(VML)),], Radium[c(6:nrow(Radium)),])

Metsum_t2_names <- rownames(Metsum_t2)
#Metsum_t2 <- as.data.frame(Metsum_t2[c(6:nrow(Metsum_t2)),])

Metsum_t2 <- sapply(Metsum_t2,as.numeric)


rownames(Metsum_t2) <- Metsum_t2_names #

Radium_2 <- as.data.frame(Radium[c(6:nrow(Radium)),])
Radium_2 <- sapply(Radium_2,as.numeric)
rownames(Radium_2) <- rownames(Radium[c(6:nrow(Radium)),])
colnames(Radium_2) <- c(-6:-1, 1:98,c("+1","+2", "+3","+4", "+5", "+6", "+7", "+8", "+9"))

VML_2 <- sapply(VML[c(6:nrow(VML)),], as.numeric)
rownames(VML_2) <- rownames(VML[c(6:nrow(VML)),])



Metsum_ROI_long <- Metsum_ROI_long[,c(2:21)]

# Data frame with column annotations.
mat_col <- data.frame(Feature = Metsum_ROI_long$feature)
rownames(mat_col) <- colnames(Metsum_t)
row_anno <- read.csv("Row_annotation.csv", header = T, row.names = 1)
row_anno <- row_anno[c(2:nrow(row_anno)),]

library(dplyr)
row_anno_VML <- row_anno %>% filter(Datset =="VML")
row_anno_VML$Datset <- NULL
row_anno_Radium <- row_anno %>% filter(Datset =="Radium")
row_anno_Radium$Datset <- NULL

ann_colors = list(
  PyroSeq = c(Unmethylated = "grey", Methylated = "black"),
  NPexon1 = c(Unmethylated = "grey", Methylated = "black"),
  Feature = c(Promoter = "yellow", Exon1 = "blue", Intron1 = "green")
)

ann_colors_single = list(
  PyroSeq = c(Unmethylated = "grey", Methylated = "black"),
  NPexon1 = c(Unmethylated = "grey", Methylated = "black")
)


library(corrplot)   

# correlation og CpGs
corrplot(cor(Radium_2), 
         method = "color",
         order = "original",
         #hclust.method = "ward.D",
         tl.cex = 0.5,
         tl.col = 'black')


# correlation of samples
corrplot(cor(t(Radium_2)), 
         method = "circle",
         order = "hclust",
         hclust.method = "ward.D",
         addrect = 2,
         tl.cex = 0.7,
         tl.col = 'black')

library(pheatmap)
library(RColorBrewer)
library(viridis)

pheatmap(t(Metsum_t2))

pheatmap(Metsum_t2)

###
pheatmap(Radium_2, 
         cluster_rows = T, 
         cluster_cols = F,
         clustering_method = "ward.D", #"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
         cutree_rows = 2,
         scale = "none",
         border_color = NA)

###


mat_breaks <- seq(0:10)
legend_breaks <- seq(0,100)

pheatmap(t(Metsum_t2), 
         cluster_rows = F, 
         cluster_cols = T,
         clustering_method = "ward.D",
         cutree_cols = 3,
         scale = "none",
         border_color = NA,
         color             = inferno(10),
         fontsize_col = 10,
         fontsize_row = 6,
         legend = T,
         legend_breaks = seq(0,1,0.1),
         legend_labels = c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"),
         angle_col = 45,
         annotation_row = mat_col,
         )

### Heatmap for all data

pheatmap(Metsum_t2, 
         cluster_rows = T, 
         cluster_cols = F,
         clustering_method = "ward.D",
         cutree_rows = 2,
         scale = "none",
         border_color = NA,
         #color = inferno(10),
         #color = rev(brewer.pal(n = 11, name = "RdBu")),
         fontsize_col = 6,
         fontsize_row = 10,
         legend = T,
         legend_breaks = seq(0,1,0.1),
         legend_labels = c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"),
        # angle_col = 315,
         annotation_col = mat_col,
        annotation_row = row_anno,
        annotation_colors = ann_colors
)

### Heatmap for each dataset, VML_2 or Radium_2. must change annotation_row

pheatmap(Radium_2, 
         cluster_rows = T, 
         cluster_cols = F,
         clustering_method = "ward.D", #"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
         cutree_rows = 2,
         scale = "none",
         border_color = NA,
         #color = inferno(10),
         #color = rev(brewer.pal(n = 11, name = "RdBu")),
         fontsize_col = 6,
         fontsize_row = 10,
         legend = T,
         legend_breaks = seq(0,1,0.1),
         legend_labels = c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"),
         # angle_col = 315,
         annotation_col = mat_col,
         #annotation_row = row_anno_VML,
         annotation_row = row_anno_Radium,
         annotation_colors = ann_colors,
         main = "Radium samples"
)
