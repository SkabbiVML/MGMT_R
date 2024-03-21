library(dplyr)
library(purrr)

## Also load MGMT.Rdata from github
## Make new All_MGMT

files <- list.files(path="../../../MGMT/Export/Recall_modkit_BedMethyl/",  pattern = "methyl.bed", full.names = T)

All_MGMT <- map_df(files, ~read.delim2(.x, header = F) %>% mutate(File = basename(.x)))

All_MGMT$File <- gsub(".cpg.methyl.bed","",All_MGMT$File)

All_MGMT <- All_MGMT %>% select(2,10:19) %>% relocate(11)

names(All_MGMT) <- c("SampleID", "Pos", "Valid_cov", "Methylation_percent", "N_mod", "N_canon", "N_other_mod", "N_delete", "N_fail", "N_diff", "N_nocall")

All_MGMT$Methylation_percent <- as.numeric(All_MGMT$Methylation_percent)

###############
###  Create a single data frame from all the bedmethyl files subsamples to 6 reads (median cov should be 5)
############


#files <- list.files(pattern = '\\.FIVE_reads.methyl.bed$', full.names = F)
files <- list.files(path="../../../MGMT/Export/SubSampling/6.reads.modbed/",  pattern = "methyl.bed", full.names = T) # for stepwise subsampling
#files <- list.files(path="E:/6.reads.modbed/",  pattern = "methyl.bed", full.names = T)

Subsample5 <- map_df(files, ~read.delim2(.x, header = F) %>% mutate(File = basename(.x)))

#Subsample$File <- gsub(".FIVE_reads.methyl.bed","",Subsample$File)
Subsample5$File <- gsub(".6_reads.methyl.bed","",Subsample5$File) # for stepwise subsampling


Subsample5 <- Subsample5 %>% select(2,10:19) %>% relocate(11)

names(Subsample5) <- c("SampleID", "Pos", "Valid_cov", "Methylation_percent", "N_mod", "N_canon", "N_other_mod", "N_delete", "N_fail", "N_diff", "N_nocall")

Subsample5$Methylation_percent <- as.numeric(Subsample5$Methylation_percent)

# Extract the 98 CpG sites in the MGMT CpG island
Subsample_Island_5 <- Subsample5 %>% dplyr::filter(between(Pos, 129466683, 129467448))

# Calculate mean and median coverage of each CpG site
cov5 <- Subsample_Island_5 %>% 
  group_by(SampleID) %>% 
  summarise(AVG_cov = mean(Valid_cov), MED_cov = median(Valid_cov), Avg_meth = mean(Methylation_percent)) 

###############
####3 Repeat for subsample 11

files <- list.files(path="../../../MGMT/Export/SubSampling/11.reads.modbed",  pattern = "methyl.bed", full.names = T) # for stepwise subsampling
#files <- list.files(path="E:/11.reads.modbed/",  pattern = "methyl.bed", full.names = T) # for stepwise subsampling


Subsample10 <- map_df(files, ~read.delim2(.x, header = F) %>% mutate(File = basename(.x)))

#Subsample$File <- gsub(".FIVE_reads.methyl.bed","",Subsample$File)
Subsample10$File <- gsub(".11_reads.methyl.bed","",Subsample10$File) # for stepwise subsampling


Subsample10 <- Subsample10 %>% select(2,10:19) %>% relocate(11)

names(Subsample10) <- c("SampleID", "Pos", "Valid_cov", "Methylation_percent", "N_mod", "N_canon", "N_other_mod", "N_delete", "N_fail", "N_diff", "N_nocall")

Subsample10$Methylation_percent <- as.numeric(Subsample10$Methylation_percent)

# Extract the 98 CpG sites in the MGMT CpG island
Subsample_Island_10 <- Subsample10 %>% dplyr::filter(between(Pos, 129466683, 129467448))

# Calculate mean and median coverage of each CpG site
cov10 <- Subsample_Island_10 %>% 
  group_by(SampleID) %>% 
  summarise(AVG_cov = mean(Valid_cov), MED_cov = median(Valid_cov), Avg_meth = mean(Methylation_percent)) 

######################

# pull samples that have 20 or more median coverage

cov_Full <- cov %>% filter(MED_cov > 19)

######################

#####################
# Get a list of all samples that have high readcount in the full list and are also represented in the subsampled files
#####################

List_of_subsampling <- intersect(intersect(cov_Full$SampleID, cov5$SampleID),cov10$SampleID)

All_MGMT_subsample_full <- All_MGMT_Island %>% dplyr::filter(SampleID %in% List_of_subsampling) %>%
  dplyr::select(c("SampleID", "Pos", "Methylation_percent")) %>%
  pivot_wider(names_from = Pos, values_from=Methylation_percent) %>%
  column_to_rownames("SampleID")

All_MGMT_subsample_10 <- Subsample_Island_10 %>% filter(SampleID %in% List_of_subsampling) %>%
  dplyr::select(c("SampleID", "Pos", "Methylation_percent")) %>%
  pivot_wider(names_from = Pos, values_from=Methylation_percent) %>%
  column_to_rownames("SampleID")

All_MGMT_subsample_5 <- Subsample_Island_5 %>% filter(SampleID %in% List_of_subsampling) %>%
  dplyr::select(c("SampleID", "Pos", "Methylation_percent")) %>%
  pivot_wider(names_from = Pos, values_from=Methylation_percent) %>%
  column_to_rownames("SampleID")

colnames(All_MGMT_subsample_full) <- c(1:98)
colnames(All_MGMT_subsample_10) <- c(1:98)
colnames(All_MGMT_subsample_5) <- c(1:98)

# remove rows with missing values from all data
 All_MGMT_subsample_5 <- na.omit(All_MGMT_subsample_5)
 All_MGMT_subsample_full <- All_MGMT_subsample_full[rownames(All_MGMT_subsample_5),]
 All_MGMT_subsample_10 <- All_MGMT_subsample_10[rownames(All_MGMT_subsample_5),]

##################
library(class)

MGMT_long_island <- na.omit(MGMT_long_island)

res <- pheatmap(MGMT_long_island,
                cluster_rows = T,
                cluster_cols = F,
                clustering_method = "ward.D",
                cutree_rows = 2,
                fontsize_row = 5,
                scale = "none",
                drop_levels = T,
)

MGMT_clusters <- data.frame(cluster=cutree(res$tree_row, k=2))

cl <- factor(MGMT_clusters$cluster)

class_full <- knn(MGMT_long_island, All_MGMT_subsample_full, cl, k = 2, l = 0, prob = FALSE, use.all = TRUE)

class_sub10 <- knn(MGMT_long_island, All_MGMT_subsample_10, cl, k = 2, l = 0, prob = FALSE, use.all = TRUE)

class_sub5 <- knn(MGMT_long_island, All_MGMT_subsample_5, cl, k = 2, l = 0, prob = FALSE, use.all = TRUE)

subsample_knn <- as.data.frame(cbind(class_full, class_sub10, class_sub5))

rownames(subsample_knn) <- rownames(All_MGMT_subsample_full)

###### clustering

res_full <- pheatmap(All_MGMT_subsample_full,
                cluster_rows = T,
                cluster_cols = F,
                clustering_method = "ward.D",
                cutree_rows = 3,
                scale = "none",
                drop_levels = T,
)

MGMT_clusters_full <- data.frame(cluster=cutree(res_full$tree_row, k=3))

MGMT_clusters_full$SampleID <- rownames(MGMT_clusters_full)

names(MGMT_clusters_full) <- c("cluster_full", "SampleID")


res_10 <- pheatmap(All_MGMT_subsample_10,
                     cluster_rows = T,
                     cluster_cols = F,
                     clustering_method = "ward.D",
                     cutree_rows = 3,
                     scale = "none",
                     drop_levels = T,
)

MGMT_clusters_10 <- data.frame(cluster=cutree(res_10$tree_row, k=3))

MGMT_clusters_10$SampleID <- rownames(MGMT_clusters_10)

names(MGMT_clusters_10) <- c("cluster_10", "SampleID")

res_5 <- pheatmap(All_MGMT_subsample_5,
                   cluster_rows = T,
                   cluster_cols = F,
                   clustering_method = "ward.D",
                   cutree_rows = 3,
                   scale = "none",
                   drop_levels = T,
)

MGMT_clusters_5 <- data.frame(cluster=cutree(res_5$tree_row, k=3))

MGMT_clusters_5$SampleID <- rownames(MGMT_clusters_5)

names(MGMT_clusters_5) <- c("cluster_5", "SampleID")

clust_compare <- left_join(MGMT_clusters_full,MGMT_clusters_10) %>% left_join(MGMT_clusters_5)

#####################
### FOR STEPWISE ONLY
###################

ggplot(cov,aes(x=MED_cov, y=Avg_meth, color=SampleID))+
  geom_line()+
  geom_point(size=3)+
  geom_vline(xintercept = 5)+
  scale_x_log10()+
  xlab("Median coverage in CpG-island")+
  ylab("Mean methylation (%) in CpG-island")+
  theme_bw(base_size = 16)
##########################
## Make a long frame for heatmap
######

Subsample_long_island <- Subsample_Island %>%
  select(c("SampleID", "Pos", "Methylation_percent")) %>%
  pivot_wider(names_from = Pos, values_from=Methylation_percent) %>%
  column_to_rownames("SampleID")

colnames(Subsample_long_island) <- c(1:98)

GBMs <- MGMT_RunSum %>% dplyr::select(SampleID,Diagnosis,Series) %>% dplyr::filter(Diagnosis == "Glioblastoma" ) %>% distinct()

Subsample_long_island_GBM <- Subsample_long_island[GBMs$SampleID,]

Subsample_long_island_GBM <- na.omit(Subsample_long_island_GBM)

row_anno_GBM <- row_anno[1]

#################
####### Heatmap of only GBM subsampled - only one not clustering the same
###################


pheatmap(Subsample_long_island_GBM,
         cluster_rows = T,
         cluster_cols = F,
         clustering_method = "ward.D",
         cutree_rows = 2,
         treeheight_row = 20,
         border_color = NA,
         scale = "none",
         drop_levels = T,
         color = rev(brewer.pal(n = 10, name = "Spectral")),
         fontsize_col = 8,
         fontsize_row = 8,
         fontsize = 12,
         legend = T,
         show_rownames = T,
         legend_breaks = seq(0,100,10),
         #         annotation_col = col_anno,
         annotation_row = row_anno_GBM,
         annotation_colors = ann_colors,
         legend_labels = c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")
         
)

# Pull all matching samples from MGMT_long_island that are represented in subsample long island

MGMT_matching <- MGMT_long_island[rownames(Subsample_long_island),]

MGMT_matching <- na.omit(MGMT_matching)

Subsample_long_island <- Subsample_long_island[rownames(MGMT_matching),]

##################
### Heatmaps of all matching samples - pretty much the same but first cut of tree is not the same
################

pheatmap(Subsample_long_island,
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
         fontsize = 12,
         legend = T,
         show_rownames = F,
         legend_breaks = seq(0,100,10),
         #  annotation_col = col_anno,
         annotation_row = row_anno,
         annotation_colors = ann_colors,
         main = "Subsampling to 5 reads per sample,  (n=115)",
         legend_labels = c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")
)

pheatmap(MGMT_matching,
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
         fontsize = 12,
         legend = T,
         show_rownames = F,
         legend_breaks = seq(0,100,10),
         #  annotation_col = col_anno,
         annotation_row = row_anno,
         annotation_colors = ann_colors,
         #            # main = "All samples,  (n=148)",
         legend_labels = c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")
)


##########################

FullHeatmap <- pheatmap(MGMT_matching,
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
                 fontsize_row = 6,
                 fontsize = 12,
                 legend = T,
                 show_rownames = T,
                 legend_breaks = seq(0,100,10),
                 #  annotation_col = col_anno,
                 annotation_row = row_anno,
                 annotation_colors = ann_colors,
                 #            # main = "All samples,  (n=148)",
                 legend_labels = c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")
)

FULL_clusters <- data.frame(cluster=cutree(FullHeatmap$tree_row, k=3))

FULL_clusters$SampleID <- rownames(FULL_clusters)

names(FULL_clusters) <- c("Fullcluster", "SampleID")

SubHeatmap <- pheatmap(Subsample_long_island,
                       cluster_rows = T,
                       cluster_cols = F,
                       clustering_method = "ward.D",
                       cutree_rows = 2,
                       treeheight_row = 30,
                       border_color = NA,
                       scale = "none",
                       drop_levels = T,
                       color = rev(brewer.pal(n = 10, name = "Spectral")),
                       fontsize_col = 8,
                       fontsize_row = 6,
                       fontsize = 12,
                       legend = T,
                       show_rownames = T,
                       legend_breaks = seq(0,100,10),
                       #  annotation_col = col_anno,
                       annotation_row = row_anno,
                       annotation_colors = ann_colors,
                       main = "Subsampling to 5 reads per sample,  (n=115)",
                       legend_labels = c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")
)

Sub_clusters <- data.frame(cluster=cutree(SubHeatmap$tree_row, k=3))

Sub_clusters$SampleID <- rownames(Sub_clusters)

names(Sub_clusters) <- c("Subcluster", "SampleID")

Clust_compare <- left_join(FULL_clusters,Sub_clusters)

Clust_compare <- na.omit(Clust_compare)

Clust_compare$Subcluster <- ifelse(Clust_compare$Subcluster == 1,"Unmethylated","Methylated")
Clust_compare$Fullcluster <- ifelse(Clust_compare$Fullcluster == 1,"Unmethylated","Methylated")

Clust_compare$Match <- ifelse(Clust_compare$Fullcluster == Clust_compare$Subcluster, "Match", "Mismatch")
