library(dplyr)
library(purrr)

## Also load MGMT.Rdata from github

###############
###  Create a single data frame from all the bedmethyl files
############

#Subsample <- read.delim("~/MGMT/Export/SubSampling/FIVE_reads_ALL/", header = F)

#files <- list.files(pattern = '\\.FIVE_reads.methyl.bed$', full.names = F)
files <- list.files(pattern = '\\methyl.bed$', full.names = F) # for stepwise subsampling

Subsample <- map_df(files, ~read.delim2(.x, header = F) %>% mutate(File = basename(.x)))

#Subsample$File <- gsub(".FIVE_reads.methyl.bed","",Subsample$File)
Subsample$File <- gsub(".methyl.bed","",Subsample$File) # for stepwise subsampling


Subsample <- Subsample %>% select(2,10:19) %>% relocate(11)

names(Subsample) <- c("SampleID", "Pos", "Valid_cov", "Methylation_percent", "N_mod", "N_canon", "N_other_mod", "N_delete", "N_fail", "N_diff", "N_nocall")

Subsample$Methylation_percent <- as.numeric(Subsample$Methylation_percent)

# Extract the 98 CpG sites in the MGMT CpG island
Subsample_Island <- Subsample %>% filter(between(Pos, 129466683, 129467448))

# Calculate mean and median coverage of each CpG site
cov <- Subsample_Island %>% 
  group_by(SampleID) %>% 
  summarise(AVG_cov = mean(Valid_cov), MED_cov = median(Valid_cov), Avg_meth = mean(Methylation_percent)) %>%
  separate(SampleID, c("SampleID","Sample"), sep="[.]")

#Subsample_RunSum <- inner_join(Samples,cov)

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

# Pull all matching samples from MGMT_long_Island that are represented in subsample long island

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
