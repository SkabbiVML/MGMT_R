colnames(Metsum_ROI_long) <- c(-7:-1, 1:98,c("+1","+2", "+3","+4", "+5", "+6", "+7", "+8", "+9", "+10"))
load("MGMT_R/R-Markdown/Data/Annotations.Rdata")

library(pheatmap)
library(grid)
library(RColorBrewer)


##### Removing DenStem Samples
Metsum_ROI_long_DENREMOVE <- Metsum_ROI_long[c(2:129),]

res <- pheatmap(Metsum_ROI_long_DENREMOVE,
                cluster_rows = T,
                cluster_cols = F,
                clustering_method = "ward.D",
                cutree_rows = 2,
                treeheight_row = 100,
                border_color = NA,
                scale = "none",
                drop_levels = F,
                color = rev(brewer.pal(n = 10, name = "Spectral")),
                fontsize_col = 8,
                fontsize_row = 4,
                fontsize = 14,
                legend = T,
                show_rownames = F,
                legend_breaks = seq(0,100,10),
                annotation_col = col_anno,
                annotation_row = row_anno,
                annotation_colors = ann_colors,
                # main = "All samples,  (n=148)",
                legend_labels = c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")
)


##### Only GBMs

GBMs <- MGMT_RunSum %>% dplyr::select(SampleID,Diagnosis,Series) %>% filter(Diagnosis == "Glioblastoma" & Series != "DenStem") %>% distinct()

Metsum_ROI_long_GBM <- na.omit(Metsum_ROI_long[GBMs$SampleID,])

row_anno_GBM <- row_anno %>% dplyr::select(Known_status)

res <- pheatmap(Metsum_ROI_long_GBM,
                cluster_rows = T,
                cluster_cols = F,
                clustering_method = "ward.D",
                cutree_rows = 3,
                treeheight_row = 100,
                border_color = NA,
                scale = "none",
                drop_levels = F,
                color = rev(brewer.pal(n = 10, name = "Spectral")),
                fontsize_col = 8,
                fontsize_row = 4,
                fontsize = 14,
                legend = T,
                show_rownames = F,
                legend_breaks = seq(0,100,10),
                annotation_col = col_anno,
                annotation_row = row_anno_GBM,
                annotation_colors = ann_colors,
                # main = "All samples,  (n=148)",
                legend_labels = c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")
)


res <-  pheatmap(Metsum_ROI_long_GBM,
         cluster_rows = T,
         cluster_cols = F,
         clustering_method = "ward.D",
         cutree_rows = 3,
         treeheight_row = 100,
         border_color = NA,
         scale = "none",
         kmeans_k = 4,
         drop_levels = F,
         color = rev(brewer.pal(n = 10, name = "Spectral")),
         fontsize_col = 8,
         fontsize_row = 4,
         fontsize = 14,
         legend = T,
         show_rownames = F,
         legend_breaks = seq(0,100,10),
         annotation_col = col_anno,
        # annotation_row = row_anno_GBM,
         annotation_colors = ann_colors,
         # main = "All samples,  (n=148)",
         legend_labels = c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")
)


res$kmeans

#################### Dotplot of average methylation per CpG ######################

x <- Metsum_ROI_long[c(2:129),]

colnames(x) <- c(-7:-1, 1:108)

x$SampleID <- rownames(x)

y2 <- MGMT_RunSum %>% dplyr::select(SampleID,Known_status)

x2 <- left_join(x, y2) %>% na.omit() %>% gather(key = "CpG", value = "percentMeth", -c(Known_status,SampleID))

x3 <- x2 %>% group_by(Known_status, CpG) %>% summarise(Average_Methylation = mean(percentMeth), SD_Methylation = sd(percentMeth))


x3$CpG <- as.integer(x3$CpG)

ggplot(x3, aes(x=CpG, y=Average_Methylation, group =  Known_status))+
  geom_point(aes(color=Known_status), position = position_dodge(width = 0.9), size = 2)+
#  guides(color = guide_legend(override.aes = list(size=2)))+
  geom_linerange(aes(ymin=Average_Methylation-SD_Methylation, ymax=Average_Methylation+SD_Methylation),position = position_dodge(width = 0.9))+
  scale_color_brewer(palette = "Set1")+
  geom_line(aes(color=Known_status))+
  geom_smooth(aes(color=Known_status), show.legend = FALSE)+
  theme_bw()

############################################### Dendrogram

