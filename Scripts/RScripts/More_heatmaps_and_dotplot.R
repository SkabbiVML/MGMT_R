colnames(Metsum_ROI_long) <- c(-7:-1, 1:98,c("+1","+2", "+3","+4", "+5", "+6", "+7", "+8", "+9", "+10"))
load("MGMT_R/R-Markdown/Data/Annotations.Rdata")

library(pheatmap)
library(grid)
library(RColorBrewer)


##### Removing DenStem Samples
Metsum_ROI_long_DENREMOVE <- Metsum_ROI_long[c(2:128),]

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

GBMs <- MGMT_RunSum %>% dplyr::select(SampleID,Diagnosis,Series) %>% filter(Diagnosis == "Glioblastoma") %>% distinct()

Metsum_ROI_long_GBM <- na.omit(Metsum_ROI_long[GBMs$SampleID,])

row_anno_GBM <- row_anno %>% dplyr::select(Known_status)

res <- pheatmap(Metsum_ROI_long_GBM,
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
                annotation_row = row_anno_GBM,
                annotation_colors = ann_colors,
                # main = "All samples,  (n=148)",
                legend_labels = c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")
)

X <- colSums(Metsum_ROI_long)/145

X2 <- rbind(colnames(Metsum_ROI_long),X)

X3 <- as.data.frame(t(X2))

names(X3) <- c("CpG", "Averag_Meth")

X3$Averag_Meth <- as.numeric(X3$Averag_Meth)
X3$CpG <- factor(X3$CpG, levels = unique(X3$CpG))
X3$CpG <- as.integer(X3$CpG)

ggplot(X3,aes(x=CpG,y=Averag_Meth))+
  geom_point()+
  geom_line()
