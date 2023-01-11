ggplot(MGMT_RunSum) +
  geom_bar(aes(x=Diagnosis,fill=Pyro_Methylation_Status,colour="black"))+
  scale_fill_viridis_d(name="Mehtylation Status")+
  theme_bw()
######################

res <- pheatmap(Metsum_ROI_long_ISLAND,
                cluster_rows = T, 
                cluster_cols = F,
                clustering_method = "ward.D",
                cutree_rows = 4,
                scale = "none",
                border_color = NA,
                color = rev(brewer.pal(n = 10, name = "Spectral")),
                fontsize = 6,
                fontsize_col = 4,
                fontsize_row = 6,
                #  kmeans_k = 6,
                show_rownames = T,
                treeheight_row = 20,
                legend = T,
                legend_breaks = seq(0,100,10),
                annotation_col = mat_col,
                # angle_col = 45,
                annotation_row = row_anno,
                annotation_colors = ann_colors,
                #main = "Glioblastoma samples ONLY (n=91)",
                legend_labels = c("0%","10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")
)

####################################



Patient_meta <- read.csv("Data/Patient_Survival.csv")

MGMT_clusters <- data.frame(cluster=cutree(res$tree_row, k=2))

MGMT_clusters$Sample_ID <- rownames(MGMT_clusters)

cluster_surv <- inner_join(MGMT_clusters, Patient_meta, by="Sample_ID")

cluster_surv <- cluster_surv %>% filter(Diagnosis == "GBM")

cluster_surv$months_to_event <- as.numeric(cluster_surv$months_to_event)

cluster_surv_filt <- cluster_surv %>% filter(Reopr == 0)

####################### Survival
library(survival)
library(survminer)
library(lubridate)

# Unfiltered
fit1 <- survfit( Surv(cluster_surv$months_to_event, cluster_surv$censor) ~ cluster_surv$Pyro_state)
fit2 <- survfit( Surv(cluster_surv$months_to_event, cluster_surv$censor) ~ cluster_surv$cluster)


# Reopr filtered out
fit1_filt <- survfit( Surv(cluster_surv_filt$months_to_event, cluster_surv_filt$censor) ~ cluster_surv_filt$Pyro_state)
fit2_filt <- survfit( Surv(cluster_surv_filt$months_to_event, cluster_surv_filt$censor) ~ cluster_surv_filt$cluster)




# Kaplan-Meier plot

################ Unfiltered
# By pyrosequencing results
ggsurvplot(
  fit = fit1,
  data = cluster_surv,
  palette = c("cornflowerblue","darkorange3"),
  xlab = "OS, Months", 
  ylab = "Overall survival probability",
  legend.title = "Methylation state",
  legend.labs = c("Methylated", "Unmethylated"),
  pval = TRUE,
  pval.coord = c(0, 0.1),
  font.x = c(14, "bold", "black"),
  font.y = c(14, "bold", "black"),
  risk.table = TRUE)

# By nanopore clustering
ggsurvplot(
  fit = fit2,
  data = cluster_surv,
  palette = c("darkorange3","cornflowerblue"),
  xlab = "OS, Months", 
  ylab = "Overall survival probability",
  legend.title = "Nanopore cluster",
  #legend.labs = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"),
  legend.labs = c("Cluster 1", "Cluster 2"),
  pval = TRUE,
  pval.coord = c(0, 0.1),
  font.x = c(14, "bold", "black"),
  font.y = c(14, "bold", "black"),
  risk.table = TRUE)

############################## Filtered
# By pyrosequencing results
ggsurvplot(
  fit = fit1_filt,
  data = cluster_surv_filt,
  palette = c("cornflowerblue","darkorange3"),
  xlab = "OS, Months", 
  ylab = "Overall survival probability",
  legend.title = "Classification by Pyrosequencing",
  legend.labs = c("Methylated", "Unmethylated"),
  pval = TRUE,
  pval.coord = c(0, 0.1),
  font.x = c(14, "bold", "black"),
  font.y = c(14, "bold", "black"),
  risk.table = TRUE)

# By nanopore clustering
ggsurvplot(
  fit = fit2_filt,
  data = cluster_surv_filt,
  
  palette = c("darkorange3","cornflowerblue"),
  xlab = "OS, Months", 
  ylab = "Overall survival probability",
  legend.title = "Classification by Nanopore cluster",
  #legend.labs = c("Cluster 1", "Cluster 2", "Cluster 3"),
  legend.labs = c("Cluster 1", "Cluster 2"),
  pval = TRUE,
  pval.coord = c(0, 0.1),
  font.x = c(14, "bold", "black"),
  font.y = c(14, "bold", "black"),
  risk.table = TRUE)

##################### 

Metsum_ROI_long_DMR2 <- Metsum_ROI_long[c(2:nrow(Metsum_ROI_long)),c(76:79)]

