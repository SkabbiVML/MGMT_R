library(mixOmics)

dim(MGMT_long_island)

############### PCA

tune.pca.mgmt <- tune.pca(MGMT_long_island, scale = T)

plot(tune.pca.mgmt)


pca.mgmt <- pca(MGMT_long_island, ncomp = 3, center = T, scale = T)

plotIndiv(pca.mgmt,
          comp = c(2,1),
          ind.names = F,
          group = row_anno$Diagnosis,
          legend = T,
          ellipse = T,
          title = "Principal component analysis",
          legend.title = "Known status"
          )

grid.keep <- c(seq(5,60,5))

tune.spca.result <- tune.spca(MGMT_long_island, 
                              ncomp = 2,
                              folds = 5,
                              test.keepX = grid.keep,
                              nrepeat = 50)

spca.mgmt <- spca(MGMT_long_island,
                  ncomp = 2,
                  keepX = c(20,20),
                  scale = T)

plotIndiv(spca.mgmt,
          ind.names = F,
          group = row_anno$Known_status,
          legend = T
)

head(selectVar(pca.mgmt, comp = 2)$value)

plotLoadings(spca.mgmt, comp = 1)
#########
row_anno.meth <- row_anno %>% filter(Known_status=="Methylated") %>% filter(Diagnosis != "Other")
mgmt.meth <- MGMT_long_island[rownames(row_anno.meth),]

pca.mgmt <- pca(mgmt.meth, ncomp = 2, center = T, scale = T)

plotIndiv(pca.mgmt,
          ind.names = F,
          group = row_anno.meth$Diagnosis,
          legend = T
)

################### PLS-DA

plsda.mgmt <- plsda(MGMT_long_island, row_anno$Known_status)

plotIndiv(plsda.mgmt,
          ind.names = F,
          group = row_anno$Known_status,
          ellipse = F,
          legend = T,
          title = "PLS-DA, methylated only",
          legend.title = "Diagnosis")



##################

pheatmap(MGMT_long_island,
         cluster_rows = T,
         cluster_cols = T,
         clustering_method = "ward.D",
         cutree_rows = 3,
         cutree_cols = 5,
         treeheight_row = 30,
         border_color = NA,
         scale = "none",
         drop_levels = T,
         color = rev(brewer.pal(n = 10, name = "Spectral")),
         fontsize_col = 8,
         fontsize_row = 5,
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



################### PLS-DA

plsda.mgmt <- plsda(mgmt.meth, row_anno.meth$Diagnosis)

plotIndiv(plsda.mgmt,
          ind.names = F,
          group = row_anno.meth$Diagnosis,
          ellipse = T,
          legend = T)

splads.mgmt <- splsda(mgmt.meth, row_anno.meth$Diagnosis, ncomp = 2, keepX = c(20,10))

plotIndiv(splads.mgmt,
          ind.names = F,
          group = row_anno.meth$Diagnosis,
          ellipse = T,
          legend = T)

plotLoadings(splads.mgmt, comp = 1)

#####

Table <- MGMT_RunSum %>%
  dplyr::select(SampleID,Method,Diagnosis) %>%
  distinct() %>%
  dplyr::group_by(Method,Diagnosis) %>%
  tally() %>%
  spread(Diagnosis,n) %>%
  replace(is.na(.), 0) %>%
  adorn_totals(c("row", "col"))
