Metsum_ROI_long <- read.csv("Data/Metsum_ROI_long.csv", header = T, row.names = 1)
colnames(Metsum_ROI_long) <- c(-7:-1, 1:98,c("+1","+2", "+3","+4", "+5", "+6", "+7", "+8", "+9", "+10"))
load("Data/Annotations.Rdata")

library(pheatmap)
library(grid)
library(RColorBrewer)

Metsum_ROI_long_ISLAND <- Metsum_ROI_long[c(2:nrow(Metsum_ROI_long)),c(8:105)]

##### Removing DenStem Samples
Metsum_ROI_long_ISLAND <- Metsum_ROI_long_ISLAND[c(1:128),]

row_anno$SampleID <- rownames(row_anno)

Diagnosis <- distinct(MGMT_RunSum[,c("SampleID", "Diagnosis")])

row_anno <- inner_join(row_anno,Diagnosis)

row_anno <- row_anno %>% mutate(Diagnosis = case_when(Diagnosis == 'Glioblastoma' ~ 'Glioblastoma',
                                                  Diagnosis == 'Meningioma' ~ 'Meningioma',    
                                                  Diagnosis == 'Astrocytoma' ~ 'IDHglioma',
                                                  Diagnosis == '`Astrocytoma HG`' ~ 'IDHglioma',
                                                  Diagnosis == 'Oligodendroglioma' ~ 'IDHglioma',
                                                  TRUE ~ 'Other'))

rownames(row_anno) <- row_anno$SampleID

row_anno <- row_anno[,c("Known_status","Diagnosis")]

ann_colors = list(
  Known_status = c(UnMethylated = "blue", Methylated = "red"),
  Diagnosis = c(Glioblastoma = "#31688EFF", IDHglioma = "#440154FF" , Meningioma = "#FDE725FF", Other = "#35B779FF"),
  Feature = c(Promoter = "#7FC97F", Exon1 = "#377EB8", Intron1 = "#FDC086"),
  Method = c(STP27 = "black", Siller_et_al = "darkgrey", " "  = "white")
)

mat_col$Method <- ifelse(row.names(mat_col)==84, "STP27", 
                         ifelse(row.names(mat_col)==31, "STP27", 
                                ifelse(row.names(mat_col)>=74, "Siller_et_al", " ")))
mat_col[c(15,16),2] <- " "