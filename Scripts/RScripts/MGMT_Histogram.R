#Binarizing and summarizing methylation of CpGs

library(ggplot2)
library(dplyr)
library(janitor)

# Use Metsum_ROI_long_ISLAND from main R-markdown

t <- ifelse(Metsum_ROI_long_ISLAND > 10, 1,0) # If methylation of any given CpG (column), in any given sample (row) is above 10% than 1, else 0

t2 <- as.data.frame(rowSums(t))

names(t2) <-"CumMeth"

t2$SampleID <- rownames(t2)

row_anno$SampleID <- rownames(row_anno)

t3 <- left_join(t2, row_anno)

coverage <- MGMT_RunSum %>% select(c(Sample_ID,On_target_seqs))
names(coverage) <- c("Series", "Known_status", "SampleID", "Coverage")

t3 <- left_join(t3, coverage, by="SampleID")

# rough histogram
ggplot(t3, aes(x=CumMeth, fill=Known_status)) + geom_histogram(binwidth = 10, alpha=.5, position = "identity")

# rough density plot
ggplot(t3, aes(x=CumMeth, fill=Known_status)) + 
  geom_density(alpha=.5, position = "identity")+
  xlab("Cumulative methylation")+
  theme_bw()

# rough box plot
ggplot(t3, aes(y=CumMeth, x=Known_status, fill=Series)) + geom_boxplot()

# rough dotplot
ggplot(t3, aes(y=CumMeth, x=Known_status, fill=Series, color=Series, group=Series)) + 
  geom_point(position = position_dodge(width = .5))
 # geom_text(aes(label=SampleID))

# Beeswarm based on cumulative methylation

ggplot(t3, aes(x=Known_status.x, y=CumMeth, color = Series.y))+
  geom_hline(yintercept = 45)+
  geom_beeswarm(aes(size = Coverage), alpha = 0.6, cex=3)+
  scale_size_continuous(name="Depth",breaks = c(5,10,50,100), range = c(1,10) )+
  scale_color_brewer(palette = "Set1", guide="none")+
  xlab("Known methylation status")+
  ylab("Cumulative methylation per nanopore")+
  facet_grid(.~Series.y)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1), aspect.ratio = 2.5)
