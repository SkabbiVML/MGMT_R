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
  geom_point(position = position_dodge(width = .5))+
  geom_text(aes(label=SampleID))
