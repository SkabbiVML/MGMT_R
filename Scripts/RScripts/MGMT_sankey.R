# This script generates a sankey plot of pyrosequencing results compared to nanopore results
library(dplyr)
library(tidyr)
library(tibble)

FOUR_CpGs_MGMT <- read.csv("/media/vml/Ghibli/MGMT/Data/FOUR_CpGs_MGMT.csv", header = T, row.names = 1)

test <- FOUR_CpGs_MGMT %>%  
  gather(key = "Pos", value = percent_methylated, -c(1,2)) %>%
  group_by(Sample_ID, DataSet) %>%
  mutate(Nanopore_average = mean(percent_methylated)) %>%
  select(c(1,2,5)) %>%
  unique()

test$Nanopore_status <- ifelse(test$Nanopore_average > 9, "Methylated", "UnMethylated")

Q <- MGMT_RunSum %>% select(Sample_ID,Pyro_Methylation_Status)

test <- na.omit(left_join(test,Q)) %>% select(Sample_ID, Nanopore_status, Pyro_Methylation_Status)

test <- as_tibble(test)

test$Concordance <- ifelse(test$Nanopore_status == test$Pyro_Methylation_Status, "concordant", "discordant")

test$Nanopore_status <- gsub("Methylated", "Methylated ", test$Nanopore_status)

links <- test %>% select(Nanopore_status, Pyro_Methylation_Status, Concordance)

links <- links %>% unite(First_key, Pyro_Methylation_Status, Nanopore_status,  sep="->", remove = F) %>%
  unite(Second_key, Nanopore_status, Concordance, sep="->", remove = F) %>%
  select(First_key, Second_key, Concordance) %>% 
  gather(Keys, Concordance) %>%
  select(Concordance) %>% group_by(Concordance) %>% tally(sort = TRUE) %>%
  separate(Concordance, c("Source","Target"), sep="->")

#change the column names to match the Sankey inbut
colnames(links)<-c("Source", "Target", "value")

# Create a node data frame: it lists all entities involved in the Sankey
nodes <- data.frame(
  name=c(as.character(links$Source), 
         as.character(links$Target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe. Need to reformat it.
links$IDsource <- match(links$Source, nodes$name)-1 
links$IDtarget <- match(links$Target, nodes$name)-1

p <- sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", 
              NodeID = "name",
              fontSize = 22,
              fontFamily = "sans-serif",
              nodeWidth = 15,
              nodePadding = 50,
              sinksRight=F,
              margin = c(top=100,right=0, left=0,bottom=100),
              iterations = 100)


p