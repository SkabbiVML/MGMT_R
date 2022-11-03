setwd("~/")

library(dplyr)
library(tidyr)
library(data.table)
library(janitor)
library(stringr)

#read file paths

# On Ghibli
all_paths_MGMT <- list.files(path ="/MGMT/Results/", pattern = ".bed", full.names = T, recursive = T )

# On Precilla
setwd("~/")
all_paths_MGMT <- list.files(path ="MGMT/Analysis/Results/", pattern = ".bed", full.names = T, recursive = T )

#read file content

all_content_MGMT <- all_paths_MGMT %>% lapply(fread, header =F, select = c(1:2,6,10:11))

#All <- lapply(all_paths_MGMT, read.delim, header = F) %>% bind_rows()

#read file name

all_filenames <- all_paths_MGMT %>% dirname() %>% as.list()

# all_filenames_MGMT <- all_paths_MGMT %>% gsub("MGMT/Results/Areeba/","") %>% as.list()
# all_filenames_MGMT <- gsub("MGMT/Analysis/Results/Areeba/","",all_paths_MGMT) 
# all_filenames_MGMT <- gsub("/modified_bases.5mC.bed","",all_filenames_MGMT) 


#combine file content list and filename

all_lists_MGMT <- mapply(c, all_content_MGMT, all_filenames, SIMPLIFY = F)

# unlist and change column name

#all_results_MGMT <- rbind.data.frame(all_lists_MGMT, simplify = F)
all_results_MGMT <- rbindlist(all_lists_MGMT, fill = T)

colnames(all_results_MGMT) <- c("Chromosome", "Pos","Strand","Calls_num","Meth_percent","loc")
all_results_MGMT <- all_results_MGMT %>%  filter(between(Pos, 131264784, 131266094))

all_results_MGMT_extralong <- all_results_MGMT %>%  filter(between(Pos, 131264784, 131269094)) #extralong

#all_results_MGMT$loc <- gsub("/MGMT/Results/", "", all_results_MGMT$loc) # Data structure on Ghibli
all_results_MGMT$loc <- gsub("MGMT/Analysis/Results//", "", all_results_MGMT$loc) # Data structure on Precilla
all_results_MGMT$loc <- gsub("/MGMT/Results//", "", all_results_MGMT$loc) # Data structure on Ghibli

all_results_MGMT <- separate(all_results_MGMT,loc,c("DataSet","Sample"),sep = "/")

##### Fixing sample duplicates
all_results_MGMT$Sample <- gsub("2159b","2159",all_results_MGMT$Sample)
all_results_MGMT$Sample <- gsub("_A","",all_results_MGMT$Sample)
all_results_MGMT$Sample <- gsub("_2","",all_results_MGMT$Sample)

all_results_MGMT <- all_results_MGMT %>% group_by(DataSet,Sample,Pos) %>% summarise(across(c(Calls_num,  Meth_num), sum)) %>% ungroup()


# Calculate number of calls for each location

all_results_MGMT <- all_results_MGMT %>% mutate(Meth_num = round(Calls_num * Meth_percent /100))

# Shift the position of negative strands by one so that they match the positive strand

all_results_MGMT <- all_results_MGMT %>% mutate(Pos = ifelse(Strand == "-", -1 + Pos, Pos))

# Add the values from positive and negative strands
all_results_MGMT <- all_results_MGMT %>% group_by(DataSet,Sample,Pos) %>% summarise(across(c(Calls_num,  Meth_num), sum)) %>% ungroup()

# Create a new percent methylated column from the new Calls_num and Meth_num values

all_results_MGMT <- all_results_MGMT %>% mutate(Percent_Methylated = Meth_num/Calls_num * 100)

all_results_MGMT$feature <- ifelse(all_results_MGMT$Pos < 131265505, "Promoter", ifelse(all_results_MGMT$Pos > 131265560,"Intron1","Exon1"))

all_results_MGMT$Sample <- as.factor(all_results_MGMT$Sample)
all_results_MGMT$DataSet <- as.factor(all_results_MGMT$DataSet)
all_results_MGMT$feature <- as.factor(all_results_MGMT$feature)
all_results_MGMT$Pos <- as.factor(all_results_MGMT$Pos)

# continue with plotting in script "MGMT_plots"

write.csv(all_results_MGMT, "F:/MGMT_stuff/Aggregated_MGMT_05102022.csv")
