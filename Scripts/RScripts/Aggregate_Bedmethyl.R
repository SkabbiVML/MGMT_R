#Load the aggregated bed.methyl file

Bed <- read.delim("/encrypteddata/MGMT/mgmt_retest_01082022/results/t21-240/modified_bases.5mC.bed")

# Remove duplicate columns and the ones I don't understand

Bed <- Bed[,c(1:2,6,10:11)]

# Give names to the columns

colnames(Bed) <- c("Chromosome","Pos","Strand","Calls_num","Meth_percent")

Bed <- as_tibble(Bed)

# Filter out irrelevant locations
library(dplyr)

Bed2 <- Bed %>% filter(Chromosome == "chr10") %>% filter(between(Pos,131264000, 131270000))

# Calculate number of calls for each location

Bed2 <- Bed2 %>% mutate(Meth_num = round(Calls_num * Meth_percent /100))

# Shift the position of negative strands by one so that they match the positive strand

Bed2 <- Bed2 %>% mutate(Pos = ifelse(Strand == "-", -1 + Pos, Pos))

# Add the values from positive and negative strands
Bed2 <- Bed2 %>% group_by(Pos) %>% summarise(across(c(Calls_num,  Meth_num), sum))

# Create a new percent methylated column from the new Calls_num and Meth_num values

Bed2 <- Bed2 %>% mutate(Percent_Methylated = Meth_num/Calls_num * 100)

##### Tidy version

Bed2 <- Bed %>% 
  filter(Chromosome == "chr10") %>% 
  filter(between(Pos,131264000, 131270000)) %>%
  mutate(Meth_num = round(Calls_num * Meth_percent /100)) %>% 
  mutate(Pos = ifelse(Strand == "-", -1 + Pos, Pos)) %>% 
  group_by(Pos) %>% 
  summarise(across(c(Calls_num,  Meth_num), sum)) %>%
  mutate(Percent_Methylated = Meth_num/Calls_num * 100)

###### If Megalodon results have been compiled with Modbam2Bed

Bed <- read.csv("/encrypteddata/MGMT/mgmt_retest_01082022/results/t21-240/5hmc/Aggregated_calls.csv", sep = "\t")

# Remove duplicate columns and the ones I don't understand

Bed <- Bed[,c(1:2,6,10,12)]

# Give names to the columns

colnames(Bed) <- c("Chromosome","Pos","Strand","Calls_num","Meth_num")

Bed <- as_tibble(Bed)

Bed3 <- Bed %>% 
  mutate(Pos = ifelse(Strand == "-", -1 + Pos, Pos)) %>% 
  group_by(Pos) %>% 
  summarise(across(c(Calls_num,  Meth_num), sum)) %>%
  mutate(Percent_Methylated = Meth_num/Calls_num * 100)

                       