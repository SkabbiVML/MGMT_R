all_results_MGMT <- read.csv("R-markdown/Data/Aggregated_MGMT_24102022.csv")

# Fetch the positions for the two STP27 CpGs
STP27_CpGs <- all_results_MGMT %>% filter(Pos == 131265208 | Pos == 131265574 )

# make a column with the official epic array CpG IDs
STP27_CpGs$CpG_ID <- ifelse(STP27_CpGs$Pos == 131265208, "cg12434587", "cg12981137")

write.csv(STP27_CpGs, "R-markdown/Data/STP27_All_samples.csv")

