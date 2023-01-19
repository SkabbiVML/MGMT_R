# T-tests

# sequence depth of Methylated versus unmethylated in all samples
t.test(On_target_seqs~Pyro_Methylation_Status, data = MGMT_RunSum, var.equal=TRUE )

# sequence depth of Methylated versus unmethylated in Adaptive
Adaptive <- MGMT_RunSum %>% 
  select(Comment,Pyro_Methylation_Status, On_target_seqs) %>%
  filter(Comment == "AdaptiveSampling")

t.test(On_target_seqs~Pyro_Methylation_Status, data = Adaptive, var.equal=TRUE)

# sequence depth of Methylated versus unmethylated in Barcoded
Barcoded <- MGMT_RunSum %>% 
  select(Comment,Pyro_Methylation_Status, On_target_seqs) %>%
  filter(Comment == "Barcoded")

t.test(On_target_seqs~Pyro_Methylation_Status, data = Barcoded, var.equal=TRUE)

# sequence depth of Methylated versus unmethylated in Barcoded
Single <- MGMT_RunSum %>% 
  select(Comment,Pyro_Methylation_Status, On_target_seqs) %>%
  filter(Comment == "Single")

t.test(On_target_seqs~Pyro_Methylation_Status, data = Single, var.equal=TRUE)

# Sequence depth between Barcoded and Single in Radium
Radium <- MGMT_RunSum %>% 
  select(Comment,Pyro_Methylation_Status, On_target_seqs) %>%
  filter(Series == "Radium")

t.test(On_target_seqs~Comment, data = Radium, var.equal=TRUE)

# Sequence depth between Barcoded and Single in DenStem
DenStem <- MGMT_RunSum %>% 
  select(Comment,Pyro_Methylation_Status, On_target_seqs) %>%
  filter(Series == "DenStem")

t.test(On_target_seqs~Comment, data = DenStem, var.equal=TRUE)
