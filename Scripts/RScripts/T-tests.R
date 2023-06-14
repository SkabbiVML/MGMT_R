# T-tests

# sequence depth of Methylated versus unmethylated in all samples
t.test(Depth~Known_status, data = MGMT_RunSum, var.equal=TRUE)

# sequence depth of Methylated versus unmethylated in Adaptive
Adaptive <- MGMT_RunSum %>% 
  select(Method,Known_status, Depth) %>%
  filter(Method == "AdaptiveSampling")

t.test(Depth~Known_status, data = Adaptive, var.equal=TRUE)

# sequence depth of Methylated versus unmethylated in Barcoded
Barcoded <- MGMT_RunSum %>% 
  select(Method,Known_status, Depth) %>%
  filter(Method == "Barcoded")

t.test(Depth~Known_status, data = Barcoded, var.equal=TRUE)

# sequence depth of Methylated versus unmethylated in Barcoded
Single <- MGMT_RunSum %>% 
  select(Method,Known_status, Depth) %>%
  filter(Method == "Single")

t.test(Depth~Known_status, data = Single, var.equal=TRUE)

# Sequence depth between Barcoded and Single in Radium
Radium <- MGMT_RunSum %>% 
  select(Method,Known_status, Method) %>%
  filter(Series == "Radium")

t.test(Method~Method, data = Radium, var.equal=TRUE)

# Sequence depth between Barcoded and Single in DenStem
DenStem <- MGMT_RunSum %>% 
  select(Method,Known_status, Method) %>%
  filter(Series == "DenStem")

t.test(Method~Method, data = DenStem, var.equal=TRUE)


# sequence depth of concordant versus discordant samples
t.test(Depth~Nano_Pyro_Concordance, data=p, var.equal=T)
