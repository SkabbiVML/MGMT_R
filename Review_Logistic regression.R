### Use the Radium samples to make a training model for MGMT-Pyro classification

#make a random split of Radium samples to create training and testing data

dat <- FOUR_CpGs_Group %>% select(Average_Nano,Known_status)

split1<- sample(c(rep(0, 43), rep(1, 19)))
Train <- dat[split1 == 0, ] 
Train <- Train[order(Train$Average_Nano),]
Train$P <- ifelse(Train$Known_status=="Methylated",1,0)

#make test data, either the remaining Radium Samples or all other samples
testRad <- dat[split1 == 1, ] 
testRad$P <- ifelse(testRad$Known_status=="Methylated",1,0)
testRad <- testRad[order(testRad$Average_Nano),]

testOther <- All_MGMT_Pyro %>% 
  group_by(SampleID) %>% 
  summarise(Average_Nano=mean(Methylation_percent), 
            Average_Cov=mean(Valid_cov)) %>%
  inner_join(Samples) %>%
  filter(Average_Cov>2) %>%
  filter(Series != "Radium") %>%
  select(Average_Nano,Known_status)

testOther$P <- ifelse(testOther$Known_status=="Methylated",1,0)
testOther <- testOther[order(testOther$Average_Nano),]

#make the training model

model <- glm(P ~ Average_Nano, family = "binomial", data = Train)

summary(model)

#predict the other samples
predictedRad <- predict(model, testRad, type = "response")
predictedOther <- predict(model, testOther, type = "response")

# plot the training data en fitted curve from model

ggplot(Train, aes(x=Average_Nano, y=P))+
  geom_line(aes(x=Average_Nano, y=model$fitted.values))+
  geom_point(aes(color=Known_status),size=3, alpha=0.7)+
  geom_vline(xintercept = 21)+
  scale_color_brewer(name="Known Status", palette = "Set1")+
  ylab("")+
  xlab("Nanopore sequencing\n(% methylated)")+
  theme_bw(base_size = 16)

#### make the ROC cur
library(pROC)

ROC.Rad <- roc(testRad$P, 
               predictedRad , 
                 legacy.axes = T, 
                 percent = T, 
                 print.auc = T)

ROC.Other <- roc(testOther$P, 
                 predictedOther,
                   legacy.axes = T, 
                   percent = T, 
                   print.auc = T)

rocs_Pyro <- list()

rocs_Pyro[["Retrospective nCATs"]] <- ROC.Rad
rocs_Pyro[["Other"]] <- ROC.Other

##### MAKE ROC curves

#ggroc(ROC.Train) # training set
#ggroc(ROC.Predict) # prediction set

auc(rocs_Pyro$`Retrospective nCATs`)
auc(rocs_Pyro$Other)

ggroc(rocs_Pyro, aes = c("colour","linetype"),size=1.5, alpha=0.8)+
  scale_color_manual("",
                     values = c("black", "darkgrey"),
                     labels=c("Retrospective nCATs\nAUC = 100 %",
                              "\nOther\nAUC = 95.08 %"))+
  scale_linetype_manual("", 
                        values=c("solid","dashed"), 
                        labels=c("Retrospective nCATs\nAUC = 100 %",
                                 "\nOther\nAUC = 95.08 %"))+
  theme_bw(base_size = 16)+
  theme(legend.position = c(.95, .65),
        legend.justification = c("right", "top"),
        legend.margin = margin(6, 6, 6, 6))

#############################
# Do the same for STP-27

###################### STP27 regression model
STP27 <- All_MGMT %>% 
  filter(Pos == 129466944 | Pos == 129467310) %>% 
  select(1:2,4) %>% 
  pivot_wider(names_from = Pos, values_from=Methylation_percent) %>%
  inner_join(MGMT_RunSum) %>%
  filter(AVG_cov>2) %>%
  na.omit()

colnames(STP27) <- c("SampleID", "CpG_31", "CpG_84", "Series", "Method", "Known_status", "Diagnosis","AVG_cov","MED_cov" )

STP27$Series <- ifelse(STP27$Series == "Rapid-CNS", "Rapid-CNS", "Other")
STP27 <- STP27 %>% mutate(STPMean = rowMeans(select(., starts_with("CpG")), na.rm = TRUE))
STP27$P <- ifelse(STP27$Known_status=="Methylated",1,0)

dat2 <- STP27 %>% filter(Series == "Rapid-CNS") 

split1<- sample(c(rep(0, 45), rep(1, 20)))
Stp27.train <- dat2[split1 == 0, ] 


#make test data, either the remaining Radium Samples or all other samples
testRAP <- dat2[split1 == 1, ] 
testOther <- STP27 %>% filter(Series == "Other")


stp.model <- glm(P ~ CpG_31 + CpG_84, family = "binomial", data = Stp27.train)

stp.predicted.RAP <- predict(stp.model, testRAP, type = "response")
stp.predicted.Other <- predict(stp.model, testOther, type = "response")



##### Stp ROC curves

#### make the ROC cur
library(pROC)

stp.ROC.Train <- roc(testRAP$P, 
                     stp.predicted.RAP, 
                     legacy.axes = T, 
                     percent = T, 
                     print.auc = T)

stp.ROC.Predict <- roc(testOther$P, 
                       stp.predicted.Other,
                       legacy.axes = T, 
                       percent = T, 
                       print.auc = T)

auc(stp.ROC.Train)
auc(stp.ROC.Predict)

rocs.STP <- list()

rocs.STP[["Rapid-CNS"]] <- stp.ROC.Train
rocs.STP[["Other"]] <- stp.ROC.Predict

ggroc(rocs.STP, aes = c("colour","linetype"),size=1.5, alpha = 0.7)+
  scale_color_manual("",
                     values = c("black", "darkgrey"),
                     labels=c("Rapid-CNS\nAUC = 97.9 %",
                              "\nOther\nAUC = 93.2 %"))+
  scale_linetype_manual("", 
                        values=c("solid","dashed"), 
                        labels=c("Rapid-CNS\nAUC = 97.9 %",
                                 "\nOther\nAUC = 93.2 %"))+
  theme_bw(base_size = 16)+
  theme(legend.position = c(.95, .65),
        legend.justification = c("right", "top"),
        legend.margin = margin(6, 6, 6, 6))