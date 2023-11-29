### Use the Radium samples to make a training model 

#make training data
Train <- FOUR_CpGs_Group %>% select(Average_Nano,Known_status)
Train <- Train[order(Train$Average_Nano),]
Train$P <- ifelse(Train$Known_status=="Methylated",1,0)

#make test data
test <- All_MGMT_Pyro %>% 
  group_by(SampleID) %>% 
  summarise(Average_Nano=mean(Methylation_percent), 
            Average_Cov=mean(Valid_cov)) %>%
  inner_join(Samples) %>%
  filter(Average_Cov>2) %>%
  filter(Series != "Radium") %>%
  select(Average_Nano,Known_status)

test$P <- ifelse(test$Known_status=="Methylated",1,0)
test <- test[order(test$Average_Nano),]

#make the training model

model <- glm(P ~ Average_Nano, family = "binomial", data = Train)

#predict the other samples
predicted <- predict(model, test, type = "response")

# plot the training data en fitted curve from model

ggplot(Train, aes(x=Average_Nano, y=P))+
  geom_line(aes(x=Average_Nano, y=model$fitted.values))+
  geom_point(aes(color=Known_status),size=3, alpha=0.7)+
  geom_vline(xintercept = 21.8)+
  scale_color_brewer(name="Known Status", palette = "Set1")+
  ylab("")+
  xlab("Nanopore sequencing\n(% methylated)")+
  theme_bw(base_size = 16)

#### make the ROC cur
library(pROC)

ROC.Train <- roc(Train$P, 
    model$fitted.values, 
    legacy.axes = T, 
    percent = T, 
    print.auc = T)

ROC.Predict <- roc(test$P, 
                   predicted,
                   legacy.axes = T, 
                   percent = T, 
                   print.auc = T)

rocs <- list()

rocs[["Radium"]] <- ROC.Train
rocs[["Other"]] <- ROC.Predict

##### MAKE ROC curves

#ggroc(ROC.Train) # training set
#ggroc(ROC.Predict) # prediction set

auc(rocs$Radium)
auc(rocs$Other)

ggroc(rocs, size=1.5, alpha = 0.7)+
  scale_color_brewer(palette = "Dark2", name="Samples",labels=c("Radium\nAUC = 99.15%","\nOther\nAUC = 95.08"))+
  theme_bw(base_size = 16)+
  theme(legend.position = c(.95, .45),
        legend.justification = c("right", "top"),
        legend.margin = margin(6, 6, 6, 6))

################
###################### STP27 regression model

Stp27.train <- STP27 %>% filter(Series == "Rapid-CNS") %>% select(2,3,6)
Stp27.train$P <- ifelse(Stp27.train$Known_status=="Methylated",1,0)
Stp27.train <- Stp27.train[order(Stp27.train$P),]

stp.model <- glm(P ~ cg12434587 + cg12981137, family = "binomial", data = Stp27.train)

stp.ROC.Train <- roc(Stp27.train$P, 
                 stp.model$fitted.values, 
                 legacy.axes = T, 
                 percent = T, 
                 print.auc = T)
auc(stp.ROC.Train)
