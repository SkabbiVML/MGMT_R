df <- FOUR_CpGs_Group %>% select(Average_Nano,Known_status)

df$Known_status <- ifelse(df$Known_status=="Methylated",1,0)

plot(x=df$Average_Nano , y=df$Known_status)

glm.fit <- glm(df$Known_status ~ df$Average_Nano, family = binomial)
lines(df$Average_Nano,glm.fit$fitted.values)

library(pROC)

par(pty = "s")


roc(df$Known_status, 
    glm.fit$fitted.values, 
    legacy.axes = T, 
    plot = T, 
    percent = T, 
    xlab ="False Positive percentage",
    ylab = "True Positive percentage",
    lwd = 3,
    print.auc = T)

roc.info<- roc(df$Known_status, glm.fit$fitted.values, legacy.axes = T)
roc.df <- data.frame(tpp=roc.info$sensitivities*100, 
                     fpp=(1 - roc.info$specificities)*100,
                     thresholds=roc.info$thresholds)

par(pty = "m")

ggplot(roc.df)

#########3 Sigmoid plot

LigiReg <- FOUR_CpGs_Group %>% 
  select("Average_Nano", "Known_status") 

LigiReg$Status <- ifelse(LigiReg$Known_status == "Methylated", 1,0) 
glm.fit <- glm(LigiReg$SigLine ~ LigiReg$Average_Nano, family = binomial)

LigiReg$FitLine <- glm.fit$fitted.values

ggplot(LigiReg, aes(x=Average_Nano, y=Status))+
  geom_line(aes(x=Average_Nano, y=FitLine))+
  geom_point(aes(color=Known_status),size=3, alpha=0.7)+
  theme_bw()